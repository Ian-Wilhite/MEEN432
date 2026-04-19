%% This code generates the track parameters in the path structur

% FOR PURE PURSUIT INITIALIZATION
path.pure_pursuit_lookaheaddist = 10; % m - look ahead distance

path.radius = 200; % assigning the radius of the track into a structure called "path"
path.l_st = 900; % assigning the straight length of the track into a structure called "path"

path.width = 15; % assigning the width of the track into a structure called "track"

radius = path.radius; % assigning radius from struct
    
l_st = path.l_st; % assigning straightaway from struct

l_curve = pi * radius; % assigning the length of curve (circumference divided by two since its only half)

path.total_length = 2 * l_st + 2 * l_curve; % the total length of the track

delta_s = 10;
npts = round(path.total_length/delta_s); % number of data points

delta_s = path.total_length/(npts-1); % the step size for the total length of track for 1 - number of points
 
delta_theta = delta_s / radius; % the step size for the curve portions

path.xpath = zeros(npts,1); % an array of zeros for the x component center path
path.ypath = zeros(npts,1); % an array of zeros for the y component center path
path.tpath = zeros(npts,1); % an array of zeros for the theta component of the overall path
path.xinpath = zeros(npts,1); % an array of zeros for the x component inside path
path.yinpath = zeros(npts,1); % an array of zeros for the y component inside path
path.xoutpath = zeros(npts,1); % an array of zeros for the x component outside path
path.youtpath = zeros(npts,1); % an array of zeros for the y component outside path



i = 1; % counter for what data point we're at

while i < npts % while the counter is less than the number of points
path.yinpath(1) = 7.5; % initializing where inside border starts 
path.youtpath(1) = -7.5; % initializing where outside border starts

    if path.xpath(i) < l_st % if the x component at the index of the counter is less than the straightaway length

        if path.xpath(i) >= 0 % if that x component is more than or equal to zero

            if path.ypath(i) < radius % if the y component at the counter index is less than the radius length

                % at this condition, we are at the start point of the
                % track, so the first x and y will be at 0, y will continue
                % to be 0 since we are at the first straightaway section

                path.xpath(i+1) = path.xpath(i) + delta_s; % at the next index, the x component is the previous plus the step size 
                path.ypath(i+1) = path.ypath(i); % since the y component has not reached the curve yet, it will stay the same point at the straight away

                path.xinpath(i+1) = path.xinpath(i) + delta_s; 
                path.yinpath(i+1) = path.yinpath(i);

                path.xoutpath(i+1) = path.xoutpath(i) + delta_s;
                path.youtpath(i+1) = path.youtpath(i);

                path.tpath(i+1) = 0; % straightaway is at theta of 0, think of it like a unit circle
            else

                % this is for the second straightaway section
                path.xpath(i+1) = path.xpath(i) - delta_s; % since second straightaway is going in the negative x direction, subtract by delta_s
                path.ypath(i+1) = path.ypath(i); % y component will stay the same a the second straightaway (400 meters in y axis)

                path.xinpath(i+1) = path.xinpath(i) - delta_s; 
                path.yinpath(i+1) = path.yinpath(i);

                path.xoutpath(i+1) = path.xoutpath(i) - delta_s; 
                path.youtpath(i+1) = path.youtpath(i);
                path.tpath(i+1) = pi; % second straightaway is at theta of pi, remember unit circle
            end
        else

            % for the second curved section
            cx = 0; cy = radius;
            rx = path.xpath(i) - cx; ry = path.ypath(i) - cy;
            rxin = path.xinpath(i) - cx; ryin = path.yinpath(i) - cy;
            rxout = path.xoutpath(i) - cx; ryout = path.youtpath(i) - cy;

            tt = rotate([rx;ry], delta_theta); 
            ttin = rotate([rxin;ryin], delta_theta); 
            ttout = rotate([rxout;ryout], delta_theta); 

            path.xpath(i+1) = tt(1) + cx; 
            path.ypath(i+1) = tt(2) + cy;

            path.xinpath(i+1) =  ttin(1) + cx; 
            path.yinpath(i+1) = ttin(2) + cy;

            path.xoutpath(i+1) = ttout(1) + cx; 
            path.youtpath(i+1) = ttout(2) + cy;
            path.tpath(i+1) = path.tpath(i) + delta_theta;
        end
    else

        % once the x component of the path is larger than the length of the
        % straightaway, it has reached the first curved portion of the track

        cx = l_st; cy = radius;                 % assigning the length of straightaway and the radius as x and y components respectively
        rx = path.xpath(i) - cx; ry = path.ypath(i) - cy; % subtracting the x and y path components by their original lengths
        rxin = path.xinpath(i) - cx; ryin = path.yinpath(i) - cy;
        rxout = path.xoutpath(i) - cx; ryout = path.youtpath(i) - cy;

        tt = rotate([rx;ry],delta_theta);    % rotating the rx and ry components by delta theta using p2_rotate function
        ttin = rotate([rxin;ryin], delta_theta); 
        ttout = rotate([rxout;ryout], delta_theta);

        path.xpath(i+1) = tt(1) + cx;                % adding the transformed x component to the original length of straightaway
        path.ypath(i+1) = tt(2) + cy;                % adding the transformed y component to the original radius 

        path.xinpath(i+1) =  ttin(1) + cx; 
        path.yinpath(i+1) = ttin(2) + cy;

        path.xoutpath(i+1) = ttout(1) + cx;
        path.youtpath(i+1) = ttout(2) + cy;
        path.tpath(i+1) = path.tpath(i) + delta_theta; 
    end
    i = i + 1;    
end
% === Elevation Profile ===
% Height(x) = exp((x-450)/50) / (1 + exp((x-450)/50)) * 10  [m]
% Applied on straight segments only; curved sections assumed flat.
% Signed grade: positive = uphill in the direction of travel.
path.hpath       = zeros(npts, 1);  % height above datum (m)
path.grade       = zeros(npts, 1);  % dh/ds along direction of travel (dimensionless)
path.theta_grade = zeros(npts, 1);  % grade angle (rad), + = uphill

for i = 1:npts
    x   = path.xpath(i);
    psi = path.tpath(i);

    on_start  = abs(psi)       < 0.01;   % heading ~0   → moving in +x (uphill)
    on_return = abs(psi - pi)  < 0.01;   % heading ~pi  → moving in -x (downhill)

    if on_start || on_return
        sig  = exp((x - 450)/50) / (1 + exp((x - 450)/50));
        dhdx = (10/50) * sig * (1 - sig);   % d(Height)/dx

        path.hpath(i) = sig * 10;

        if on_start
            path.grade(i) =  dhdx;   % vehicle climbs as x increases
        else
            path.grade(i) = -dhdx;   % vehicle descends as x decreases
        end
        path.theta_grade(i) = atan(path.grade(i));
    end
end

% Cumulative arc-length along path (used by Simulink lookup table)
ds_vec   = sqrt(diff(path.xpath).^2 + diff(path.ypath).^2);
path.cumS = [0; cumsum(ds_vec)];

% Make track struct available to Simulink blocks and post-processing
assignin('base','path',path);

% Also export a 'track' struct with the field names expected by the
% provided raceStat.m (track.radius, track.width, track.l_straightaways)
track.radius         = path.radius;
track.width          = path.width;
track.l_straightaways = path.l_st;
assignin('base','track',track);
function v_rot = rotate(v, theta)
    % Rotate 2D vector v by angle theta (radians)
    R = [cos(theta), -sin(theta);
         sin(theta),  cos(theta)];
    v_rot = R * v;
end
