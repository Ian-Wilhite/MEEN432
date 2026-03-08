% From a Kinematic Model
simout = sim("Project_2_Kinematic_Model.slx");
car_X = simout.X.Data;
car_Y = simout.Y.Data;
car_psi = simout.psi.Data;
car_time = simout.tout;

path_X = path.xpath;
path_Y = path.ypath;
% Animating the patch
fh = figure();
fh.WindowState = 'maximized';
hold on
plot(path.xpath,path.ypath,'--r'); axis equal; % plots center line
plot(path.xinpath, path.yinpath, 'b'); axis equal; % plots inside border 
plot(path.xoutpath, path.youtpath,'b'); axis equal; % plots outside border
axis([min(path.xoutpath)-100 , max(path.xoutpath)+100 , min(path.youtpath)-100 , max(path.youtpath)+100]) % limits axis for better viewing
xlabel('X Distance (m)') % x axis label
ylabel('Y Distance (m)') % y axis label
title('Project 2 Track') % title
grid
h = animatedline; % handle of animatedline

% 
L = 15;
width = 5;
for i = 1:length(car_X)
    x = car_X(i);
    y = car_Y(i);
    psi = car_psi(i); % heading angle (radians)

    % Car vertices (centered at origin)
    car = [-L/2 -width/2;
           -L/2  width/2;
            L/2  width/2;
            L/2 -width/2];

    % Rotation matrix
    R = [cos(psi) -sin(psi);
         sin(psi)  cos(psi)];

    % Rotate and translate
    rcar = (R * car')';
    rcar = rcar + [x, y];

    % Plot
%    addpoints(h,x,y);
    a = polyshape(rcar);
    ap = plot(a, 'FaceColor', 'k');
    drawnow limitrate;
    pause(0.05);
    if i ~= length(car_X)
        delete(ap);
    end
end
race = raceStat(car_X, car_Y, car_time, path);
hold off;
function stats = raceStat(car_X, car_Y, car_time, path)
    % Compute race statistics and check path adherence

    %========================================================
    % 
    % Usage: rs = raceStat(X,Y,t,path)
    %
    % Inputs:   X, Y are coordinates from your vehicle simulations. 
    %           t is the set times corresponding to X and Y
    %           path is a structure of with fields "width" (width of the 
    %                  track), "l_st" (length of the straight away), and 
    %                  "radius" (radius of the curved section)
    %
    % Outputs: raceStats is a structure with the following fields:
    %
    %   loops - this is an integer that tells how many loops around
    %              the track the vehicle has gone around
    %   tloops - this is an array that tells you the time(s) that the
    %              start line was crossed 
    %   lefftTrack - this has the X,Y and t values when the vehicle
    %              went outside the track
    %   
    %========================================================
    
    prev_section = 6;
    loops = -1;
    j = 0;
    k = 0;
    Xerr = zeros(1,100000);
    Yerr = zeros(1,100000);
    terr = zeros(1,100000);
    tloops = zeros(1,100000);
    for i = 1:length(car_X)
        if car_X(i) < path.l_st
            if car_X(i) >= 0
                if car_Y(i) < path.radius
                    section = 1;
                else
                    section = 4;
                end
            else
                if car_Y(i) < path.radius
                    section = 6;
                else
                    section = 5;
                end
            end
        else
            if car_Y(i) < path.radius
                section = 2;
            else
                section = 3;
            end
        end
        if ((prev_section == 6) && (section == 1))
            loops = loops  + 1;
            j = j+1;
            tloops(j) = car_time(i);
        end
        prev_section = section;
        if ~insideTrack(car_X(i),car_Y(i),section,path)
            k = k+1;
            Xerr(k) = car_X(i);
            Yerr(k) = car_Y(i);
            terr(k) = car_time(i);
        end
    end

    % 1. Total Distance Traveled
    dx = diff(car_X);
    dy = diff(car_Y);
    total_distance = sum(sqrt(dx.^2 + dy.^2));

    % 2. Total Time and Average Speed
    total_time = car_time(end) - car_time(1);
    avg_speed = total_distance / total_time;


    %for i = 1:length(car_X)
        % Find the closest centerline point to the car
        %dx = path.xpath - car_X(i);
        %dy = path.ypath - car_Y(i);
        %dist = sqrt(dx.^2 + dy.^2);
        %[min_dist, ~] = min(dist);
    %end

    % 5. Package into output struct
    stats = struct();
    stats.total_distance = total_distance;
    stats.total_time = total_time;
    stats.avg_speed = avg_speed;


    % 6. Display basic info
    fprintf('--- Race Statistics ---\n');
    fprintf('Total Distance Traveled: %.2f m\n', total_distance);
    fprintf('Total Time: %.2f s\n', total_time);
    fprintf('Average Speed: %.2f m/s\n', avg_speed);
    fprintf('Loops made: %\n', loops);
    %for i = 2:length(tloops)
        %fprintf('Lap %d time: %.2f s\n', i-1, tloops(i));
    %end
    if k > 0
        fprintf('Car left the track.\n');
    else
        fprintf('Car did not leave the track.\n');
    end

end
function yesorno = insideTrack(x,y,section,path)
    switch section
        case 1
            if ((y < (0.0 + path.width)) && (y > (0.0 - path.width))) 
                yesorno = 1;
            else
                yesorno = 0;
            end
        case {2, 3}
            rad = sqrt((x - path.l_st)^2 + (y - path.radius)^2);
            if ((rad < path.radius + path.width) && ...
                    (rad > path.radius - path.width))
                yesorno = 1;
            else
                yesorno = 0;
            end
        case 4
            if ((y < (2 * path.radius + path.width)) && ...
                    (y > (2 * path.radius - path.width))) 
                yesorno = 1;
            else
                yesorno = 0;
            end        
        case {5, 6}
            rad = sqrt((x - 0.0)^2 + (y - path.radius)^2);
            if ((rad < path.radius + path.width) && ...
                    (rad > path.radius - path.width))
                yesorno = 1;
            else
                yesorno = 0;
            end
        otherwise
            print("error");
    end
end
