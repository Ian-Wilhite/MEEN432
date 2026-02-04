% Animate the vehicle path produced by the Simulink model
try
    simout = sim("Project_2_Kinematic_Model.slx");
catch simErr
    error("Simulink run failed: %s", simErr.message);
end

car_X = getSimSignal(simout,"X");
car_Y = getSimSignal(simout,"Y");
car_psi = getSimSignal(simout,"psi");
car_time = simout.tout;

if isempty(car_X) || isempty(car_Y) || isempty(car_psi)
    error("Missing simulation outputs: expected signals X, Y, and psi in simout.");
end

fh = figure();
fh.WindowState = 'maximized';
hold on
plot(path.xpath,path.ypath,'--r'); axis equal; % center line
plot(path.xinpath, path.yinpath, 'b'); axis equal; % inside border 
plot(path.xoutpath, path.youtpath,'b'); axis equal; % outside border
axis([min(path.xoutpath)-100 , max(path.xoutpath)+100 , min(path.youtpath)-100 , max(path.youtpath)+100])
xlabel('X Distance (m)')
ylabel('Y Distance (m)')
title('Project 2 Track')
grid
h = animatedline('Color',[0 0.4 0]); % trail

L = 15;
width = 5;
for i = 1:length(car_X)
    x = car_X(i);
    y = car_Y(i);
    psi = car_psi(i);

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
    addpoints(h,x,y);
    a = polyshape(rcar);
    ap = plot(a, 'FaceColor', 'k');
    drawnow limitrate;
    pause(0.02);
    if i ~= length(car_X)
        delete(ap);
    end
end

race = raceStat(car_X, car_Y, car_time, path); %#ok<NASGU>
hold off;

function sig = getSimSignal(simout,name)
    if isprop(simout,name)
        data = simout.(name);
        if isprop(data,"Data")
            sig = data.Data;
        else
            sig = data;
        end
    else
        sig = [];
    end
end
