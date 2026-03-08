% Animate the vehicle path produced by the Simulink model
%
% Week 3 note: Outputs must come from "To Workspace" blocks.
% This script supports either:
%   - out.X, out.Y, out.psi (Timeseries or Structure with time)
%   - SimulationOutput signals (fallback)

try
    simout = sim("Project_2_Kinematic_Model.slx");
catch simErr
    error("Simulink run failed: %s", simErr.message);
end

% Prefer To Workspace outputs (required by assignment)
[car_X, car_time] = getWorkspaceSignal('out','X');
[car_Y, ~]        = getWorkspaceSignal('out','Y');
[car_psi, ~]      = getWorkspaceSignal('out','psi');

% Fallback: try SimulationOutput if To Workspace variables not found
if isempty(car_X) || isempty(car_Y)
    car_X = getSimSignal(simout,"X");
    car_Y = getSimSignal(simout,"Y");
    car_psi = getSimSignal(simout,"psi");
    if isprop(simout,'tout')
        car_time = simout.tout;
    elseif isprop(simout,'tout')
        car_time = simout.tout;
    else
        car_time = [];
    end
end

if isempty(car_X) || isempty(car_Y)
    error(["Missing simulation outputs. Add To Workspace blocks for out.X and out.Y ", ...
           "(Timeseries recommended)."]);
end

if isempty(car_psi)
    % Heading is required by the project statement, but raceStat doesn't use it.
    % If you haven't exported it yet, set it to zeros so animation still runs.
    car_psi = zeros(size(car_X));
end

if isempty(car_time)
    % As a last resort, create a sample index as time.
    car_time = (0:numel(car_X)-1)';
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

function [sig, t] = getWorkspaceSignal(containerVar, fieldName)
% Get signal + time from a To Workspace output.
% Supports Timeseries and Structure with time.
    sig = [];
    t   = [];
    if ~evalin('base', sprintf("exist('%s','var')", containerVar))
        return;
    end
    outv = evalin('base', containerVar);
    if ~isstruct(outv) || ~isfield(outv, fieldName)
        return;
    end
    v = outv.(fieldName);
    if isa(v,'timeseries')
        sig = v.Data;
        t   = v.Time;
        return;
    end
    if isstruct(v) && isfield(v,'signals') && isfield(v,'time')
        sig = v.signals.values;
        t   = v.time;
        return;
    end
    if isnumeric(v)
        sig = v(:);
        return;
    end
end
