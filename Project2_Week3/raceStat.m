function stats = raceStat(car_X, car_Y, car_time, path)
%RACESTAT Basic lap statistics and track-out checks
%   stats = raceStat(X,Y,t,path)

    % Allow calling raceStat with no arguments:
    if nargin < 4
        [car_X, car_Y, car_time, path] = raceStat_fromBase();
    end

    % Allow calling raceStat with no inputs (pulls X,Y,t,path from base workspace)
    if nargin == 0
        [car_X, car_Y, car_time, path] = raceStat_fromBase();
    elseif nargin < 4
        error('raceStat requires inputs: stats = raceStat(X, Y, t, path)');
    end

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
    
function [car_X, car_Y, car_time, path] = raceStat_fromBase()
    % Prefer variables that already exist from your sim/scripts
    if evalin('base','exist("car_X","var") && exist("car_Y","var") && exist("car_time","var")')
        car_X    = evalin('base','car_X');
        car_Y    = evalin('base','car_Y');
        car_time = evalin('base','car_time');
    else
        % Fall back to To Workspace style: out.X / out.Y
        if ~evalin('base','exist("out","var")')
            error(['Could not find car_X/car_Y/car_time or out.X/out.Y in base workspace. ' ...
                   'Run the simulation first (run.m), or add To Workspace blocks exporting out.X and out.Y.']);
        end
        out = evalin('base','out');
        if ~isfield(out,'X') || ~isfield(out,'Y')
            error('Found variable "out" but missing out.X or out.Y. Check To Workspace block settings.');
        end

        Xts = out.X; Yts = out.Y;
        if isstruct(Xts) && isfield(Xts,'signals') && isfield(Xts,'time')
            car_X = Xts.signals.values;
            car_Y = Yts.signals.values;
            car_time = Xts.time;
        else
            car_X = Xts.Data;
            car_Y = Yts.Data;
            car_time = Xts.Time;
        end
    end

    % Path/track dims
    if evalin('base','exist("path","var")')
        path = evalin('base','path');
    elseif evalin('base','exist("track","var")')
        path = evalin('base','track');  % allow either name
    else
        error('Could not find "path" (or "track") in base workspace.');
    end
end

    % 1. Total Distance Traveled
    dx = diff(car_X);
    dy = diff(car_Y);
    total_distance = sum(sqrt(dx.^2 + dy.^2));

    % 2. Total Time and Average Speed
    total_time = car_time(end) - car_time(1);
    avg_speed = total_distance / total_time;

    % 5. Package into output struct
    stats = struct();
    stats.total_distance = total_distance;
    stats.total_time = total_time;
    stats.avg_speed = avg_speed;
    stats.loops = loops;
    stats.loop_times = tloops(1:j);
    stats.offtrack = struct('X', Xerr(1:k), 'Y', Yerr(1:k), 't', terr(1:k));

    % 6. Display basic info
    fprintf('--- Race Statistics ---\n');
    fprintf('Total Distance Traveled: %.2f m\n', total_distance);
    fprintf('Total Time: %.2f s\n', total_time);
    fprintf('Average Speed: %.2f m/s\n', avg_speed);
    fprintf('Loops made: %d\n', loops);
    for n = 2:length(stats.loop_times)
        fprintf('Lap %d time: %.2f s\n', n-1, stats.loop_times(n) - stats.loop_times(n-1));
    end
    if k > 0
        fprintf('Car left the track %d times.\n', k);
    else
        fprintf('Car did not leave the track.\n');
    end
end

function yesorno = insideTrack(x,y,section,path)
    switch section
        case 1
            yesorno = ((y < (0.0 + path.width)) && (y > (0.0 - path.width)));
        case {2, 3}
            rad = sqrt((x - path.l_st)^2 + (y - path.radius)^2);
            yesorno = ((rad < path.radius + path.width) && ...
                    (rad > path.radius - path.width));
        case 4
            yesorno = ((y < (2 * path.radius + path.width)) && ...
                    (y > (2 * path.radius - path.width)));
        case {5, 6}
            rad = sqrt((x - 0.0)^2 + (y - path.radius)^2);
            yesorno = ((rad < path.radius + path.width) && ...
                    (rad > path.radius - path.width));
        otherwise
            yesorno = false;
    end
end

function [car_X, car_Y, car_time, path] = raceStat_fromBase()
% Pull required inputs from base workspace (supports out.X timeseries or X vectors)

    % Track/path
    if evalin('base',"exist('path','var')")
        path = evalin('base',"path");
    elseif evalin('base',"exist('track','var')")
        tr = evalin('base',"track");
        path = struct();
        path.radius = tr.radius;
        path.width  = tr.width;
        if isfield(tr,'l_straightaways')
            path.l_st = tr.l_straightaways;
        elseif isfield(tr,'l_st')
            path.l_st = tr.l_st;
        else
            error("track struct found but missing straightaway length field.");
        end
    else
        error("No track struct found in base workspace. Run gentrack.m first.");
    end

    % Signals
    if evalin('base',"exist('out','var')")
        out = evalin('base',"out");
        if isstruct(out) && isfield(out,'X') && isfield(out,'Y')
            [car_X, car_time] = unpack(out.X);
            [car_Y, ~]        = unpack(out.Y);
            return;
        end
    end

    if evalin('base',"exist('X','var')") && evalin('base',"exist('Y','var')")
        [car_X, car_time] = unpack(evalin('base',"X"));
        [car_Y, ~]        = unpack(evalin('base',"Y"));
        if isempty(car_time) && evalin('base',"exist('t','var')")
            car_time = evalin('base',"t");
        end
        return;
    end

    error("Could not find X and Y in base workspace. Ensure To Workspace blocks export out.X and out.Y (Timeseries).");

    function [sig, t] = unpack(v)
        sig = [];
        t = [];
        if isa(v,'timeseries')
            sig = v.Data; t = v.Time; return;
        end
        if isstruct(v) && isfield(v,'signals') && isfield(v,'time')
            sig = v.signals.values; t = v.time; return;
        end
        if isnumeric(v)
            sig = v;
        end
    end
end
