function stats = raceStat(car_X, car_Y, car_time, path)
%RACESTAT Basic lap statistics and track-out checks
%   stats = raceStat(X,Y,t,path)
%     X, Y : car coordinates
%     t    : timestamps matching X, Y
%     path : struct with width, l_st, radius

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
