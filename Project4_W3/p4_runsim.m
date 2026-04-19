%% p4_week3_runsim.m
% Project 4 Week 3 run script
% Runs the P4 model for the full 60-minute event and checks deliverable constraints

clc; clear; clf;

% Initialize model/workspace
p4_init;

fprintf('\nRunning p4_model.slx for Week 3...\n');
simout = sim('p4_model.slx', 'StopTime', num2str(sim_stop_time));
fprintf('Simulation complete.\n\n');

% -----------------------
% Read outputs from sim
% -----------------------
t = simout.tout;

X_ts = simout.X;
Y_ts = simout.Y;
if isa(X_ts,'timeseries')
    X = X_ts.Data;
    Y = Y_ts.Data;
else
    X = X_ts;
    Y = Y_ts;
end

vx_ts  = simout.vx_out;
soc_ts = simout.soc_out;
if isa(vx_ts,'timeseries')
    vx  = vx_ts.Data;
    soc = soc_ts.Data;
else
    vx  = vx_ts;
    soc = soc_ts;
end

X   = X(:);
Y   = Y(:);
vx  = vx(:);
soc = soc(:);
t   = t(:);

% -----------------------
% Metrics / constraints
% -----------------------
results = evaluate_week3_run(X, Y, vx, soc, path);

fprintf('========== WEEK 3 SUMMARY ==========\n');
fprintf('Target speed                : %.1f mph\n', carData.vxd_mph);
fprintf('Lookahead distance          : %.1f m\n', path.pure_pursuit_lookaheaddist);
fprintf('Simulation time             : %.1f s\n', t(end));
fprintf('Laps completed              : %.2f\n', results.laps);
fprintf('Average speed               : %.2f mph\n', mean(vx)*2.23694);
fprintf('Peak speed                  : %.2f mph\n', max(vx)*2.23694);
fprintf('Initial SOC                 : %.2f %%\n', soc(1)*100);
fprintf('Final SOC                   : %.2f %%\n', soc(end)*100);
fprintf('Minimum SOC                 : %.2f %%\n', min(soc)*100);
fprintf('Maximum SOC                 : %.2f %%\n', max(soc)*100);
fprintf('Max cross-track error       : %.2f m\n', results.max_cte);
fprintf('Track half-width            : %.2f m\n', path.width/2);
fprintf('Stayed on track?            : %s\n', tf2str(results.on_track));
fprintf('SOC stayed within 10-95%%?   : %s\n', tf2str(results.soc_ok));
fprintf('Overall Week 3 valid?       : %s\n', tf2str(results.week3_valid));
fprintf('====================================\n');

% -----------------------
% Plots
% -----------------------
figure(1); clf;
plot(path.xpath, path.ypath, 'k-', 'LineWidth', 1.5); hold on;
plot(path.xinpath, path.yinpath, 'b--', 'LineWidth', 1.0);
plot(path.xoutpath, path.youtpath, 'r--', 'LineWidth', 1.0);
plot(X, Y, 'm-', 'LineWidth', 1.5);
legend('Track Centre','Inner Edge','Outer Edge','Vehicle Path','Location','best');
xlabel('X (m)'); ylabel('Y (m)');
title('Project 4 Week 3 - Full Vehicle Trajectory');
axis equal; grid on;

figure(2); clf;
plot(t, vx*2.23694, 'b-', 'LineWidth', 1.5); hold on;
yline(carData.vxd_mph, 'r--', 'Desired Speed', 'LineWidth', 1.2);
xlabel('Time (s)'); ylabel('Speed (mph)');
title('Project 4 Week 3 - Vehicle Speed vs Time');
grid on; legend('Actual Speed','Desired Speed','Location','best');

figure(3); clf;
plot(t, soc*100, 'g-', 'LineWidth', 1.5); hold on;
yline(10, 'r--', 'Min 10%', 'LineWidth', 1.2);
yline(95, 'b--', 'Max 95%', 'LineWidth', 1.2);
xlabel('Time (s)'); ylabel('State of Charge (%)');
title('Project 4 Week 3 - Battery SOC vs Time');
ylim([0 100]); grid on;

figure(4); clf;
plot(t, results.cte, 'k-', 'LineWidth', 1.2); hold on;
yline(path.width/2, 'r--', 'Track Limit', 'LineWidth', 1.2);
yline(-path.width/2, 'r--', 'Track Limit', 'LineWidth', 1.2);
xlabel('Time (s)'); ylabel('Cross-track error (m)');
title('Project 4 Week 3 - Cross-Track Error');
grid on;

save('week3_results.mat', 't', 'X', 'Y', 'vx', 'soc', 'results');

% =============================================================
function results = evaluate_week3_run(X, Y, vx, soc, path)
    [laps, arc, L, cte] = compute_laps_and_error(X, Y, path);

    results.laps        = laps;
    results.arc         = arc;
    results.track_length = L;
    results.cte         = cte;
    results.max_cte     = max(abs(cte));
    results.on_track    = all(abs(cte) <= path.width/2);
    results.soc_ok      = all(soc >= 0.10) && all(soc <= 0.95);
    results.week3_valid = results.on_track && results.soc_ok;
end

function [laps, arc, L, cte] = compute_laps_and_error(X, Y, path)
    cx = path.xpath(:);
    cy = path.ypath(:);
    nv = length(X);

    ds   = sqrt(diff(cx).^2 + diff(cy).^2);
    cumS = [0; cumsum(ds)];
    L    = cumS(end);

    arc = zeros(nv,1);
    cte = zeros(nv,1);

    for i = 1:nv
        d2 = (cx - X(i)).^2 + (cy - Y(i)).^2;
        [d2min, idx] = min(d2);
        arc(i) = cumS(idx);
        cte(i) = sqrt(d2min);
    end

    delta_arc = diff(arc);
    delta_arc(delta_arc < -L/2) = delta_arc(delta_arc < -L/2) + L;
    delta_arc(delta_arc >  L/2) = delta_arc(delta_arc >  L/2) - L;

    total_forward_arc = sum(max(delta_arc, 0));
    laps = total_forward_arc / L;
end

function out = tf2str(flag)
    if flag
        out = 'YES';
    else
        out = 'NO';
    end
end