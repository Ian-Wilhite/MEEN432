%% p4_runsim.m
% Project 4 Week 2 run script
% Runs the P4 model for the full 60-minute event and checks deliverable constraints.

clc; clear; clf;

% Initialize model/workspace
p4_init;

fprintf('\nRunning p4_model.slx for Week 2...\n');
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
results = evaluate_week2_run(X, Y, vx, soc, path);

fprintf('========== WEEK 2 SUMMARY ==========\n');
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
fprintf('Reached at least 5 laps?    : %s\n', tf2str(results.five_laps_ok));
fprintf('Overall Week 2 valid?       : %s\n', tf2str(results.week2_valid));
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
title('Project 4 Week 2 - Full Vehicle Trajectory');
axis equal; grid on;

figure(2); clf;
plot(t, vx*2.23694, 'b-', 'LineWidth', 1.5); hold on;
yline(carData.vxd_mph, 'r--', 'Desired Speed', 'LineWidth', 1.2);
xlabel('Time (s)'); ylabel('Speed (mph)');
title('Project 4 Week 2 - Vehicle Speed vs Time');
grid on; legend('Actual Speed','Desired Speed','Location','best');

figure(3); clf;
plot(t, soc*100, 'g-', 'LineWidth', 1.5); hold on;
yline(10, 'r--', 'Min 10%', 'LineWidth', 1.2);
yline(95, 'b--', 'Max 95%', 'LineWidth', 1.2);
xlabel('Time (s)'); ylabel('State of Charge (%)');
title('Project 4 Week 2 - Battery SOC vs Time');
ylim([0 100]); grid on;

figure(4); clf;
plot(t, results.cte, 'k-', 'LineWidth', 1.2); hold on;
yline(path.width/2, 'r--', 'Track Limit', 'LineWidth', 1.2);
yline(-path.width/2, 'r--', 'Track Limit', 'LineWidth', 1.2);
xlabel('Time (s)'); ylabel('Cross-track error (m)');
title('Project 4 Week 2 - Cross-Track Error');
grid on;

% --- Elevation profile figure ---
figure(5); clf;
yyaxis left
plot(path.cumS, path.hpath, 'b-', 'LineWidth', 1.5);
ylabel('Height (m)'); ylim([-1 12]);
yyaxis right
plot(path.cumS, path.theta_grade * 180/pi, 'r-', 'LineWidth', 1.2);
ylabel('Grade angle (deg)');
xlabel('Arc-length position (m)');
title('Track Elevation and Grade Profile');
xline(path.l_st,                    'k--', 'End Straight 1', 'LabelVerticalAlignment','bottom');
xline(path.l_st + pi*path.radius,   'k--', 'End Curve 1',    'LabelVerticalAlignment','bottom');
xline(2*path.l_st + pi*path.radius, 'k--', 'End Straight 2', 'LabelVerticalAlignment','bottom');
grid on; legend('Height','Grade angle','Location','best');

% Save summary — use elevated_results.mat once grade is wired into Simulink
if exist('grade_lut', 'var') && any(path.theta_grade ~= 0)
    save_file = 'elevated_results.mat';
else
    save_file = 'week2_results.mat';
end
save(save_file, 't', 'X', 'Y', 'vx', 'soc', 'results');
fprintf('Results saved to %s\n', save_file);

% =============================================================
function results = evaluate_week2_run(X, Y, vx, soc, path)
    [laps, arc, L, cte] = compute_laps_and_error(X, Y, path);

    results.laps         = laps;
    results.arc          = arc;
    results.track_length = L;
    results.cte          = cte;
    results.max_cte      = max(abs(cte));
    results.on_track     = all(abs(cte) <= path.width/2);
    results.soc_ok       = all(soc >= 0.10) && all(soc <= 0.95);
    results.five_laps_ok = laps >= 5.0;
    results.week2_valid  = results.on_track && results.soc_ok && results.five_laps_ok;
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