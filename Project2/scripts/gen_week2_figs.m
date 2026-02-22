% gen_week2_figs.m  — headless Week 2 figure generation
% Run from Project2/ directory:
%   run('scripts/gen_week2_figs.m')

close all; clc;
addpath(pwd);

init;      % vehicle parameters → carData, carDataBus
gentrack;  % track geometry   → path

figdir = fullfile(pwd, 'figures');
if ~exist(figdir, 'dir'); mkdir(figdir); end

%% Run Simulink model (headless)
fprintf('Running Simulink model ...\n');
simout = sim('Project_2_Kinematic_Model.slx');

% Extract signals
function sig = getsig(simout, name)
    if isprop(simout, name)
        d = simout.(name);
        if isprop(d, 'Data'); sig = d.Data; else; sig = d; end
    else
        sig = [];
    end
end

car_X   = getsig(simout, 'X');
car_Y   = getsig(simout, 'Y');
car_psi = getsig(simout, 'psi');
t       = simout.tout;

if isempty(car_X) || isempty(car_Y)
    error('Simulation did not output X or Y. Check To Workspace block names.');
end

%% ---- Figure 1: Track + simulated vehicle path -------------------------
f1 = figure('Visible','off','Position',[100 100 900 500]);
hold on; axis equal; grid on;
plot(path.xoutpath, path.youtpath, 'b-',  'LineWidth', 1.2, 'DisplayName','Track boundary');
plot(path.xinpath,  path.yinpath,  'b-',  'LineWidth', 1.2, 'HandleVisibility','off');
plot(path.xpath,    path.ypath,    '--r', 'LineWidth', 0.8, 'DisplayName','Center line');
plot(car_X, car_Y, 'g-', 'LineWidth', 1.8, 'DisplayName','Vehicle path');
plot(car_X(1), car_Y(1), 'go', 'MarkerSize', 8, 'MarkerFaceColor','g', 'DisplayName','Start');

% Draw vehicle patch at a few waypoints along the simulated path
L = 15; W = 5;
n = length(car_X);
poses = unique(round(linspace(1, n, 6)));
for k = poses
    x   = car_X(k);
    y   = car_Y(k);
    psi = car_psi(k);
    corners = [-L/2 -W/2; -L/2 W/2; L/2 W/2; L/2 -W/2];
    R = [cos(psi) -sin(psi); sin(psi) cos(psi)];
    rc = (R * corners')' + [x y];
    fill(rc(:,1), rc(:,2), [0.2 0.2 0.8], 'FaceAlpha', 0.4, 'EdgeColor','k', 'HandleVisibility','off');
end

xlabel('X (m)'); ylabel('Y (m)');
title(sprintf('Vehicle Path — v_x = %.0f m/s (%.0f km/h)', carData.vxd, carData.vxd*3.6));
legend('Location','southoutside','Orientation','horizontal','NumColumns',3);
saveas(f1, fullfile(figdir, 'vehicle_path.png'));
fprintf('Saved vehicle_path.png\n');

%% ---- Figure 2: Heading and steering over time -------------------------
f2 = figure('Visible','off','Position',[100 100 900 400]);

subplot(2,1,1);
plot(t, car_psi * 180/pi, 'b-', 'LineWidth', 1.2);
ylabel('\psi (deg)');
title('Vehicle Heading and Steering Angle vs. Time');
grid on;

% Try to pull steering angle if available
delta = getsig(simout, 'delta_f');
if ~isempty(delta) && length(delta) == length(t)
    subplot(2,1,2);
    plot(t, delta * 180/pi, 'r-', 'LineWidth', 1.2);
    ylabel('\delta_f (deg)');
else
    subplot(2,1,2);
    text(0.5, 0.5, 'delta\_f not exported', 'Units','normalized','HorizontalAlignment','center');
end
xlabel('Time (s)');
grid on;

saveas(f2, fullfile(figdir, 'heading_steering.png'));
fprintf('Saved heading_steering.png\n');

%% ---- Figure 3: Track overview (static, for reference) -----------------
f3 = figure('Visible','off','Position',[100 100 900 400]);
hold on; axis equal; grid on;
plot(path.xoutpath, path.youtpath, 'b-',  'LineWidth', 1.5, 'DisplayName','Outer boundary');
plot(path.xinpath,  path.yinpath,  'b--', 'LineWidth', 1.5, 'DisplayName','Inner boundary');
plot(path.xpath,    path.ypath,    'r-',  'LineWidth', 0.8, 'DisplayName','Center line');
xlabel('X (m)'); ylabel('Y (m)');
title('Oval Track: 900 m straights, R = 200 m, Width = 15 m');
legend('Location','southoutside','Orientation','horizontal');
saveas(f3, fullfile(figdir, 'track_overview.png'));
fprintf('Saved track_overview.png\n');

fprintf('\nAll figures saved to %s\n', figdir);
