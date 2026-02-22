% gen_figs_week2.m — run from Project2 directory
clc;
addpath(pwd);

init;
gentrack;   % creates 'path' variable

figdir = fullfile(pwd, 'figures');
if ~exist(figdir,'dir'); mkdir(figdir); end

fprintf('Running Simulink model...\n');
% Full lap ≈ 3057 m at 25 m/s ≈ 122 s; run 160 s to capture at least one complete lap
simout = sim('Project_2_Kinematic_Model.slx', 'StopTime', '160');
fprintf('Simulation complete: %d time steps\n', length(simout.tout));

% Extract signals from simout
signames = {'X','Y','psi'};
car_X = []; car_Y = []; car_psi = [];
if isprop(simout,'X')
    d = simout.X; if isprop(d,'Data'); car_X = d.Data; else; car_X = d; end
end
if isprop(simout,'Y')
    d = simout.Y; if isprop(d,'Data'); car_Y = d.Data; else; car_Y = d; end
end
if isprop(simout,'psi')
    d = simout.psi; if isprop(d,'Data'); car_psi = d.Data; else; car_psi = d; end
end
t = simout.tout;

fprintf('X: %d points, Y: %d points, psi: %d points\n', ...
    length(car_X), length(car_Y), length(car_psi));

%% Figure 1: Track + simulated vehicle path
f1 = figure('Visible','off','Position',[100 100 1000 500]);
hold on; axis equal; grid on;
plot(path.xoutpath, path.youtpath, 'b-',  'LineWidth',1.5, 'DisplayName','Track boundary');
plot(path.xinpath,  path.yinpath,  'b-',  'LineWidth',1.5, 'HandleVisibility','off');
plot(path.xpath,    path.ypath,    '--r', 'LineWidth',0.8, 'DisplayName','Center line');
if ~isempty(car_X) && ~isempty(car_Y)
    plot(car_X, car_Y, 'g-', 'LineWidth',2.0, 'DisplayName','Vehicle path');
    plot(car_X(1), car_Y(1), 'go', 'MarkerSize',10, 'MarkerFaceColor','g', 'DisplayName','Start');
    % Draw vehicle patches at several time steps
    if ~isempty(car_psi)
        L = 15; W = 5;
        npts = length(car_X);
        poses = unique(round(linspace(1, npts, 8)));
        for k = poses
            xk = car_X(k); yk = car_Y(k); pk = car_psi(k);
            corners = [-L/2 -W/2; -L/2 W/2; L/2 W/2; L/2 -W/2];
            Rk = [cos(pk) -sin(pk); sin(pk) cos(pk)];
            rc = (Rk * corners')' + [xk yk];
            fill(rc(:,1), rc(:,2), [0.1 0.3 0.8], 'FaceAlpha',0.5, ...
                 'EdgeColor','k', 'HandleVisibility','off');
        end
    end
end
xlabel('X (m)'); ylabel('Y (m)');
title(sprintf('Simulated Vehicle Path — v_x = %.0f m/s (%.0f km/h)', ...
    carData.vxd, carData.vxd * 3.6));
legend('Location','southoutside','Orientation','horizontal','NumColumns',3);
saveas(f1, fullfile(figdir,'vehicle_path.png'));
fprintf('Saved vehicle_path.png\n');

%% Figure 2: Heading over time
if ~isempty(car_psi) && ~isempty(t)
    f2 = figure('Visible','off','Position',[100 100 900 350]);
    plot(t, car_psi * 180/pi, 'b-', 'LineWidth',1.5);
    xlabel('Time (s)'); ylabel('\psi (deg)');
    title('Vehicle Heading \psi vs. Time');
    grid on;
    saveas(f2, fullfile(figdir,'heading.png'));
    fprintf('Saved heading.png\n');
end

%% Figure 3: Track overview (static)
f3 = figure('Visible','off','Position',[100 100 1000 400]);
hold on; axis equal; grid on;
plot(path.xoutpath, path.youtpath, 'b-',  'LineWidth',1.5, 'DisplayName','Outer boundary');
plot(path.xinpath,  path.yinpath,  'b--', 'LineWidth',1.5, 'DisplayName','Inner boundary');
plot(path.xpath,    path.ypath,    'r-',  'LineWidth',0.8, 'DisplayName','Center line');
xlabel('X (m)'); ylabel('Y (m)');
title('Oval Track: 900 m straights, R = 200 m, Width = 15 m');
legend('Location','southoutside','Orientation','horizontal');
saveas(f3, fullfile(figdir,'track_overview.png'));
fprintf('Saved track_overview.png\n');

fprintf('\nAll figures saved to: %s\n', figdir);
