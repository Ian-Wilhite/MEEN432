% Generate and save Week 1 figures directly from starter MATLAB code.
% Run inside MATLAB from the Project2 directory:
%   run('scripts/save_figs_matlab.m');
% Requires only gentrack.m (no Simulink).

close all; clc;
addpath(pwd); % ensure gentrack is on path

% Build track using starter code
clear path;
gentrack; % populates 'path' and assigns to base; uses 900 m / 200 m / 15 m

figdir = fullfile(pwd,'figures');
if ~exist(figdir,'dir'); mkdir(figdir); end

%% Figure 1: track overview
f1 = figure('Visible','off');
hold on;
plot(path.xpath, path.ypath, '--r', 'DisplayName','Center line');
plot(path.xinpath, path.yinpath, 'b', 'DisplayName','Inner border');
plot(path.xoutpath, path.youtpath, 'b', 'DisplayName','Outer border');
axis equal;
xlabel('X (m)'); ylabel('Y (m)');
title('Project 2 Track: 900 m straights, 200 m radius, 15 m width');
legend('Location','southoutside','Orientation','horizontal');
grid on;
tightfig(f1);
saveas(f1, fullfile(figdir,'track_overview.png'));

%% Figure 2: sample vehicle poses
f2 = figure('Visible','off');
plot(path.xpath, path.ypath, '--k', 'Color',[0 0 0 0.4]);
hold on; axis equal; grid on;
xlabel('X (m)'); ylabel('Y (m)');
title('Sample vehicle poses along center line');

L = 15; W = 5;
sample_idx = round(linspace(1, numel(path.xpath), 8));
for idx = sample_idx
    x = path.xpath(idx); y = path.ypath(idx);
    nxt = mod(idx, numel(path.xpath)) + 1;
    psi = atan2(path.ypath(nxt) - y, path.xpath(nxt) - x);
    car = [-L/2 -W/2; -L/2 W/2; L/2 W/2; L/2 -W/2; -L/2 -W/2];
    R = [cos(psi) -sin(psi); sin(psi) cos(psi)];
    rcar = (R * car')' + [x y];
    plot(rcar(:,1), rcar(:,2), 'r-');
end

tightfig(f2);
saveas(f2, fullfile(figdir,'vehicle_poses.png'));

fprintf('Saved figures to %s\n', figdir);

function tightfig(fig)
% Minimal padding helper
set(fig,'Units','inches');
pos = get(fig,'Position');
tight = pos;
tight(3) = pos(3); tight(4) = pos(4);
set(fig,'PaperPositionMode','auto','InvertHardcopy','off');
end
