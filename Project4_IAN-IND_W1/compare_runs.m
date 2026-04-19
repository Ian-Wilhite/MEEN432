%% compare_runs.m
% Generates comparison figures between the flat-track baseline and the
% elevated-track individual run.
%
% Expected files:
%   week2_results.mat    – flat track (team Week 2 baseline)
%   elevated_results.mat – elevated track (P4 Individual)
%
% Both files must contain: t, X, Y, vx, soc, results
% (same format as saved by p4_runsim.m)

clc;

% ── Load data ────────────────────────────────────────────────────────────
flat = load('week2_results.mat');

if ~exist('elevated_results.mat', 'file')
    error('elevated_results.mat not found. Run p4_runsim with the elevated model first.');
end
elev = load('elevated_results.mat');

L_flat = flat.results.track_length;
L_elev = elev.results.track_length;

% ── Figure 1: SOC vs Time ────────────────────────────────────────────────
figure(10); clf;

plot(flat.t, flat.soc * 100, 'b-',  'LineWidth', 1.6, 'DisplayName', ...
    sprintf('Flat track  (%.2f laps)', flat.results.laps));
hold on;
plot(elev.t, elev.soc * 100, 'r-',  'LineWidth', 1.6, 'DisplayName', ...
    sprintf('Elevated track  (%.2f laps)', elev.results.laps));

yline(10, 'k--', 'Min 10%', 'LineWidth', 1.0, 'HandleVisibility', 'off');
yline(95, 'k--', 'Max 95%', 'LineWidth', 1.0, 'HandleVisibility', 'off');

soc_drop_flat = (flat.soc(1)  - flat.soc(end))  * 100;
soc_drop_elev = (elev.soc(1)  - elev.soc(end))  * 100;
text(200, flat.soc(1)*100 - 0.3, sprintf('Drop: %.2f%%', soc_drop_flat), ...
    'Color', 'b', 'FontSize', 9);
text(200, elev.soc(1)*100 - 0.6, sprintf('Drop: %.2f%%', soc_drop_elev), ...
    'Color', 'r', 'FontSize', 9);

xlabel('Time (s)');
ylabel('State of Charge (%)');
title('Battery SOC vs Time — Flat vs Elevated Track');
legend('Location', 'southwest');
grid on;
ylim([min([flat.soc; elev.soc])*100 - 2, max([flat.soc(1); elev.soc(1)])*100 + 2]);

% ── Figure 2: Speed vs Arc-Length (first lap) ────────────────────────────
figure(11); clf;

% Extract first-lap indices (arc resets after each full loop)
arc_flat = flat.results.arc(:);
arc_elev = elev.results.arc(:);

idx_flat = first_lap_indices(arc_flat, L_flat);
idx_elev = first_lap_indices(arc_elev, L_elev);

vx_flat_mph = flat.vx(idx_flat) * 2.23694;
vx_elev_mph = elev.vx(idx_elev) * 2.23694;
s_flat      = arc_flat(idx_flat);
s_elev      = arc_elev(idx_elev);

plot(s_elev, vx_elev_mph, 'r-',  'LineWidth', 2.0, 'DisplayName', 'Elevated track');
hold on;
plot(s_flat, vx_flat_mph, 'b--', 'LineWidth', 1.4, 'DisplayName', 'Flat track');

% Mark straight / curve boundaries on the shared x-axis
% (uses path struct if available in workspace, otherwise estimates)
if exist('path', 'var')
    l_st   = path.l_st;
    l_curv = pi * path.radius;
else
    l_st   = 900;
    l_curv = pi * 200;
end
seg_ends = [l_st, l_st+l_curv, 2*l_st+l_curv, 2*l_st+2*l_curv];
seg_labels = {'Curve 1','Straight 2','Curve 2','(lap end)'};
for k = 1:numel(seg_ends)
    xline(seg_ends(k), 'k:', seg_labels{k}, 'LabelVerticalAlignment', 'bottom', ...
          'LineWidth', 0.8, 'HandleVisibility', 'off');
end

xlabel('Arc-length position (m)');
ylabel('Speed (mph)');
title('Vehicle Speed vs Arc-Length — First Lap Comparison');
legend('Location', 'best');
grid on;

% ── Save figures ─────────────────────────────────────────────────────────
out_dir = 'P4IND_Graphs';
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

exportgraphics(figure(10), fullfile(out_dir, 'compare_soc.png'),   'Resolution', 150);
exportgraphics(figure(11), fullfile(out_dir, 'compare_speed.png'), 'Resolution', 150);
fprintf('Figures saved to %s/\n', out_dir);

% ── Helper ───────────────────────────────────────────────────────────────
function idx = first_lap_indices(arc, L)
    % Return indices covering the first complete traversal (0 → L).
    % arc wraps back toward 0 when the vehicle crosses the start/finish.
    delta = diff(arc);
    wrap  = find(delta < -L/2, 1, 'first');
    if isempty(wrap)
        idx = 1:length(arc);   % fewer than one full lap — use all
    else
        idx = 1:wrap;
    end
end
