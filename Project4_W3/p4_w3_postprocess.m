%% p4_w3_postprocess.m
% Load week3_results.mat and generate README figures as PNGs

clc; close all;

% Load saved simulation data
load('week3_results.mat', 't', 'X', 'Y', 'vx', 'soc', 'results');

% Rebuild path struct
p4_init;

out_dir = 'P4W3_Graphs';
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

% ---- helpers ----
cx  = path.xpath(:);
cy  = path.ypath(:);
ds  = sqrt(diff(cx).^2 + diff(cy).^2);
cumS_path = [0; cumsum(ds)];
L   = cumS_path(end);

% Arc-length and CTE along run
nv  = length(X);
arc = zeros(nv,1);
cte = zeros(nv,1);
for i = 1:nv
    d2 = (cx - X(i)).^2 + (cy - Y(i)).^2;
    [d2min, idx] = min(d2);
    arc(i) = cumS_path(idx);
    cte(i) = sqrt(d2min) .* sign(1); % magnitude; add sign below
end

% Signed CTE (left = positive using cross product heuristic)
for i = 1:nv
    d2 = (cx - X(i)).^2 + (cy - Y(i)).^2;
    [~, idx] = min(d2);
    if idx < length(cx)
        tang_x = cx(idx+1) - cx(idx);
        tang_y = cy(idx+1) - cy(idx);
    else
        tang_x = cx(idx) - cx(idx-1);
        tang_y = cy(idx) - cy(idx-1);
    end
    cross = tang_x*(Y(i)-cy(idx)) - tang_y*(X(i)-cx(idx));
    cte(i) = cte(i) * sign(cross);
end

% Lap boundaries: detect arc wrap
delta_arc = diff(arc);
wrap_idx  = find(delta_arc < -L/2);

% ----------------------------------------------------------------
% Fig 1 — Lap 1 trajectory
% ----------------------------------------------------------------
if ~isempty(wrap_idx)
    i_end_lap1 = wrap_idx(1);
else
    i_end_lap1 = nv;
end

fig1 = figure('Visible','off');
plot(path.xpath, path.ypath, 'k-', 'LineWidth', 1.5); hold on;
plot(path.xinpath, path.yinpath, 'b--', 'LineWidth', 1.0);
plot(path.xoutpath, path.youtpath, 'r--', 'LineWidth', 1.0);
plot(X(1:i_end_lap1), Y(1:i_end_lap1), 'm-', 'LineWidth', 2.0);
scatter(X(1), Y(1), 60, 'g', 'filled', 'DisplayName', 'Start');
legend('Track Centre','Inner Edge','Outer Edge','Vehicle Path','Start','Location','best');
xlabel('X (m)'); ylabel('Y (m)');
title('Project 4 Week 3 – Lap 1 Trajectory');
axis equal; grid on;
saveas(fig1, fullfile(out_dir, 'fig1_lap1_path.png'));
fprintf('Saved fig1_lap1_path.png\n');

% ----------------------------------------------------------------
% Fig 2 — All laps, coloured by time
% ----------------------------------------------------------------
fig2 = figure('Visible','off');
plot(path.xpath, path.ypath, 'k-', 'LineWidth', 1.5); hold on;
plot(path.xinpath, path.yinpath, 'b--', 'LineWidth', 0.8);
plot(path.xoutpath, path.youtpath, 'r--', 'LineWidth', 0.8);
scatter(X, Y, 2, t, 'filled');
c = colorbar; c.Label.String = 'Time (s)';
colormap(jet);
xlabel('X (m)'); ylabel('Y (m)');
title('Project 4 Week 3 – Full Run Trajectory (coloured by time)');
axis equal; grid on;
saveas(fig2, fullfile(out_dir, 'fig2_all_laps_path.png'));
fprintf('Saved fig2_all_laps_path.png\n');

% ----------------------------------------------------------------
% Fig 3 — Speed vs Time (full run)
% ----------------------------------------------------------------
fig3 = figure('Visible','off');
plot(t, vx*2.23694, 'b-', 'LineWidth', 1.2); hold on;
yline(carData.vxd_mph, 'r--', sprintf('Target: %.0f mph', carData.vxd_mph), 'LineWidth', 1.2);
xlabel('Time (s)'); ylabel('Speed (mph)');
title('Project 4 Week 3 – Vehicle Speed vs Time');
grid on; legend('Actual Speed','Target Speed','Location','best');
saveas(fig3, fullfile(out_dir, 'fig3_speed.png'));
fprintf('Saved fig3_speed.png\n');

% ----------------------------------------------------------------
% Fig 3b — Speed vs Time (zoom: first 250 s)
% ----------------------------------------------------------------
mask250 = t <= 250;
fig3b = figure('Visible','off');
plot(t(mask250), vx(mask250)*2.23694, 'b-', 'LineWidth', 1.5); hold on;
yline(carData.vxd_mph, 'r--', sprintf('Target: %.0f mph', carData.vxd_mph), 'LineWidth', 1.2);
xlabel('Time (s)'); ylabel('Speed (mph)');
title('Project 4 Week 3 – Vehicle Speed vs Time (first 250 s)');
xlim([0 250]); grid on;
legend('Actual Speed','Target Speed','Location','best');
saveas(fig3b, fullfile(out_dir, 'fig3b_speed_zoom.png'));
fprintf('Saved fig3b_speed_zoom.png\n');

% ----------------------------------------------------------------
% Fig 4 — Battery SOC vs Time
% ----------------------------------------------------------------
fig4 = figure('Visible','off');
plot(t, soc*100, 'g-', 'LineWidth', 1.5); hold on;
yline(10, 'r--', 'Min 10%', 'LineWidth', 1.2);
yline(95, 'b--', 'Max 95%', 'LineWidth', 1.2);
xlabel('Time (s)'); ylabel('State of Charge (%)');
title('Project 4 Week 3 – Battery SOC vs Time');
ylim([0 100]); grid on;
legend('SOC','Min limit','Max limit','Location','best');
saveas(fig4, fullfile(out_dir, 'fig4_soc.png'));
fprintf('Saved fig4_soc.png\n');

% ----------------------------------------------------------------
% Fig 5 — Lateral position (CTE) vs arc-length within each lap
%   Use 10 m arc-length bins (mean CTE per bin) to eliminate
%   nearest-point lookup noise
% ----------------------------------------------------------------
bin_width = 10;
bin_edges = 0 : bin_width : L + bin_width;
bin_ctrs  = bin_edges(1:end-1) + bin_width/2;

fig5 = figure('Visible','off'); hold on;
colors  = lines(ceil(results.laps)+1);
lap_num = 0;
i_start = 1;

wrap_pts = [wrap_idx; nv];
for w = 1:length(wrap_pts)
    i_end_seg = wrap_pts(w);
    lap_num   = lap_num + 1;
    seg_arc   = arc(i_start:i_end_seg);
    seg_cte   = cte(i_start:i_end_seg);

    % Bin by arc position and take mean CTE per bin
    [~, ~, bin_id] = histcounts(seg_arc, bin_edges);
    n_bins    = length(bin_ctrs);
    cte_binned = nan(n_bins, 1);
    for b = 1:n_bins
        pts = seg_cte(bin_id == b);
        if ~isempty(pts); cte_binned(b) = mean(pts); end
    end
    valid = ~isnan(cte_binned);

    ci = mod(lap_num-1, size(colors,1)) + 1;
    plot(bin_ctrs(valid), cte_binned(valid), '-', ...
         'Color', colors(ci,:), 'LineWidth', 1.5, ...
         'DisplayName', sprintf('Lap %d', lap_num));
    i_start = i_end_seg + 1;
    if i_start > nv; break; end
end
yline(path.width/2,  'r--', '+Track limit', 'LineWidth', 1.2, 'HandleVisibility','off');
yline(-path.width/2, 'r--', '−Track limit', 'LineWidth', 1.2, 'HandleVisibility','off');
plot(nan, nan, 'r--', 'DisplayName', 'Track limit (±7.5 m)');
xlabel('Arc-length position within lap (m)');
ylabel('Cross-track error (m)');
title('Project 4 Week 3 – Lateral Position vs Arc-Length (per lap)');
grid on; legend('Location','best');
saveas(fig5, fullfile(out_dir, 'fig5_lateral_laps.png'));
fprintf('Saved fig5_lateral_laps.png\n');

fprintf('\nAll figures saved to %s/\n', out_dir);
