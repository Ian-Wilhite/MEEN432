%% gen_all_figs.m  — figure generator
% Loads figures/simdata.mat and regenerates all figures + GIF.
% No simulation required — run p4_runsim.m (or run_p4.sh) first.

clc;

p4_init;   % rebuild path, carData

figdir   = fullfile(pwd, 'figures');
datafile = fullfile(figdir, 'simdata.mat');
assert(exist(datafile,'file')==2, 'simdata.mat not found — run p4_runsim first.');

d = load(datafile);
t = d.t; X = d.X; Y = d.Y; vx = d.vx; soc = d.soc; psi = d.psi;
if isfield(d, 'lap_count')
    lap_count = d.lap_count;
else
    lap_count = compute_laps(X, Y, path);
end

first_lap_idx = find_first_lap_end(X, Y, path);

fprintf('Generating figures  (%.2f laps, %.0f s) ...\n', lap_count, t(end));

% ================================================================
%  FIG 1 — Lap 1 path
% ================================================================
fig1 = figure('Visible','off','Position',[100 100 900 500],'Color','w');
plot(path.xpath,    path.ypath,    'k-',  'LineWidth',1.5); hold on;
plot(path.xinpath,  path.yinpath,  'b--', 'LineWidth',1.0);
plot(path.xoutpath, path.youtpath, 'r--', 'LineWidth',1.0);
plot(X(1:first_lap_idx), Y(1:first_lap_idx), 'm-', 'LineWidth',2);
plot(X(1), Y(1), 'go', 'MarkerSize',10, 'MarkerFaceColor','g');
plot(X(first_lap_idx), Y(first_lap_idx), 'rs', 'MarkerSize',10, 'MarkerFaceColor','r');
legend('Track Centre','Inner Edge','Outer Edge','Vehicle Path','Start','Lap 1 End','Location','best');
xlabel('X (m)'); ylabel('Y (m)');
title(sprintf('P4 – Lap 1 Trajectory  (t = 0 – %.1f s)', t(first_lap_idx)));
axis equal; grid on;
save_fig(fig1, figdir, 'fig1_lap1_path');

% ================================================================
%  FIG 2 — All-laps path
% ================================================================
fig2 = figure('Visible','off','Position',[100 100 900 500],'Color','w');
plot(path.xpath,    path.ypath,    'k-',  'LineWidth',1.5); hold on;
plot(path.xinpath,  path.yinpath,  'b--', 'LineWidth',1.0);
plot(path.xoutpath, path.youtpath, 'r--', 'LineWidth',1.0);
scatter(X, Y, 2, t, 'filled');
cb = colorbar; cb.Label.String = 'Time (s)';
plot(X(1), Y(1), 'go', 'MarkerSize',10, 'MarkerFaceColor','g');
legend('Track Centre','Inner Edge','Outer Edge','Vehicle Path','Start','Location','best');
xlabel('X (m)'); ylabel('Y (m)');
title(sprintf('P4 – Full Run Trajectory  (%.2f laps in %.0f s)', lap_count, t(end)));
axis equal; grid on;
save_fig(fig2, figdir, 'fig2_all_laps_path');

% ================================================================
%  FIG 3 — Speed vs time (full)
% ================================================================
fig3 = figure('Visible','off','Position',[100 100 800 400],'Color','w');
plot(t, vx*2.23694, 'b-', 'LineWidth',1.5);
yline(carData.vxd*2.23694, 'r--', 'Desired', 'LineWidth',1.2);
xlabel('Time (s)'); ylabel('Speed (mph)');
title('Vehicle Speed vs Time');
grid on; legend('Actual Speed','Desired Speed');
save_fig(fig3, figdir, 'fig3_speed');

% ================================================================
%  FIG 3b — Speed vs time (first 250 s zoom)
% ================================================================
mask = t <= 250;
fig3b = figure('Visible','off','Position',[100 100 800 400],'Color','w');
plot(t(mask), vx(mask)*2.23694, 'b-', 'LineWidth',1.5);
yline(carData.vxd*2.23694, 'r--', 'Desired', 'LineWidth',1.2);
xlabel('Time (s)'); ylabel('Speed (mph)');
title('Vehicle Speed vs Time — First 250 s');
xlim([0 250]); grid on; legend('Actual Speed','Desired Speed');
save_fig(fig3b, figdir, 'fig3b_speed_zoom');

% ================================================================
%  FIG 4 — Battery SOC
% ================================================================
fig4 = figure('Visible','off','Position',[100 100 800 400],'Color','w');
plot(t, soc*100, 'g-', 'LineWidth',1.5);
yline(10, 'r--', 'Min 10%', 'LineWidth',1.2);
yline(95, 'b--', 'Max 95%', 'LineWidth',1.2);
xlabel('Time (s)'); ylabel('State of Charge (%)');
title('Battery SOC vs Time');
ylim([0 100]); grid on;
save_fig(fig4, figdir, 'fig4_soc');

% ================================================================
%  FIG 5 — Lateral offset vs lap position
% ================================================================
fig5 = plot_lateral_laps(X, Y, path, lap_count, figdir);
fprintf('  fig5_lateral_laps\n');

% ================================================================
%  GIF — animation
% ================================================================
fprintf('Generating animation GIF ...\n');
gifpath = fullfile(figdir, 'animation.gif');
n = length(X); nframe = 100;
idx = unique(round(linspace(1, n, nframe)));

L_car = 4.5; W_car = 2.0;
corners = [-L_car/2 -W_car/2; -L_car/2 W_car/2; L_car/2 W_car/2; L_car/2 -W_car/2];
ax_lim = [min(path.xoutpath)-30 max(path.xoutpath)+30 ...
          min(path.youtpath)-30 max(path.youtpath)+30];

fgif = figure('Visible','off','Position',[100 100 900 480],'Color','w');
for k = 1:length(idx)
    i = idx(k);
    clf(fgif); hold on; axis equal; box on; grid on;
    plot(path.xoutpath, path.youtpath, 'b-', 'LineWidth',1.5);
    plot(path.xinpath,  path.yinpath,  'b-', 'LineWidth',1.5);
    plot(path.xpath,    path.ypath,    'k--','LineWidth',0.8);
    plot(X(1:i), Y(1:i), 'm-', 'LineWidth',1.2);
    Rm = [cos(psi(i)) -sin(psi(i)); sin(psi(i)) cos(psi(i))];
    rc = (Rm*corners')' + [X(i) Y(i)];
    fill(rc(:,1), rc(:,2), [0.1 0.3 0.8], 'FaceAlpha',0.9, 'EdgeColor','k','LineWidth',1.2);
    axis(ax_lim);
    xlabel('X (m)'); ylabel('Y (m)');
    title(sprintf('Speed: %.1f mph  |  SOC: %.1f%%  |  t = %.1f s', vx(i)*2.23694, soc(i)*100, t(i)));
    drawnow;
    frame = getframe(fgif);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 128);
    if k == 1
        imwrite(imind, cm, gifpath, 'gif', 'Loopcount',inf, 'DelayTime',0.06);
    else
        imwrite(imind, cm, gifpath, 'gif', 'WriteMode','append', 'DelayTime',0.06);
    end
end
close(fgif);
fprintf('  animation.gif (%d frames)\n', length(idx));

fprintf('\nAll figures saved to %s\n', figdir);
ls_out = dir(fullfile(figdir,'*.png'));
for i = 1:length(ls_out)
    fprintf('  %s\n', ls_out(i).name);
end

% ================================================================
%  Helpers
% ================================================================
function save_fig(fig, figdir, name)
    exportgraphics(fig, fullfile(figdir, [name '.pdf']));
    exportgraphics(fig, fullfile(figdir, [name '.png']), 'Resolution',150);
    fprintf('  %s\n', name);
end

function fig = plot_lateral_laps(X, Y, path, lap_count, figdir)
    cx = path.xpath; cy = path.ypath;
    nv = length(X);
    ds_path = sqrt(diff(cx).^2 + diff(cy).^2);
    cumS = [0; cumsum(ds_path)]; L = cumS(end);
    half_w = path.width / 2;
    oval_cx = mean(cx); oval_cy = mean(cy);

    idx_ds = 1:5:nv;
    Xd = X(idx_ds); Yd = Y(idx_ds); nd = length(Xd);
    arc_pos = zeros(nd,1); lat_off = zeros(nd,1);

    for i = 1:nd
        d2 = (cx-Xd(i)).^2+(cy-Yd(i)).^2; [~,k] = min(d2);
        arc_pos(i) = cumS(k);

        dvx = Xd(i)-cx(k); dvy = Yd(i)-cy(k);

        % Path tangent at nearest centreline point
        if k < length(cx); tx=cx(k+1)-cx(k); ty=cy(k+1)-cy(k);
        else;               tx=cx(k)-cx(k-1); ty=cy(k)-cy(k-1); end
        tn=sqrt(tx^2+ty^2)+eps; tx=tx/tn; ty=ty/tn;

        % Left-facing normal to tangent
        nx=-ty; ny=tx;

        % True perpendicular CTE: project offset onto normal (removes along-track noise)
        cte_perp = dvx*nx + dvy*ny;

        % Sign convention: positive = toward outer wall (away from oval centre)
        owx=cx(k)-oval_cx; owy=cy(k)-oval_cy;
        own=sqrt(owx^2+owy^2)+eps;
        lat_off(i) = cte_perp * sign(nx*(owx/own) + ny*(owy/own));
    end

    % Detect lap boundaries: large backward jumps in arc_pos = crossed start/finish
    da = diff(arc_pos);
    lap_breaks = find(da < -L/2);              % indices just before each lap wrap
    seg_starts = [1;         lap_breaks+1];
    seg_ends   = [lap_breaks; nd          ];
    n_laps = length(seg_starts);

    s_c1s = path.l_st; s_c1e = path.l_st+pi*path.radius;
    s_c2s = s_c1e+path.l_st; s_c2e = min(s_c2s+pi*path.radius,L);

    cmap = lines(max(n_laps,1));
    fig = figure('Visible','off','Position',[100 100 1100 450],'Color','w');
    hold on;
    patch([0 L L 0],[ half_w  half_w  20  20],[0.8 0.2 0.2],'FaceAlpha',0.25,'EdgeColor','none');
    patch([0 L L 0],[-half_w -half_w -20 -20],[0.2 0.2 0.8],'FaceAlpha',0.25,'EdgeColor','none');
    for reg = [s_c1s s_c1e; s_c2s s_c2e]'
        patch([reg(1) reg(2) reg(2) reg(1)],[-20 -20 20 20], ...
              [0.95 0.95 0.7],'FaceAlpha',0.3,'EdgeColor','none');
    end
    yline( half_w,'r-','LineWidth',2.0);
    yline(-half_w,'b-','LineWidth',2.0);
    yline(0,'k--','LineWidth',1.0);
    for lap = 1:n_laps
        si = seg_starts(lap); ei = seg_ends(lap);
        if (ei-si) < 3; continue; end
        plot(arc_pos(si:ei), lat_off(si:ei), '-', ...
             'Color',[cmap(lap,:) 0.6], 'LineWidth',0.9);
    end
    text(path.l_st/2,         8.5,'Straight 1','HorizontalAlignment','center','FontSize',8,'Color',[0.35 0.35 0.35]);
    text((s_c1s+s_c1e)/2,     8.5,'Corner 1',  'HorizontalAlignment','center','FontSize',8,'Color',[0.5 0.3 0]);
    text((s_c1e+s_c2s)/2,     8.5,'Straight 2','HorizontalAlignment','center','FontSize',8,'Color',[0.35 0.35 0.35]);
    text((s_c2s+s_c2e)/2,     8.5,'Corner 2',  'HorizontalAlignment','center','FontSize',8,'Color',[0.5 0.3 0]);
    text(L*0.01, half_w+0.4, 'Outer wall (+7.5m)','FontSize',8,'Color',[0.75 0 0]);
    text(L*0.01,-half_w-0.9, 'Inner wall (-7.5m)','FontSize',8,'Color',[0 0 0.75]);
    text(L*0.01, 0.4,        'Centreline',         'FontSize',8,'Color',[0.2 0.2 0.2]);
    xlabel('Position along lap (m)');
    ylabel('Lateral offset from centreline (m)  [+ve = outer wall]');
    title(sprintf('Lateral Position vs Lap Position  (%.1f laps  |  track width = %.0f m)',lap_count,path.width));
    xlim([0 L]); ylim([-12.5 12.5]); grid on; box on;
    save_fig(fig, figdir, 'fig5_lateral_laps');
end

function laps = compute_laps(X, Y, path)
    cx = path.xpath; cy = path.ypath;
    nv = length(X);
    ds = sqrt(diff(cx).^2+diff(cy).^2);
    cumS = [0; cumsum(ds)]; L = cumS(end);
    arc = zeros(nv,1);
    for i = 1:nv
        d2 = (cx-X(i)).^2+(cy-Y(i)).^2; [~,k] = min(d2);
        arc(i) = cumS(k);
    end
    da = diff(arc); da(da<-L/2) = da(da<-L/2)+L;
    laps = sum(da(da>0)) / L;
end

function idx = find_first_lap_end(X, Y, path)
    cx = path.xpath; cy = path.ypath;
    nv = length(X);
    ds = sqrt(diff(cx).^2+diff(cy).^2);
    cumS = [0; cumsum(ds)]; L = cumS(end);
    arc = zeros(nv,1);
    for i = 1:nv
        d2 = (cx-X(i)).^2+(cy-Y(i)).^2; [~,k] = min(d2);
        arc(i) = cumS(k);
    end
    da = diff(arc); da(da<-L/2) = da(da<-L/2)+L;
    cum_arc = [0; cumsum(max(da,0))];
    idx = find(cum_arc>=L,1,'first');
    if isempty(idx); idx = nv; end
end
