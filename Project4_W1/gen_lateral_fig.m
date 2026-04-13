%% gen_lateral_fig.m
% Regenerate fig5_lateral_laps from saved simdata.mat without re-running the sim.
% Run from Project4/New/ directory.

clc;
p4_init;   % rebuild path struct

figdir = fullfile(pwd, 'figures');
data   = load(fullfile(figdir, 'simdata.mat'));

lap_count = compute_laps(data.X, data.Y, path);
fprintf('Lap count from saved data: %.2f\n', lap_count);

fig = plot_lateral_laps(data.X, data.Y, path, lap_count, figdir);
fprintf('Saved fig5_lateral_laps\n');

% ----------------------------------------------------------------
function fig = plot_lateral_laps(X, Y, path, lap_count, figdir)
    cx = path.xpath; cy = path.ypath;
    nv = length(X);
    ds_path = sqrt(diff(cx).^2 + diff(cy).^2);
    cumS    = [0; cumsum(ds_path)];
    L       = cumS(end);
    half_w  = path.width / 2;   % 7.5 m

    % Oval centroid — used for consistent outward sign convention
    oval_cx = mean(cx);
    oval_cy = mean(cy);

    % Downsample for speed
    idx_ds = 1:5:nv;
    Xd = X(idx_ds); Yd = Y(idx_ds);
    nd = length(Xd);

    arc_pos = zeros(nd,1);
    lat_off = zeros(nd,1);

    for i = 1:nd
        d2 = (cx - Xd(i)).^2 + (cy - Yd(i)).^2;
        [~,k] = min(d2);
        arc_pos(i) = cumS(k);

        % Cross-track error: perpendicular distance to centreline,
        % signed positive toward outer wall.
        dvx = Xd(i)-cx(k); dvy = Yd(i)-cy(k);
        if k < length(cx)
            tx = cx(k+1)-cx(k); ty = cy(k+1)-cy(k);
        else
            tx = cx(k)-cx(k-1); ty = cy(k)-cy(k-1);
        end
        tn = sqrt(tx^2+ty^2)+eps; tx=tx/tn; ty=ty/tn;
        nx = -ty; ny = tx;
        cte = dvx*nx + dvy*ny;
        owx = cx(k)-oval_cx; owy = cy(k)-oval_cy;
        own = sqrt(owx^2+owy^2)+eps;
        lat_off(i) = cte * sign(nx*(owx/own) + ny*(owy/own));
    end

    % Cumulative arc to assign lap numbers
    da = diff(arc_pos);
    da(da < -L/2) = da(da < -L/2) + L;
    cum_arc   = [0; cumsum(max(da, 0))];
    lap_num   = floor(cum_arc / L);
    arc_inlap = mod(cum_arc, L);
    n_laps    = max(lap_num) + 1;

    % ---- Plot ----
    cmap = lines(n_laps);
    fig  = figure('Visible','off','Position',[100 100 1100 450],'Color','w');
    hold on;

    % Wall shading
    patch([0 L L 0], [ half_w  half_w  half_w+2  half_w+2], [0.8 0.2 0.2],'FaceAlpha',0.25,'EdgeColor','none');
    patch([0 L L 0], [-half_w -half_w -half_w-2 -half_w-2], [0.2 0.2 0.8],'FaceAlpha',0.25,'EdgeColor','none');

    % Corner region shading (yellow tint)
    s_c1s = path.l_st;
    s_c1e = path.l_st + pi*path.radius;
    s_c2s = s_c1e + path.l_st;
    s_c2e = min(s_c2s + pi*path.radius, L);
    for reg = [s_c1s s_c1e; s_c2s s_c2e]'
        patch([reg(1) reg(2) reg(2) reg(1)], [-half_w-2 -half_w-2 half_w+2 half_w+2], ...
              [0.95 0.95 0.7],'FaceAlpha',0.3,'EdgeColor','none');
    end

    % Boundaries and centreline
    yline( half_w, 'r-', 'LineWidth', 2.0);
    yline(-half_w, 'b-', 'LineWidth', 2.0);
    yline(0,       'k--','LineWidth', 1.0);

    % One line per lap
    for lap = 0:n_laps-1
        mask = lap_num == lap;
        if sum(mask) < 3; continue; end
        plot(arc_inlap(mask), lat_off(mask), '-', ...
            'Color', [cmap(lap+1,:) 0.6], 'LineWidth', 0.9);
    end

    % Labels
    text(path.l_st/2,            half_w+1.4, 'Straight 1','HorizontalAlignment','center','FontSize',8,'Color',[0.35 0.35 0.35]);
    text((s_c1s+s_c1e)/2,        half_w+1.4, 'Corner 1',  'HorizontalAlignment','center','FontSize',8,'Color',[0.5 0.3 0]);
    text((s_c1e+s_c2s)/2,        half_w+1.4, 'Straight 2','HorizontalAlignment','center','FontSize',8,'Color',[0.35 0.35 0.35]);
    text((s_c2s+s_c2e)/2,        half_w+1.4, 'Corner 2',  'HorizontalAlignment','center','FontSize',8,'Color',[0.5 0.3 0]);
    text(L*0.01,  half_w+0.35, 'Outer wall','FontSize',8,'Color',[0.75 0 0]);
    text(L*0.01, -half_w-0.85, 'Inner wall','FontSize',8,'Color',[0 0 0.75]);
    text(L*0.01,  0.35,        'Centreline','FontSize',8,'Color',[0.2 0.2 0.2]);

    xlabel('Position along lap (m)');
    ylabel('Lateral offset from centreline (m)  [+ve = outer wall]');
    title(sprintf('Lateral Position vs Lap Position  (%.1f laps  |  track width = %.0f m)', ...
          lap_count, path.width));
    xlim([0 L]);  ylim([-12.5 12.5]);
    grid on; box on;

    exportgraphics(fig, fullfile(figdir, 'fig5_lateral_laps.pdf'));
    exportgraphics(fig, fullfile(figdir, 'fig5_lateral_laps.png'), 'Resolution', 150);
end

function laps = compute_laps(X, Y, path)
    cx = path.xpath; cy = path.ypath;
    nv = length(X);
    ds = sqrt(diff(cx).^2 + diff(cy).^2);
    cumS = [0; cumsum(ds)]; L = cumS(end);
    arc = zeros(nv,1);
    for i = 1:nv
        d2 = (cx-X(i)).^2+(cy-Y(i)).^2; [~,k] = min(d2);
        arc(i) = cumS(k);
    end
    da = diff(arc); da(da<-L/2) = da(da<-L/2)+L;
    laps = sum(da(da>0)) / L;
end
