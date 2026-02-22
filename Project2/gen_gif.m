% gen_gif.m — render animation to GIF for README
% Run from Project2/ directory by typing: gen_gif

clc; addpath(pwd);
init;
gentrack;

simout = sim('Project_2_Kinematic_Model.slx', 'StopTime', '260');

car_X=[]; car_Y=[]; car_psi=[];
if isprop(simout,'X');   d=simout.X;   if isprop(d,'Data'); car_X=d.Data;   else; car_X=d;   end; end
if isprop(simout,'Y');   d=simout.Y;   if isprop(d,'Data'); car_Y=d.Data;   else; car_Y=d;   end; end
if isprop(simout,'psi'); d=simout.psi; if isprop(d,'Data'); car_psi=d.Data; else; car_psi=d; end; end

figdir = fullfile(pwd, 'figures');
if ~exist(figdir,'dir'); mkdir(figdir); end
gifpath = fullfile(figdir, 'animation.gif');

% Downsample to ~60 frames for a reasonable GIF size
n      = length(car_X);
nframe = 80;
idx    = unique(round(linspace(1, n, nframe)));

L = 15; W = 5;

f = figure('Visible','off','Position',[100 100 900 480], 'Color','w');

for k = 1:length(idx)
    i = idx(k);
    clf(f);
    hold on; axis equal; box on; grid on;

    % Track
    plot(path.xoutpath, path.youtpath, 'b-',  'LineWidth',1.5);
    plot(path.xinpath,  path.yinpath,  'b-',  'LineWidth',1.5);
    plot(path.xpath,    path.ypath,    '--r', 'LineWidth',0.8);

    % Trail
    plot(car_X(1:i), car_Y(1:i), 'g-', 'LineWidth',1.5);

    % Vehicle patch
    psi = car_psi(i);
    corners = [-L/2 -W/2; -L/2 W/2; L/2 W/2; L/2 -W/2];
    Rm = [cos(psi) -sin(psi); sin(psi) cos(psi)];
    rc = (Rm*corners')' + [car_X(i) car_Y(i)];
    fill(rc(:,1), rc(:,2), [0.1 0.3 0.8], 'FaceAlpha',0.9, 'EdgeColor','k', 'LineWidth',1.2);

    axis([-250 1150 -150 550]);
    xlabel('X (m)'); ylabel('Y (m)');
    title(sprintf('v_x = 15 m/s  |  t = %.1f s', simout.tout(i)));
    drawnow;

    % Capture frame
    frame = getframe(f);
    im    = frame2im(frame);
    [imind, cm] = rgb2ind(im, 128);

    if k == 1
        imwrite(imind, cm, gifpath, 'gif', 'Loopcount',inf, 'DelayTime',0.06);
    else
        imwrite(imind, cm, gifpath, 'gif', 'WriteMode','append', 'DelayTime',0.06);
    end
end

close(f);
fprintf('Saved animation GIF (%d frames): %s\n', length(idx), gifpath);
