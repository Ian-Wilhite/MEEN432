%% p4_runsim.m  — simulation only
% Runs p4_model.slx and saves outputs to figures/simdata.mat.
% Call gen_all_figs.m separately (or use run_p4.sh) to generate figures.

clc; clear;

p4_init;   % carData, datCar, datBat, datMotor, path, track, sim_stop_time

figdir = fullfile(pwd, 'figures');
if ~exist(figdir, 'dir'); mkdir(figdir); end

fprintf('\nRunning p4_model.slx  (stop time = %.0f s) ...\n', sim_stop_time);
simout = sim('p4_model.slx', 'StopTime', num2str(sim_stop_time));
fprintf('Simulation complete.\n\n');

% ---- Extract signals ----
t = simout.tout;

X_ts = simout.X;   Y_ts = simout.Y;
if isa(X_ts,'timeseries'); X = X_ts.Data; Y = Y_ts.Data;
else;                      X = X_ts;      Y = Y_ts; end

if isprop(simout,'psi') && ~isempty(simout.psi)
    psi_ts = simout.psi;
    if isa(psi_ts,'timeseries'); psi = psi_ts.Data;
    else;                        psi = psi_ts; end
else
    dx = gradient(X); dy = gradient(Y);
    psi = atan2(dy, dx);
end

vx_ts  = simout.vx_out;
soc_ts = simout.soc_out;
if isa(vx_ts,'timeseries'); vx = vx_ts.Data; soc = soc_ts.Data;
else;                        vx = vx_ts;      soc = soc_ts; end

% ---- Stats ----
lap_count = compute_laps(X, Y, path);
fprintf('Laps completed : %.2f\n', lap_count);
fprintf('Final SOC      : %.1f%%  (started %.1f%%)\n', soc(end)*100, soc(1)*100);
fprintf('Avg speed      : %.1f m/s  (%.1f mph)\n', mean(vx), mean(vx)*2.23694);

% ---- Save ----
datafile = fullfile(figdir, 'simdata.mat');
save(datafile, 't', 'X', 'Y', 'vx', 'soc', 'psi', 'lap_count');
fprintf('\nData saved to %s\n', datafile);

% ================================================================
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
