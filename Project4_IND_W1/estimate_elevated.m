%% estimate_elevated.m
% Post-processes the flat-track run to estimate elevated-track behaviour.
%
% The grade is mild (<3 deg) and the speed controller maintains target speed,
% so the trajectory and speed are unchanged. The main effect is on SOC:
%   - uphill: motor draws extra battery power  = F_grade * v / eta_drive
%   - downhill: regen recovers partial energy  = F_grade * v * eta_regen  (negative = reduces drain)
%
% Motor/regen efficiency taken from the operating region of the eta_val map (~0.85 each way).
% Results saved to elevated_results.mat in the same format as week2_results.mat.

clc;

%% Load flat run
flat     = load('week2_results.mat');
t        = flat.t(:);
vx       = flat.vx(:);
soc_flat = flat.soc(:);
arc      = flat.results.arc(:);

%% Build path struct with elevation data
p4_init;   % populates path.cumS, path.theta_grade, datCar, datBat

%% Battery pack energy capacity (J)
% OCV at operating SOC (~0.8): interpolate from cell OCV table
% Pack voltage = numSeries * OCV_per_cell
V_pack = datBat.numSeries * interp1(datBat.SOC, datBat.OCV, mean(soc_flat), 'linear', 'extrap');
E_pack_J = datBat.C * V_pack * 3600;  % C [Ah] * V [V] * 3600 [s/h] = J

fprintf('Pack voltage (avg)     : %.1f V\n', V_pack);
fprintf('Pack energy capacity   : %.2f kWh\n', E_pack_J / 3.6e6);

%% Motor / regen efficiencies (typical operating point from eta_val map)
eta_drive = 0.85;   % motor + inverter, traction mode
eta_regen = 0.85;   % motor + inverter, regeneration mode

%% Grade angle at each timestep via arc-length interpolation
theta_at_arc = interp1(path.cumS, path.theta_grade, arc, 'linear', 0);

%% Extra battery power due to grade  [W]
%  F_grade > 0  → uphill  → motor works harder  → battery drains faster
%  F_grade < 0  → downhill → motor works less / regen → battery drains slower
m = datCar.m;
g = datCar.g;
F_grade = m * g * sin(theta_at_arc);

dP_bat = zeros(size(F_grade));
uphill           = F_grade > 0;
downhill         = F_grade < 0;
dP_bat(uphill)   =  F_grade(uphill)   .* vx(uphill)   / eta_drive;
dP_bat(downhill) =  F_grade(downhill) .* vx(downhill) * eta_regen;  % negative

%% Cumulative extra energy and SOC correction
dE_bat   = cumtrapz(t, dP_bat);           % J
soc_elev = soc_flat - dE_bat / E_pack_J;

%% Report
fprintf('SOC drop   flat     : %.4f%%\n', (soc_flat(1) - soc_flat(end)) * 100);
fprintf('SOC drop   elevated : %.4f%%\n', (soc_elev(1) - soc_elev(end)) * 100);
fprintf('Extra SOC drain     : %.4f%%\n', ((soc_elev(1)-soc_elev(end)) - (soc_flat(1)-soc_flat(end))) * 100);
fprintf('Net extra energy    : %.2f kJ\n', dE_bat(end) / 1000);

%% Save — same struct layout as week2_results.mat
X       = flat.X;
Y       = flat.Y;
soc     = soc_elev;
results = flat.results;
save('elevated_results.mat', 't', 'X', 'Y', 'vx', 'soc', 'results');
fprintf('Saved elevated_results.mat\n');
