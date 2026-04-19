%% p4_init.m
% Project 4 Individual – Week 1 initialization
% Combined lateral dynamics (P2) and longitudinal dynamics (P3) with elevation

clc;

% ======== Longitudinal Parameters ========== %
VehicleInit;        % loads datCar
BatteryInit;        % loads datBat
ElectricMotorInit;  % loads datMotor

% Project requirement: initial battery SOC = 80%
datBat.SOC_init = 0.80;

% ======== Lateral Parameters ========== %
carData.Inertia   = datCar.I;     % kg*m^2
carData.Mass      = datCar.m;     % kg

% Initial conditions
carData.init.X0     = datCar.init.X0;
carData.init.Y0     = datCar.init.Y0;
carData.init.vx0    = datCar.init.vx0;
carData.init.vy0    = datCar.init.vy0;
carData.init.omega0 = datCar.init.omega0;
carData.init.psi0   = datCar.init.psi0;

% Tire / geometry parameters
carData.Calpha_f = datCar.Calpha_f;
carData.Calpha_r = datCar.Calpha_r;
carData.Fyfmax   = datCar.Fyfmax;
carData.Fyrmax   = datCar.Fyrmax;
carData.lr       = datCar.lr;
carData.lf       = datCar.lf;
carData.radius   = datCar.radius;
carData.maxAlpha = datCar.maxAlpha;
carData.vx_threshold1 = 0.1;

% ====== Speed target ====== %
% Start with a safer speed that is more likely to stay on track.
carData.vxd_mph = 20;
carData.vxd     = carData.vxd_mph * 0.44704;

% Candidate speeds for optional tuning sweep
carData.week2_speed_candidates_mph = [20 21 22 23 24 25];

% Constant block bus
carDataBus = Simulink.Bus.createObject(carData);

% ====== Generate Track ========= %
gentrack;   % creates 'path' and 'track' structs in base workspace

% Tighter lookahead for better curve tracking
path.pure_pursuit_lookaheaddist = 7;

% Push updated path back to base workspace for Simulink
assignin('base','path',path);

% ======= Elevation lookup table for Simulink ========= %
% 1-D lookup: arc-length position (m) → grade angle (rad)
% Use with a Simulink 1-D Lookup Table block:
%   breakpoints = grade_lut(:,1),  table data = grade_lut(:,2)
grade_lut = [path.cumS, path.theta_grade];  % [N x 2]

% Maximum grade magnitude and corresponding force for reference
datCar.g = 9.81;  % m/s^2
[max_grade_rad, idx_mg] = max(abs(path.theta_grade));
max_grade_force_N = datCar.m * datCar.g * sin(max_grade_rad);
fprintf('Max grade angle       : %.4f rad (%.2f deg)\n', ...
        max_grade_rad, max_grade_rad*180/pi);
fprintf('Max grade force       : %.1f N (%.1f%% of weight)\n', ...
        max_grade_force_N, max_grade_force_N/(datCar.m*datCar.g)*100);

% ======= Simulation settings ========= %
% Project requirement: total time of 60 minutes
sim_stop_time = 3600;

fprintf('Project 4 Individual Week 1 Initialization Complete\n');
fprintf('Vehicle mass          : %.0f kg\n', datCar.m);
fprintf('Initial battery SOC   : %.0f%%\n', datBat.SOC_init*100);
fprintf('Target speed          : %.1f m/s (%.1f mph)\n', carData.vxd, carData.vxd_mph);
fprintf('Lookahead distance    : %.1f m\n', path.pure_pursuit_lookaheaddist);
fprintf('Track length          : %.0f m per lap\n', path.total_length);
fprintf('Required simulation   : %.0f s (%.1f min)\n', sim_stop_time, sim_stop_time/60);