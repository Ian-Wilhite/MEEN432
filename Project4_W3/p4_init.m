%% p4_init.m
% Project 4 Week 3 initialization
% Combined lateral + longitudinal EV model
% Goal: maximize laps in 60 minutes while satisfying constraints

clc;

% ======== Longitudinal Parameters ========== %
VehicleInit;
BatteryInit;
ElectricMotorInit;

% Project requirement: initial battery SOC = 80%
datBat.SOC_init = 0.80;

% ======== Lateral Parameters ========== %
carData.Inertia = datCar.I;
carData.Mass    = datCar.m;

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

% ====== Week 3 target / tuning candidates ====== %
carData.vxd_mph = 20;
carData.vxd     = carData.vxd_mph * 0.44704;

% Wider sweep for Week 3 optimization
carData.week3_speed_candidates_mph = 18:1:28;

% Constant block bus
carDataBus = Simulink.Bus.createObject(carData); %#ok<NASGU>

% ====== Generate Track ========= %
gentrack;

% Keep a tighter lookahead for better curve tracking
path.pure_pursuit_lookaheaddist = 7;
assignin('base','path',path);

% ======= Simulation settings ========= %
sim_stop_time = 3600;

fprintf('Project 4 Week 3 Initialization Complete\n');
fprintf('Vehicle mass          : %.0f kg\n', datCar.m);
fprintf('Friction coefficient  : %.2f\n', datCar.tire_mu);
fprintf('Initial battery SOC   : %.0f%%\n', datBat.SOC_init*100);
fprintf('Target speed          : %.1f m/s (%.1f mph)\n', carData.vxd, carData.vxd_mph);
fprintf('Lookahead distance    : %.1f m\n', path.pure_pursuit_lookaheaddist);
fprintf('Track length          : %.0f m per lap\n', path.total_length);
fprintf('Required simulation   : %.0f s (%.1f min)\n', sim_stop_time, sim_stop_time/60);