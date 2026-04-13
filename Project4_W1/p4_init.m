%% p4_init.m
% Combined lateral dynamics (P2) and longitudinal Dynamics (P3)

clc;

% ======== Longitudinal Paramters (P3) ==========%
VehicleInit; % loads car data
BatteryInit; % loads battery data
ElectricMotorInit;  % loads motor data

% Initial Battery SOC 80%
datBat.SOC_init = 0.80;

% ======== Lateral Paramters (P2) ==========%
%Getting Car Data/Chatacteristics
carData.Inertia   = datCar.I; % kg·m^2
carData.Mass      = datCar.m; % kg

% Initial conditions
carData.init.X0     = datCar.init.X0; % m
carData.init.Y0     = datCar.init.Y0; % m
carData.init.vx0    = datCar.init.vx0; % m/s  (small positive so no /0)
carData.init.vy0    = datCar.init.vy0; % m/s
carData.init.omega0 = datCar.init.omega0; % rad/s
carData.init.psi0   = datCar.init.psi0; % rad

% Tire parameters
carData.Calpha_f = datCar.Calpha_f;
carData.Calpha_r = datCar.Calpha_r;
carData.Fyfmax   = datCar.Fyfmax;
carData.Fyrmax   = datCar.Fyrmax;
carData.lr       = datCar.lr;
carData.lf       = datCar.lf;
carData.radius   = datCar.radius;
carData.maxAlpha = datCar.maxAlpha;

% Rough Desired velocity
carData.vxd = 9.5;  % ~21 mph
carData.vx_threshold1 = 0.1;

% Constant block bus
carDataBus = Simulink.Bus.createObject(carData);

% ====== Generate Track =========%
gentrack;   % creates 'path' and 'track' structs in base workspace

% ======= SIMULATION SETTINGS =========%
% roughly 300s for single lap testing
sim_stop_time = 3600;

fprintf('Project 4 Initialization Complete\n');
fprintf('Vehicle mass  : %.0f kg\n', datCar.m);
fprintf('Battery SOC0  : %.0f%%\n',  datBat.SOC_init*100);
fprintf('Desired speed : %.1f m/s (%.1f mph)\n', carData.vxd, carData.vxd/0.44704);
fprintf('Track length  : %.0f m per lap\n', path.total_length);
fprintf('Stop time     : %.0f s\n', sim_stop_time);
