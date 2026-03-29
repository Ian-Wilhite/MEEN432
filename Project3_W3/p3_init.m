%%% PROJECT 3 PARAMETERS %%%

clearvars -except DriveData Time cycleName
clc

% ---------------- Vehicle Parameters ----------------
datCar.I = 2000;       % kg m^2
datCar.m = 1600;       % kg

% Drag coefficients
datCar.C0 = 0.0041;    % N
datCar.C1 = 0.000066;  % N/(m/s)
Rho = 1.225;           % kg/m^3
A = 2.6;               % m^2
Cd = 0.36;
datCar.C2 = 0.5 * Rho * A * Cd;

% Initial conditions
datCar.init.X0 = 0;
datCar.init.Y0 = 0;
datCar.init.vx0 = 0.001;
datCar.init.vy0 = 0;
datCar.init.omega0 = 0;
datCar.init.psi0 = 0;

% Tire / wheel data
datCar.radius = 0.3;   % m
datCar.Iw = 0.5 * 7 * (datCar.radius^2);

% Transmission / driveline
datCar.gearRatio = 10.0;   % use one single-speed ratio for EV model
datCar.FDRatio   = 7.5;
datCar.FDeta     = 0.95;

% Brake
datCar.maxBrakeTorque = 50000;   % Nm

% ---------------- Motor Data ----------------
scaleFactor = 0.75;

datMotor.rpm = [0,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,11000,12000];
datMotor.maxtorque = [ ...
    0,0,0,0,0,0,0,0,0,0,0,0,0;
    280,280,280,220,160,110,95,75,60,55,50,40,0;
    280,280,275,255,240,180,140,125,95,75,70,50,0;
    280,280,275,260,250,220,180,150,125,100,80,70,0;
    280,280,275,260,250,230,200,175,140,120,100,75,0] * scaleFactor;

datMotor.vbus = [250,350,500,600,700] * scaleFactor;

datMotor.inertia = 0.5;  % kg-m^2

% Simple assumed efficiencies for Week 3 if map is not fully used
datMotor.eta_drive = 0.92;   % motoring efficiency
datMotor.eta_regen = 0.70;   % regen efficiency

% Regen settings
datMotor.regenEnable = 1;
datMotor.maxRegenTorque = 120;   % Nm, choose a realistic capped value

% ---------------- Battery Data ----------------
datBat.SOC_bp = [0, 0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0];
datBat.OCV_bp = [0, 3.1, 3.55, 3.68, 3.74, 3.76, 3.78, 3.85, 3.9, 3.95, 4.08, 4.15]; % V/cell
datBat.Rint_cell = 0.1695;   % ohm/cell
datBat.C_Ah = 150;           % Ah total pack basis
datBat.numSeries = 96;
datBat.numParallel = 74;

% Pack-level nominal values
datBat.Q_C = datBat.C_Ah * 3600;                         % Coulombs
datBat.Rint_pack = datBat.Rint_cell * datBat.numSeries / datBat.numParallel;
datBat.SOC_init = 0.90;

% Optional safety limits
datBat.SOC_min = 0.20;
datBat.SOC_max = 1.00;

% Conversion
mph2mps = 1609.344 / 3600;