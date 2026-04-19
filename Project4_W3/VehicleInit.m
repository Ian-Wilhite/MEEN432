%%% PROJECT 4 / WEEK 3 PARAMETERS %%%
% Vehicle parameters
datCar.I = 2000; % kg m^2
datCar.m = 1600; % kg

% Vehicle Drag Coefficients
datCar.C0 = 0.0041;
datCar.C1 = 0.000066;

Rho = 1.225;
A   = 2.6;
Cd  = 0.36;
datCar.C2 = 0.5 * Rho * A * Cd;

% Vehicle Initial Conditions
datCar.init.X0 = 0;
datCar.init.Y0 = 0;
datCar.init.vx0 = 0.001;
datCar.init.vy0 = 0;
datCar.init.omega0 = 0;
datCar.init.psi0 = 0;

% Tire Data
datCar.Calpha_f = 40000;
datCar.Calpha_r = 40000;
datCar.Fyfmax   = 40000*2/180*pi;
datCar.Fyrmax   = 40000*2/180*pi;
datCar.lr       = 1.5;
datCar.lf       = 1.0;
datCar.radius   = 0.3;
datCar.maxAlpha = 2 / 180 * pi;
datCar.Iw       = 0.5*7*(datCar.radius^2);

datCar.C_lambda   = 50;
datCar.lambda_max = 0.1;

% WEEK 3 requirement
datCar.tire_mu = 0.5;

datCar.init.omega0 = datCar.init.vx0 / datCar.radius;

% Transmission Data
datCar.gearRatio1 = 10.0;
datCar.gearRatio2 = 3.0;
datCar.gearRatio3 = 1.0;

datCar.FDRatio = 7.5;
datCar.FDeta   = 0.95;

% Brake
datCar.maxBrakeTorque = 50000;