%%% PROJECT 3 PARAMETERS %%% YOU SHOULD BE USING ONLY THESE PARAMETERS FOR
% Vehicle parameters
datCar.I = 2000; %kg m^2  -  Car Inertia
datCar.m = 1600; %kg      -  Car Mass 

%Vehicle Drag Coefficients

datCar.C0 = 0.0041;         %in N
datCar.C1 = 0.000066;       %in N/(m/s)
% parameters for calculation of C2
Rho =1.225;          %Kg/m^3
A  = 2.6;           % m^2 (projected area)
Cd = 0.36;           %Aerodynamic drag coefficient
datCar.C2 = 0.5*Rho*A*Cd;   % in N/(m/s)^2

% Vehicle Initial Conditions

datCar.init.X0 = 0;       % Initial X position of the vehicle (track start)
datCar.init.Y0 = 0;       % Initial Y position of the vehicle (track start)
datCar.init.vx0 = 0.001;  % Initial velocity of vehicle in body x axis (m/s)
datCar.init.vy0 = 0;      % Initial velocity of vehilce in body y axis (m/s)
datCar.init.omega0 = 0;   % Initial yaw rate of vehicle (rad/s)
datCar.init.psi0 = 0;     % Initial yaw of vehicle (rad)

%Tire Data
datCar.Calpha_f = 40000; % N/rad  initial tire stiffness - front
datCar.Calpha_r = 40000; % N/rad  initial tire stiffness - rear
datCar.Fyfmax = 40000*2/180*pi;  % maximum front tire force achievable (N)
datCar.Fyrmax = 40000*2/180*pi;  % maximum rear tire force achievable (N)
datCar.lr = 1.5; % m  distance of center of mass to front axle
datCar.lf = 1.0; % m  distance of center of mass to rear axle
datCar.radius = 0.3; %m  tire radius
datCar.maxAlpha = 2 / 180 * pi;  % maximum alpha value to be used in calc (rad)
datCar.Iw = 0.5*7*(datCar.radius^2);

datCar.C_lambda = 50; % longitudinal stiffness (N/Kg)
datCar.lambda_max = 0.1; %maximum tire slip ratio before saturation
datCar.tire_mu = 1.0;


datCar.init.omega0 = datCar.init.vx0 / datCar.radius;

%Transmission Data
datCar.gearRatio1 = 10.0;
datCar.gearRatio2 = 3.0;
datCar.gearRatio3 = 1.0;

datCar.FDRatio = 7.5;
datCar.FDeta = 0.95;

%Brake
datCar.maxBrakeTorque = 50000;

%{
track Data
track.radius = 200; %m radius of curved sections
track.width = 15;   %m width of track
track.l_st = 900;   %m length of straightaways
%}

