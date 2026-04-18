function [h, theta, F_grade] = compute_elevation(x_pos, direction, vehicle_mass)
%COMPUTE_ELEVATION  Height, grade angle, and gravitational force for the P4 track.
%
%  The elevation formula (from project spec):
%    Height(x) = exp((x-450)/50) / (1 + exp((x-450)/50)) * 10  [m]
%
%  Applied to straight segments (x in [0, 900] m).
%  Curved sections are assumed flat (h=const, theta=0).
%
%  Inputs
%    x_pos         – x-coordinate of vehicle (m)
%    direction     – +1 for start straight (moving in +x, uphill)
%                   -1 for return straight (moving in -x, downhill)
%    vehicle_mass  – vehicle mass (kg)  [optional, default 1600]
%
%  Outputs
%    h        – height above datum (m)
%    theta    – signed grade angle (rad); positive = uphill in direction of travel
%    F_grade  – longitudinal grade force (N); positive opposes motion (uphill)

    if nargin < 3
        vehicle_mass = 1600;
    end
    g = 9.81;

    sig  = exp((x_pos - 450)/50) ./ (1 + exp((x_pos - 450)/50));
    dhdx = (10/50) .* sig .* (1 - sig);

    h       = sig * 10;
    theta   = atan(direction .* dhdx);
    F_grade = vehicle_mass * g * sin(theta);
end
