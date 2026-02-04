function [t, omega, tau, f_damping, cpu_time] = simulate_system_variable(J, b, omega0, tau_const, t_final, method, omega_freq)
%SIMULATE_SYSTEM_VARIABLE Simulates the rotational system with variable time step
%
% Inputs:
%   J - Moment of inertia (kg·m²)
%   b - Damping coefficient (N·m·s/rad)
%   omega0 - Initial angular velocity (rad/s)
%   tau_const - Constant applied torque (N·m) or [] for sinusoidal
%   t_final - Final simulation time (s)
%   method - Integration method ('ode45' or 'ode23tb')
%   omega_freq - Frequency for sinusoidal torque (rad/s), optional
%
% Outputs:
%   t - Time vector
%   omega - Angular velocity vector
%   tau - Applied torque vector
%   f_damping - Damping force vector
%   cpu_time - Computation time

tic;

% Set up ODE options
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);

% Determine if torque is sinusoidal
if isempty(tau_const)
    is_sinusoidal = true;
    if nargin < 7
        omega_freq = 1;
    end
else
    is_sinusoidal = false;
end

% Define the system dynamics
if is_sinusoidal
    % Sinusoidal torque: tau(t) = 100 * sin(omega_freq * t)
    dydt = @(t, omega) (100 * sin(omega_freq * t) - b * omega) / J;
else
    % Constant torque
    dydt = @(t, omega) (tau_const - b * omega) / J;
end

% Solve ODE
switch lower(method)
    case 'ode45'
        [t, omega] = ode45(dydt, [0 t_final], omega0, options);
        
    case 'ode23tb'
        [t, omega] = ode23tb(dydt, [0 t_final], omega0, options);
        
    otherwise
        error('Unknown ODE solver: %s', method);
end

% Transpose to match expected format
t = t';
omega = omega';

% Generate torque vector at the time points
if is_sinusoidal
    tau = 100 * sin(omega_freq * t);
else
    tau = tau_const * ones(size(t));
end

% Calculate damping force
f_damping = b * omega;

cpu_time = toc;
end
