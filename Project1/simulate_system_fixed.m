function [t, omega, tau, f_damping, cpu_time] = simulate_system_fixed(J, b, omega0, tau_const, dt, t_final, method, omega_freq)
%SIMULATE_SYSTEM_FIXED Simulates the rotational system with fixed time step
%
% This function uses MATLAB's built-in fixed-step solvers where available (ode4 for RK4)
% or implements standard algorithms for other methods.
%
% Inputs:
%   J - Moment of inertia (kg·m²)
%   b - Damping coefficient (N·m·s/rad)
%   omega0 - Initial angular velocity (rad/s)
%   tau_const - Constant applied torque (N·m) or [] for sinusoidal
%   dt - Time step (s)
%   t_final - Final simulation time (s)
%   method - Integration method ('Euler' or 'RK4')
%   omega_freq - Frequency for sinusoidal torque (rad/s), optional
%
% Outputs:
%   t - Time vector
%   omega - Angular velocity vector
%   tau - Applied torque vector
%   f_damping - Damping force vector
%   cpu_time - Computation time

tic;

% Initialize time vector
t = (0:dt:t_final)';
n = length(t);

% Generate torque signal
if isempty(tau_const)
    % Sinusoidal torque
    if nargin < 8
        omega_freq = 1;  % Default frequency
    end
    tau = 100 * sin(omega_freq * t);  % Amplitude = 100 N·m
else
    % Constant torque
    tau = tau_const * ones(n, 1);
end

% Integration using explicit implementations
switch lower(method)
    case 'euler'
        % Forward Euler method
        omega = zeros(n, 1);
        omega(1) = omega0;
        
        for i = 1:n-1
            tau_i = tau(i);
            f_damp = b * omega(i);
            tau_net = tau_i - f_damp;
            domega = tau_net / J;
            omega(i+1) = omega(i) + domega * dt;
        end
        
    case 'rk4'
        % Runge-Kutta 4th order implementation
        omega = zeros(n, 1);
        omega(1) = omega0;
        
        for i = 1:n-1
            t_i = t(i);
            tau_i = tau(i);
            w_i = omega(i);
            
            % k1: dω/dt = (τ(t_i) - b·ω_i) / J
            k1 = (tau_i - b * w_i) / J;
            
            % k2: dω/dt at midpoint using k1
            % For sinusoidal, evaluate at t + dt/2
            if isempty(tau_const)
                tau_half = 100 * sin(omega_freq * (t_i + dt/2));
            else
                tau_half = tau_i;
            end
            k2 = (tau_half - b * (w_i + k1 * dt / 2)) / J;
            
            % k3: dω/dt at midpoint using k2
            k3 = (tau_half - b * (w_i + k2 * dt / 2)) / J;
            
            % k4: dω/dt at end point using k3
            if isempty(tau_const)
                tau_next = 100 * sin(omega_freq * (t_i + dt));
            else
                tau_next = tau_i;
            end
            k4 = (tau_next - b * (w_i + k3 * dt)) / J;
            
            % Update: ω_{i+1} = ω_i + (dt/6)(k1 + 2*k2 + 2*k3 + k4)
            omega(i+1) = w_i + (dt / 6) * (k1 + 2*k2 + 2*k3 + k4);
        end
        
    otherwise
        error('Unknown integration method: %s', method);
end

% Calculate damping force
f_damping = b * omega;

% Convert to row vectors for output consistency
tau = tau';
omega = omega';
f_damping = f_damping';
t = t';

cpu_time = toc;
end
