%% MEEN432 Project 1: Inertia-Rotational Damper System Simulation
% This script simulates a rotational system consisting of an inertia element
% and a viscous damper under various conditions and integration methods.

clear all; close all; clc;

%% System Parameters
% Initial conditions
omega0_vals = [10.0, 1, 0.1, 0.0];  % rad/s

% Rotational Inertia (kg·m²)
J_vals = [100, 1, 0.01];

% Damping coefficient (N·m·s/rad)
b_vals = [10, 1, 0.1];

% Applied torque - Constant (N·m)
tau_const_vals = [0.0, 100];

% Applied torque - Sinusoidal frequency (rad/s)
omega_freq_vals = [0.1, 1, 10, 100];

% Integration parameters
dt_fixed = [0.001, 0.01, 0.1, 1];  % Fixed time steps (s)
t_final = 25;  % Total simulation time (s)

% Integration methods
fixed_methods = {'Euler', 'RK4'};
variable_methods = {'ode45', 'ode23tb'};

%% Storage for results
results = struct();
results.fixed = struct();
results.variable = struct();

%% Fixed Time Step Simulations
fprintf('Running fixed time step simulations...\n');

sim_count = 0;
for j = 1:length(J_vals)
    J = J_vals(j);
    for b_idx = 1:length(b_vals)
        b = b_vals(b_idx);
        
        % Constant torque cases
        for tau_idx = 1:length(tau_const_vals)
            tau_const = tau_const_vals(tau_idx);
            
            for w0_idx = 1:length(omega0_vals)
                omega0 = omega0_vals(w0_idx);
                
                for dt_idx = 1:length(dt_fixed)
                    dt = dt_fixed(dt_idx);
                    
                    for method_idx = 1:length(fixed_methods)
                        method = fixed_methods{method_idx};
                        sim_count = sim_count + 1;
                        
                        % Run simulation
                        [t, omega, tau, f_damping, cpu_time] = simulate_system_fixed(...
                            J, b, omega0, tau_const, dt, t_final, method);
                        
                        % Calculate theoretical response (step input)
                        % Only valid for non-zero constant torque
                        if tau_const == 0
                            % For zero torque, theoretical solution is exponential decay
                            omega_theory = omega0 * exp(-b/J * t);
                        else
                            % For constant torque step input
                            omega_theory = omega0 * exp(-b/J * t) + (tau_const/b) * (1 - exp(-b/J * t));
                        end
                        
                        % Calculate maximum error
                        max_error = max(abs(omega - omega_theory));
                        
                        % Store results
                        key = sprintf('J%.2f_b%.1f_tau%.1f_w0%.1f_dt%.4f_%s', ...
                            J, b, tau_const, omega0, dt, method);
                        key = strrep(key, '.', '_');  % Replace dots with underscores for valid field name
                        
                        results.fixed.(key).t = t;
                        results.fixed.(key).omega = omega;
                        results.fixed.(key).tau = tau;
                        results.fixed.(key).f_damping = f_damping;
                        results.fixed.(key).cpu_time = cpu_time;
                        results.fixed.(key).max_error = max_error;
                        results.fixed.(key).dt = dt;
                        results.fixed.(key).method = method;
                        results.fixed.(key).J = J;
                        results.fixed.(key).b = b;
                        results.fixed.(key).tau_const = tau_const;
                        results.fixed.(key).omega0 = omega0;
                        
                        fprintf('  Completed: J=%.2f, b=%.1f, tau=%.1f, dt=%.4f, method=%s (CPU: %.4f s)\n', ...
                            J, b, tau_const, dt, method, cpu_time);
                    end
                end
            end
        end
        
        % Sinusoidal torque cases
        for freq_idx = 1:length(omega_freq_vals)
            omega_freq = omega_freq_vals(freq_idx);
            
            for w0_idx = 1:length(omega0_vals)
                omega0 = omega0_vals(w0_idx);
                
                for dt_idx = 1:length(dt_fixed)
                    dt = dt_fixed(dt_idx);
                    
                    for method_idx = 1:length(fixed_methods)
                        method = fixed_methods{method_idx};
                        sim_count = sim_count + 1;
                        
                        % Run simulation with sinusoidal torque
                        [t, omega, tau, f_damping, cpu_time] = simulate_system_fixed(...
                            J, b, omega0, [], dt, t_final, method, omega_freq);
                        
                        % Store results
                        key = sprintf('J%.2f_b%.1f_sine_w0%.1f_freq%.2f_dt%.4f_%s', ...
                            J, b, omega0, omega_freq, dt, method);
                        key = strrep(key, '.', '_');  % Replace dots with underscores for valid field name
                        
                        results.fixed.(key).t = t;
                        results.fixed.(key).omega = omega;
                        results.fixed.(key).tau = tau;
                        results.fixed.(key).f_damping = f_damping;
                        results.fixed.(key).cpu_time = cpu_time;
                        results.fixed.(key).dt = dt;
                        results.fixed.(key).method = method;
                        results.fixed.(key).J = J;
                        results.fixed.(key).b = b;
                        results.fixed.(key).omega_freq = omega_freq;
                        results.fixed.(key).omega0 = omega0;
                        results.fixed.(key).is_sinusoidal = true;
                        
                        fprintf('  Completed: J=%.2f, b=%.1f, sine(%.2f rad/s), dt=%.4f, method=%s (CPU: %.4f s)\n', ...
                            J, b, omega_freq, dt, method, cpu_time);
                    end
                end
            end
        end
    end
end

%% Variable Time Step Simulations
fprintf('\nRunning variable time step simulations...\n');

for j = 1:length(J_vals)
    J = J_vals(j);
    for b_idx = 1:length(b_vals)
        b = b_vals(b_idx);
        
        % Constant torque cases
        for tau_idx = 1:length(tau_const_vals)
            tau_const = tau_const_vals(tau_idx);
            
            for w0_idx = 1:length(omega0_vals)
                omega0 = omega0_vals(w0_idx);
                
                for method_idx = 1:length(variable_methods)
                    method = variable_methods{method_idx};
                    
                    % Run simulation
                    [t, omega, tau, f_damping, cpu_time] = simulate_system_variable(...
                        J, b, omega0, tau_const, t_final, method);
                    
                    % Calculate theoretical response (step input)
                    % Only valid for non-zero constant torque
                    if tau_const == 0
                        % For zero torque, theoretical solution is exponential decay
                        omega_theory = omega0 * exp(-b/J * t);
                    else
                        % For constant torque step input
                        omega_theory = omega0 * exp(-b/J * t) + (tau_const/b) * (1 - exp(-b/J * t));
                    end
                    
                    % Calculate maximum error
                    max_error = max(abs(omega - omega_theory));
                    
                    % Store results
                    key = sprintf('J%.2f_b%.1f_tau%.1f_w0%.1f_%s', ...
                        J, b, tau_const, omega0, method);
                    key = strrep(key, '.', '_');  % Replace dots with underscores for valid field name
                    
                    results.variable.(key).t = t;
                    results.variable.(key).omega = omega;
                    results.variable.(key).tau = tau;
                    results.variable.(key).f_damping = f_damping;
                    results.variable.(key).cpu_time = cpu_time;
                    results.variable.(key).max_error = max_error;
                    results.variable.(key).method = method;
                    results.variable.(key).J = J;
                    results.variable.(key).b = b;
                    results.variable.(key).tau_const = tau_const;
                    results.variable.(key).omega0 = omega0;
                    
                    fprintf('  Completed: J=%.2f, b=%.1f, tau=%.1f, method=%s (CPU: %.4f s)\n', ...
                        J, b, tau_const, method, cpu_time);
                end
            end
        end
        
        % Sinusoidal torque cases
        for freq_idx = 1:length(omega_freq_vals)
            omega_freq = omega_freq_vals(freq_idx);
            
            for w0_idx = 1:length(omega0_vals)
                omega0 = omega0_vals(w0_idx);
                
                for method_idx = 1:length(variable_methods)
                    method = variable_methods{method_idx};
                    
                    % Run simulation with sinusoidal torque
                    [t, omega, tau, f_damping, cpu_time] = simulate_system_variable(...
                        J, b, omega0, [], t_final, method, omega_freq);
                    
                    % Store results
                    key = sprintf('J%.2f_b%.1f_sine_w0%.1f_freq%.2f_%s', ...
                        J, b, omega0, omega_freq, method);
                    key = strrep(key, '.', '_');  % Replace dots with underscores for valid field name
                    
                    results.variable.(key).t = t;
                    results.variable.(key).omega = omega;
                    results.variable.(key).tau = tau;
                    results.variable.(key).f_damping = f_damping;
                    results.variable.(key).cpu_time = cpu_time;
                    results.variable.(key).method = method;
                    results.variable.(key).J = J;
                    results.variable.(key).b = b;
                    results.variable.(key).omega_freq = omega_freq;
                    results.variable.(key).omega0 = omega0;
                    results.variable.(key).is_sinusoidal = true;
                    
                    fprintf('  Completed: J=%.2f, b=%.1f, sine(%.2f rad/s), method=%s (CPU: %.4f s)\n', ...
                        J, b, omega_freq, method, cpu_time);
                end
            end
        end
    end
end

%% Save results to workspace
save('project1_results.mat', 'results');
fprintf('\nResults saved to project1_results.mat\n');

%% Generate comparison plots
fprintf('\nGenerating plots...\n');

% Collect data for constant torque cases (step input) from fixed methods
const_torque_fixed_data = struct();
keys = fieldnames(results.fixed);
for k = 1:length(keys)
    if ~isfield(results.fixed.(keys{k}), 'is_sinusoidal')
        method = results.fixed.(keys{k}).method;
        dt = results.fixed.(keys{k}).dt;
        if ~isfield(const_torque_fixed_data, method)
            const_torque_fixed_data.(method).dts = [];
            const_torque_fixed_data.(method).errors = [];
            const_torque_fixed_data.(method).cpu_times = [];
        end
        const_torque_fixed_data.(method).dts = [const_torque_fixed_data.(method).dts, dt];
        const_torque_fixed_data.(method).errors = [const_torque_fixed_data.(method).errors, results.fixed.(keys{k}).max_error];
        const_torque_fixed_data.(method).cpu_times = [const_torque_fixed_data.(method).cpu_times, results.fixed.(keys{k}).cpu_time];
    end
end

% Collect data for constant torque cases from variable methods
const_torque_var_data = struct();
keys = fieldnames(results.variable);
for k = 1:length(keys)
    if ~isfield(results.variable.(keys{k}), 'is_sinusoidal')
        method = results.variable.(keys{k}).method;
        if ~isfield(const_torque_var_data, method)
            const_torque_var_data.(method).cpu_times = [];
            const_torque_var_data.(method).errors = [];
        end
        const_torque_var_data.(method).cpu_times = [const_torque_var_data.(method).cpu_times, results.variable.(keys{k}).cpu_time];
        const_torque_var_data.(method).errors = [const_torque_var_data.(method).errors, results.variable.(keys{k}).max_error];
    end
end

% Plot 1: Max error vs time step for fixed methods (using medians)
figure('Name', 'Max Error vs Time Step');
hold on;
methods_data = fieldnames(const_torque_fixed_data);
for m = 1:length(methods_data)
    method = methods_data{m};
    
    % Group by timestep and calculate median error
    unique_dts = unique(const_torque_fixed_data.(method).dts);
    median_errors = [];
    
    for dt_idx = 1:length(unique_dts)
        dt_val = unique_dts(dt_idx);
        idx = const_torque_fixed_data.(method).dts == dt_val;
        median_errors = [median_errors, median(const_torque_fixed_data.(method).errors(idx))];
    end
    
    plot(unique_dts, median_errors, 'o-', 'DisplayName', method, 'LineWidth', 2, 'MarkerSize', 8);
end
xlabel('Time Step (s)');
ylabel('Max Error (rad/s)');
title('Maximum Simulation Error vs Time Step (Median Values)');
legend('Location', 'best');
grid on;
set(gca, 'XScale', 'log', 'YScale', 'log');
saveas(gcf, 'assets/plot1_error_vs_timestep.png');
saveas(gcf, 'assets/plot1_error_vs_timestep.fig');

% Plot 2: CPU time vs time step for fixed methods (using medians)
figure('Name', 'CPU Time vs Time Step');
hold on;
methods_data = fieldnames(const_torque_fixed_data);
for m = 1:length(methods_data)
    method = methods_data{m};
    
    % Group by timestep and calculate median CPU time
    unique_dts = unique(const_torque_fixed_data.(method).dts);
    median_times = [];
    
    for dt_idx = 1:length(unique_dts)
        dt_val = unique_dts(dt_idx);
        idx = const_torque_fixed_data.(method).dts == dt_val;
        median_times = [median_times, median(const_torque_fixed_data.(method).cpu_times(idx))];
    end
    
    plot(unique_dts, median_times, 'o-', 'DisplayName', method, 'LineWidth', 2, 'MarkerSize', 8);
end
xlabel('Time Step (s)');
ylabel('CPU Time (s)');
title('CPU Time vs Time Step (Median Values)');
legend('Location', 'best');
grid on;
set(gca, 'XScale', 'log', 'YScale', 'log');
saveas(gcf, 'assets/plot2_cputime_vs_timestep.png');
saveas(gcf, 'assets/plot2_cputime_vs_timestep.fig');

% Plot 3: Max error vs CPU time (all methods)
figure('Name', 'Max Error vs CPU Time');
hold on;
colors = [0 0 1; 1 0 0; 0 1 0; 1 0 1];

% Fixed time step methods
fixed_methods_list = fieldnames(const_torque_fixed_data);
for m = 1:length(fixed_methods_list)
    method = fixed_methods_list{m};
    scatter(const_torque_fixed_data.(method).cpu_times, const_torque_fixed_data.(method).errors, ...
        100, colors(m, :), 'filled', 'DisplayName', method, 'LineWidth', 1.5);
end

% Variable time step methods
var_methods_list = fieldnames(const_torque_var_data);
for m = 1:length(var_methods_list)
    method = var_methods_list{m};
    scatter(const_torque_var_data.(method).cpu_times, const_torque_var_data.(method).errors, ...
        100, colors(m+2, :), 'o', 'DisplayName', method, 'LineWidth', 1.5);
end

xlabel('CPU Time (s)');
ylabel('Max Error (rad/s)');
title('Maximum Simulation Error vs CPU Time');
legend('Location', 'best');
saveas(gcf, 'assets/plot3_error_vs_cputime.png');
saveas(gcf, 'assets/plot3_error_vs_cputime.fig');
grid on;
set(gca, 'XScale', 'log', 'YScale', 'log');

% Plot 4: Contour plot of constant system eigenvalues
% System eigenvalue = -b/J (characterizes system dynamics)
figure('Name', 'Contour: System Eigenvalue');

% Collect data for constant torque cases with system eigenvalues
cpu_data = [];
err_data = [];
z_data = [];

keys = fieldnames(results.fixed);
for k = 1:length(keys)
    if ~isfield(results.fixed.(keys{k}), 'is_sinusoidal')
        J = results.fixed.(keys{k}).J;
        b = results.fixed.(keys{k}).b;
        eigenval = -b/J;
        cpu = results.fixed.(keys{k}).cpu_time;
        error = results.fixed.(keys{k}).max_error;
        
        cpu_data = [cpu_data; cpu];
        err_data = [err_data; error];
        z_data = [z_data; eigenval];
    end
end

% Add variable step data
keys = fieldnames(results.variable);
for k = 1:length(keys)
    if ~isfield(results.variable.(keys{k}), 'is_sinusoidal')
        J = results.variable.(keys{k}).J;
        b = results.variable.(keys{k}).b;
        eigenval = -b/J;
        cpu = results.variable.(keys{k}).cpu_time;
        error = results.variable.(keys{k}).max_error;
        
        cpu_data = [cpu_data; cpu];
        err_data = [err_data; error];
        z_data = [z_data; eigenval];
    end
end

% Create contour plot using Dalys' style
% Remove inf/nan rows for gridding
good = isfinite(cpu_data) & isfinite(err_data) & isfinite(z_data) & err_data > 0 & cpu_data > 0;
cpu_data = cpu_data(good);
err_data = err_data(good);
z_data = z_data(good);

if ~isempty(cpu_data) && length(cpu_data) > 3
    % Build grid in log space (better spread)
    x = log10(cpu_data);
    y = log10(err_data);
    xi = linspace(min(x), max(x), 50);
    yi = linspace(min(y), max(y), 50);
    [X, Y] = meshgrid(xi, yi);
    
    Z = griddata(x, y, z_data, X, Y, "natural");
    
    contourf(10.^X, 10.^Y, Z, 12);
    set(gca, "XScale", "log", "YScale", "log");
    colorbar;
    ylabel(colorbar, 'Eigenvalue \lambda = -b/J');
    grid on;
    xlabel('CPU time (s)');
    ylabel('Max error (rad/s)');
    title('Contour: System Eigenvalue over (CPU, Error)');
    saveas(gcf, 'assets/plot4_contour_eigenvalue.png');
    saveas(gcf, 'assets/plot4_contour_eigenvalue.fig');
end

% Plot 5: Contour plot of input frequencies (sinusoidal cases only)
figure('Name', 'Contour: Input Frequency');

% Collect data for sinusoidal cases
cpu_data_f = [];
err_data_f = [];
z_data_f = [];

keys = fieldnames(results.fixed);
for k = 1:length(keys)
    if isfield(results.fixed.(keys{k}), 'is_sinusoidal')
        freq = results.fixed.(keys{k}).omega_freq;
        cpu = results.fixed.(keys{k}).cpu_time;
        
        % For sinusoidal inputs, calculate error vs RK4 reference (dt=0.001)
        % Find reference solution
        ref_key = sprintf('J%.2f_b%.1f_sine_w0%.1f_freq%.2f_dt%.4f_RK4', ...
            results.fixed.(keys{k}).J, results.fixed.(keys{k}).b, ...
            results.fixed.(keys{k}).omega0, freq, 0.001);
        ref_key = strrep(ref_key, '.', '_');
        
        if isfield(results.fixed, ref_key)
            omega_ref = results.fixed.(ref_key).omega;
            omega_curr = results.fixed.(keys{k}).omega;
            
            % Interpolate to same time points if needed
            if length(omega_curr) ~= length(omega_ref)
                t_curr = results.fixed.(keys{k}).t;
                t_ref = results.fixed.(ref_key).t;
                omega_curr_interp = interp1(t_curr, omega_curr, t_ref, 'linear', 'extrap');
                error = max(abs(omega_curr_interp - omega_ref));
            else
                error = max(abs(omega_curr - omega_ref));
            end
            
            cpu_data_f = [cpu_data_f; cpu];
            err_data_f = [err_data_f; error];
            z_data_f = [z_data_f; freq];
        end
    end
end

% Add variable step sinusoidal data
keys = fieldnames(results.variable);
for k = 1:length(keys)
    if isfield(results.variable.(keys{k}), 'is_sinusoidal')
        freq = results.variable.(keys{k}).omega_freq;
        cpu = results.variable.(keys{k}).cpu_time;
        
        % Find RK4 reference solution
        ref_key = sprintf('J%.2f_b%.1f_sine_w0%.1f_freq%.2f_dt%.4f_RK4', ...
            results.variable.(keys{k}).J, results.variable.(keys{k}).b, ...
            results.variable.(keys{k}).omega0, freq, 0.001);
        ref_key = strrep(ref_key, '.', '_');
        
        if isfield(results.fixed, ref_key)
            omega_ref = results.fixed.(ref_key).omega;
            omega_curr = results.variable.(keys{k}).omega;
            
            % Interpolate to same time points
            t_curr = results.variable.(keys{k}).t;
            t_ref = results.fixed.(ref_key).t;
            omega_curr_interp = interp1(t_curr, omega_curr, t_ref, 'linear', 'extrap');
            error = max(abs(omega_curr_interp - omega_ref));
            
            cpu_data_f = [cpu_data_f; cpu];
            err_data_f = [err_data_f; error];
            z_data_f = [z_data_f; freq];
        end
    end
end

% Create contour plot using Dalys' style
% Remove inf/nan rows for gridding
good_f = isfinite(cpu_data_f) & isfinite(err_data_f) & isfinite(z_data_f) & err_data_f > 0 & cpu_data_f > 0;
cpu_data_f = cpu_data_f(good_f);
err_data_f = err_data_f(good_f);
z_data_f = z_data_f(good_f);

if ~isempty(cpu_data_f) && length(cpu_data_f) > 3
    % Build grid in log space (better spread)
    x_f = log10(cpu_data_f);
    y_f = log10(err_data_f);
    xi_f = linspace(min(x_f), max(x_f), 50);
    yi_f = linspace(min(y_f), max(y_f), 50);
    [X_f, Y_f] = meshgrid(xi_f, yi_f);
    
    Z_f = griddata(x_f, y_f, z_data_f, X_f, Y_f, "natural");
    
    contourf(10.^X_f, 10.^Y_f, Z_f, 12);
    set(gca, "XScale", "log", "YScale", "log");
    colorbar;
    ylabel(colorbar, 'Input frequency w (rad/s)');
    grid on;
    xlabel('CPU time (s)');
    ylabel('Max error (rad/s)');
    title('Contour: Input Frequency over (CPU, Error)');
    saveas(gcf, 'assets/plot5_contour_frequency.png');
    saveas(gcf, 'assets/plot5_contour_frequency.fig');
end

fprintf('Plots generated successfully!\n');
fprintf('Plots saved to assets/ folder\n');

%% Additional Plots: CPU Time Comparisons and Simulation Results

% Plot 6: CPU Time Comparison for all methods (constant torque)
figure('Name', 'CPU Time Comparison');
bar_data = [];
bar_labels = {};

fixed_method_names = fieldnames(const_torque_fixed_data);
for m = 1:length(fixed_method_names)
    method = fixed_method_names{m};
    bar_data = [bar_data; mean(const_torque_fixed_data.(method).cpu_times)];
    bar_labels{m} = method;
end

var_method_names = fieldnames(const_torque_var_data);
for m = 1:length(var_method_names)
    method = var_method_names{m};
    bar_data = [bar_data; mean(const_torque_var_data.(method).cpu_times)];
    bar_labels{length(fixed_method_names) + m} = method;
end

bar(bar_data);
set(gca, 'XTickLabel', bar_labels);
ylabel('Mean CPU Time (s)');
title('Average CPU Time by Integration Method');
set(gca, 'YScale', 'log');
grid on;
saveas(gcf, 'assets/plot6_cpu_comparison.png');
saveas(gcf, 'assets/plot6_cpu_comparison.fig');

% Plot 7: Sample Simulation Results - Constant Torque Case
figure('Name', 'Sample Constant Torque Response');
keys = fieldnames(results.fixed);
% Find a representative case: J=100, b=10, tau=100, w0=10, dt=0.001, Euler
sample_key = '';
for k = 1:length(keys)
    if isfield(results.fixed.(keys{k}), 'tau_const') && ...
            results.fixed.(keys{k}).J == 100 && results.fixed.(keys{k}).b == 10 && ...
            results.fixed.(keys{k}).tau_const == 100 && results.fixed.(keys{k}).omega0 == 10 && ...
            results.fixed.(keys{k}).dt == 0.001 && strcmp(results.fixed.(keys{k}).method, 'Euler')
        sample_key = keys{k};
        break;
    end
end

if ~isempty(sample_key)
    subplot(2,2,1);
    plot(results.fixed.(sample_key).t, results.fixed.(sample_key).omega, 'LineWidth', 2);
    xlabel('Time (s)'); ylabel('Angular Velocity (rad/s)');
    title('Angular Velocity Response');
    grid on;
    
    subplot(2,2,2);
    plot(results.fixed.(sample_key).t, results.fixed.(sample_key).tau, 'LineWidth', 2);
    xlabel('Time (s)'); ylabel('Torque (N·m)');
    title('Applied Torque');
    grid on;
    
    subplot(2,2,3);
    plot(results.fixed.(sample_key).t, results.fixed.(sample_key).f_damping, 'LineWidth', 2);
    xlabel('Time (s)'); ylabel('Damping Force (N·m)');
    title('Damping Torque');
    grid on;
    
    subplot(2,2,4);
    plot(results.fixed.(sample_key).t, results.fixed.(sample_key).omega, 'b-', 'LineWidth', 2);
    hold on;
    % Calculate and plot theoretical
    if results.fixed.(sample_key).tau_const == 0
        omega_theory = results.fixed.(sample_key).omega0 * exp(-results.fixed.(sample_key).b/results.fixed.(sample_key).J * results.fixed.(sample_key).t);
    else
        omega_theory = results.fixed.(sample_key).omega0 * exp(-results.fixed.(sample_key).b/results.fixed.(sample_key).J * results.fixed.(sample_key).t) + ...
                       (results.fixed.(sample_key).tau_const/results.fixed.(sample_key).b) * (1 - exp(-results.fixed.(sample_key).b/results.fixed.(sample_key).J * results.fixed.(sample_key).t));
    end
    plot(results.fixed.(sample_key).t, omega_theory, 'r--', 'LineWidth', 2);
    xlabel('Time (s)'); ylabel('Angular Velocity (rad/s)');
    title('Simulated vs Theoretical');
    legend('Simulated (Euler)', 'Theoretical');
    grid on;
    
    sgtitle('Sample Simulation: Constant Torque (J=100, b=10, tau=100, w0=10, dt=0.001)');
end
saveas(gcf, 'assets/plot7_sample_constant_torque.png');
saveas(gcf, 'assets/plot7_sample_constant_torque.fig');

% Plot 8: Sample Simulation Results - Sinusoidal Torque Case
figure('Name', 'Sample Sinusoidal Torque Response');
sample_sin_key = '';
for k = 1:length(keys)
    if isfield(results.fixed.(keys{k}), 'is_sinusoidal') && results.fixed.(keys{k}).is_sinusoidal && ...
            results.fixed.(keys{k}).J == 100 && results.fixed.(keys{k}).b == 10 && ...
            results.fixed.(keys{k}).omega_freq == 0.1 && results.fixed.(keys{k}).omega0 == 10 && ...
            results.fixed.(keys{k}).dt == 0.001 && strcmp(results.fixed.(keys{k}).method, 'RK4')
        sample_sin_key = keys{k};
        break;
    end
end

if ~isempty(sample_sin_key)
    subplot(2,2,1);
    plot(results.fixed.(sample_sin_key).t, results.fixed.(sample_sin_key).omega, 'LineWidth', 2);
    xlabel('Time (s)'); ylabel('Angular Velocity (rad/s)');
    title('Angular Velocity Response');
    grid on;
    
    subplot(2,2,2);
    plot(results.fixed.(sample_sin_key).t, results.fixed.(sample_sin_key).tau, 'LineWidth', 2);
    xlabel('Time (s)'); ylabel('Torque (N·m)');
    title('Applied Sinusoidal Torque');
    grid on;
    
    subplot(2,2,3);
    plot(results.fixed.(sample_sin_key).t, results.fixed.(sample_sin_key).f_damping, 'LineWidth', 2);
    xlabel('Time (s)'); ylabel('Damping Force (N·m)');
    title('Damping Torque');
    grid on;
    
    subplot(2,2,4);
    plot(results.fixed.(sample_sin_key).t, results.fixed.(sample_sin_key).omega, 'LineWidth', 2);
    xlabel('Time (s)'); ylabel('Angular Velocity (rad/s)');
    title('Angular Velocity Over Time');
    grid on;
    
    sgtitle('Sample Simulation: Sinusoidal Torque (J=100, b=10, freq=0.1 rad/s, w0=10, dt=0.001, RK4)');
end
saveas(gcf, 'assets/plot8_sample_sinusoidal_torque.png');
saveas(gcf, 'assets/plot8_sample_sinusoidal_torque.fig');

% Plot 9: Method Comparison for Same Case
figure('Name', 'Integration Method Comparison');
case_J = 100; case_b = 10; case_tau = 100; case_w0 = 10; case_dt = 0.001;

keys_to_plot = {};
for k = 1:length(keys)
    if isfield(results.fixed.(keys{k}), 'tau_const') && ...
            results.fixed.(keys{k}).J == case_J && results.fixed.(keys{k}).b == case_b && ...
            results.fixed.(keys{k}).tau_const == case_tau && results.fixed.(keys{k}).omega0 == case_w0 && ...
            results.fixed.(keys{k}).dt == case_dt
        keys_to_plot = [keys_to_plot; keys{k}];
    end
end

hold on;
colors = {'b-', 'r-', 'g-', 'k-'};
for i = 1:length(keys_to_plot)
    plot(results.fixed.(keys_to_plot{i}).t, results.fixed.(keys_to_plot{i}).omega, colors{i}, 'LineWidth', 2, ...
         'DisplayName', results.fixed.(keys_to_plot{i}).method);
end

xlabel('Time (s)');
ylabel('Angular Velocity (rad/s)');
title(sprintf('Integration Method Comparison (J=%g, b=%g, tau=%g, w0=%g, dt=%g)', ...
    case_J, case_b, case_tau, case_w0, case_dt));
legend('Location', 'best');
grid on;
saveas(gcf, 'assets/plot9_method_comparison.png');
saveas(gcf, 'assets/plot9_method_comparison.fig');

% Plot 10: Error vs J and b (System Properties)
figure('Name', 'Error by System Properties');
J_errors = [];
b_errors = [];
errors = [];

keys = fieldnames(results.fixed);
for k = 1:length(keys)
    if ~isfield(results.fixed.(keys{k}), 'is_sinusoidal') && isfield(results.fixed.(keys{k}), 'max_error')
        error = results.fixed.(keys{k}).max_error;
        if ~isinf(error) && ~isnan(error)
            J_errors = [J_errors; results.fixed.(keys{k}).J];
            b_errors = [b_errors; results.fixed.(keys{k}).b];
            errors = [errors; error];
        end
    end
end

scatter(J_errors, b_errors, 100, errors, 'filled');
colorbar;
xlabel('Moment of Inertia J (kg·m²)');
ylabel('Damping Coefficient b (N·m·s/rad)');
title('Max Simulation Error by System Properties');
colormap('jet');
set(gca, 'XScale', 'log', 'YScale', 'log');
grid on;
saveas(gcf, 'assets/plot10_error_by_properties.png');
saveas(gcf, 'assets/plot10_error_by_properties.fig');

fprintf('Additional plots generated and saved!\n');
