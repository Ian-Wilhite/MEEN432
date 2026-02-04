clear;
clc;
%{
function [x, y, failed] = runSim(mdl)
    try
        output = sim(mdl);
        x = output.scopeData1(:, 1);
        y = output.ScopeData1(:, 2);
        failed = false;
    catch
        x = [];
        y = [];
        failed = true;
    end
end
%}

%---------PARAMTERS---------%
w_0 = [10, 0]; % Initial Conditionas, rad/s
J = [100, 0.01]; % Rotational inertia, kg-m^2
b = [10, 0.1]; % Damping coefficient, N-m-s/rad

% Applied torque
A = [0, 100]; % Constant, N-m
w = [0.1, 100]; % Sinusoidal Freqeuncy Conditions rad/s

% Integration Methods
    %Fixed Time Step
dt = [0.001, 0.1, 1]; % Seconds
solvers = {'ode1', 'ode4'};
solver_name = {'Euler', 'Runge Kutta'};

    %Variable Time Step
var_solvers = {'ode45', 'ode23tb'};

% Loading Simulink Model
load_system('Part1_sim');


%------------------FIXED TIME STEPS--------------------%
%%{
for s = 1: length(solvers)
    for t = 1: length(dt)
        % Setting solver parameters
        set_param('Part1_sim', ...
            'SolverType', 'Fixed-step', ...
            'Solver', solvers{s}, ...
            'FixedStep', num2str(dt(t)));

        for ini = 1: length(w_0) %initial conditions loop
            w_sim = w_0(ini);
            for rot = 1: length(J) %rotational inertia loop
                J_sim = J(rot);
                for dam = 1: length(b) % damping cefficient loop
                    b_sim = b(dam);
        
                %-------CONSTANT TORQUE SIMS--------%
                    for con = 1: length(A)
                        A_sim = A(con);
        
                        set_param('Part1_sim/Constant', 'Value', num2str(A_sim));
                        set_param('Part1_sim/Integrator2', 'InitialCondition', num2str(w_sim));
                        set_param('Part1_sim/Gain2', 'Gain', num2str(J_sim));
                        set_param('Part1_sim/Gain3', 'Gain', num2str(b_sim));
                        
                        try
                            output = sim("Part1_sim");
                        
                            x = output.ScopeData1(:, 1);
                            y = output.ScopeData1(:, 2);
                        
                            sim_failed = false;
                        catch ME
                            sim_failed = true;
                            fprintf("FAILED SIM\n");
                        end
%{
                        if ~sim_failed
                            figure;
                            plot(x, y);
                            title({['Integration Method: ', solver_name{s}, ' with dt = ', num2str(dt(t))]; ...
                                ['Constant (A = ', num2str(A_sim), '): w0 = ', num2str(w_sim), ', J = ', num2str(J_sim), ', b = ', num2str(b_sim)]});
                            xlabel('Time (s)');
                            ylabel('Amplitude');
                        end
%}
                    end
        
                %-------Sinusoidal Sims--------%
                    for sin = 1: length(w)
                        f_sim = w(sin);
                        
                        set_param('Part1_sim/Sine Wave', 'Frequency', num2str(f_sim));
                        set_param('Part1_sim/Integrator', 'InitialCondition', num2str(w_sim));
                        set_param('Part1_sim/Gain', 'Gain', num2str(J_sim));
                        set_param('Part1_sim/Gain1', 'Gain', num2str(b_sim));
                        
                        try
                            output = sim("Part1_sim");
                        
                            x = output.ScopeData1(:, 1);
                            y = output.ScopeData1(:, 2);
                        
                            sim_failed = false;
                        catch ME
                            sim_failed = true;
                            fprintf("FAILED SIM\n");
                        end
%{
                        if ~sim_failed
                            figure;
                            plot(x, y);
                            title({['Integration Method: ', solver_name{s}, ' with dt = ', num2str(dt(t))];['Sinusoid (f = ', num2str(f_sim), '): w0 = ', num2str(w_sim), ', J = ', num2str(J_sim), ', b = ', num2str(b_sim)]});
                            xlabel('Time (s)');
                            ylabel('Amplitude');
                        end 
%}
                    end
                end
            end
        end
    end
end
%}

%----------------VARIABLE TIME STEPS-----------------%
for v = 1: length(var_solvers)
    %setting solver parameters
    set_param('Part1_sim', ...
        'SolverType', 'Variable-step', ...
        'Solver', var_solvers{v});

    for ini = 1: length(w_0) %initial conditions loop
        w_sim = w_0(ini);
        for rot = 1: length(J) %rotational inertia loop
            J_sim = J(rot);
            for dam = 1: length(b) % damping cefficient loop
                b_sim = b(dam);
    
            %-------Constant Sims--------%
                for con = 1: length(A)
                    A_sim = A(con);
    
                    set_param('Part1_sim/Constant', 'Value', num2str(A_sim));
                    set_param('Part1_sim/Integrator2', 'InitialCondition', num2str(w_sim));
                    set_param('Part1_sim/Gain2', 'Gain', num2str(J_sim));
                    set_param('Part1_sim/Gain3', 'Gain', num2str(b_sim));
                    
                    output = sim("Part1_sim");
                    x = output.ScopeData1(:, 1);
                    y = output.ScopeData1(:, 2);

%{
                    figure;
                    plot(x, y);
                    title({['Integration Method: ', var_solvers{v}]; ['Constant (A = ', num2str(A_sim), '): w0 = ', num2str(w_sim), ', J = ', num2str(J_sim), ', b = ', num2str(b_sim)]});
                    xlabel('Time (s)');
                    ylabel('Amplitude');
%}
                end
    
            %-------Sinusoidal Sims--------%
                for sin = 1: length(w)
                    f_sim = w(sin);
                    
                    set_param('Part1_sim/Sine Wave', 'Frequency', num2str(f_sim));
                    set_param('Part1_sim/Integrator', 'InitialCondition', num2str(w_sim));
                    set_param('Part1_sim/Gain', 'Gain', num2str(J_sim));
                    set_param('Part1_sim/Gain1', 'Gain', num2str(b_sim));
                    
                    output = sim("Part1_sim");
                    x = output.ScopeData(:, 1);
                    y = output.ScopeData(:, 2);
                
%{
                    figure;
                    plot(x, y);
                    title({['Integration Method: ', var_solvers{v}];['Sinusoid (f = ', num2str(f_sim), '): w0 = ', num2str(w_sim), ', J = ', num2str(J_sim), ', b = ', num2str(b_sim)]});
                    xlabel('Time (s)');
                    ylabel('Amplitude');
%}
                end
            end
        end
    end
end


