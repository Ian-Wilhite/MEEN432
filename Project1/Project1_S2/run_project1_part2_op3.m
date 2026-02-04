%% Project 1 - Part 2 (Option 3: S2 integrates combined rigid model)
clear; clc;

mdl = "Project1_part2_opt3_pass_to_s2";   
load_system(mdl);
set_param(mdl, "StopTime", "25");

% Parameters
J1 = 100; b1 = 1;
J2 = 1;   b2 = 1;

omega0 = 0;
theta0 = 0;

A_list = [1 100];

% Push to base workspace
assignin("base","J1",J1); assignin("base","b1",b1);
assignin("base","J2",J2); assignin("base","b2",b2);
assignin("base","omega0",omega0);
assignin("base","theta0",theta0);

cfgs = struct( ...
    "solver",  {"ode1","ode1","ode4","ode4","ode45"}, ...
    "isFixed", { true,  true,  true,  true,  false}, ...
    "dt",      { 0.1,   1.0,   0.1,   1.0,   NaN } );

rows = struct([]);
r = 0;

fprintf("Running Part 2 Option 3 sweep...\n");

for A = A_list
    assignin("base","A",A);

    for c = 1:numel(cfgs)
        solver  = cfgs(c).solver;
        isFixed = cfgs(c).isFixed;
        dt      = cfgs(c).dt;

        if isFixed
            set_param(mdl, "SolverType","Fixed-step", "Solver",solver, "FixedStep",num2str(dt));
        else
            set_param(mdl, "SolverType","Variable-step", "Solver",solver);
        end

        tStart = tic;
        simOut = sim(mdl, "ReturnWorkspaceOutputs","on");
        cpu_s  = toc(tStart);

        % These MUST match your To Workspace block names
        w  = simOut.get("omega_out");
        th = simOut.get("theta_out");

        r = r + 1;
        rows(r).A = A;
        rows(r).solver = string(solver);
        rows(r).dt = dt;
        rows(r).cpu_s = cpu_s;
        rows(r).omega = w;
        rows(r).theta = th;
    end
end

resultsTbl = struct2table(rows);

disp("CPU time summary (Option 3):");
disp(resultsTbl(:,["A","solver","dt","cpu_s"]));

%% Plots: omega vs time (one plot per A, compare solvers)
for A = A_list
    figure("Name", sprintf("Option 3 (S2 integrates), A=%g", A));
    hold on; grid on;

    m = resultsTbl.A==A;

    for i = find(m)'
        t = resultsTbl.omega(i).Time;
        y = resultsTbl.omega(i).Data;

        if resultsTbl.solver(i) == "ode45"
            lab = "ode45";
        else
            lab = sprintf("%s dt=%g", resultsTbl.solver(i), resultsTbl.dt(i));
        end

        plot(t, y, "DisplayName", lab);
    end

    xlabel("Time (s)");
    ylabel("\omega (rad/s)");
    title(sprintf("Option 3 (inertia passed to S2): A=%g", A));
    legend("Location","best");
end

save("part2_option3_results.mat","resultsTbl");
fprintf("Saved results to part2_option3_results.mat\n");
