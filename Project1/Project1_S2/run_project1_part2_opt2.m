%% Project 1 - Part 2 (Option 2: combined inertia/damping) runner
clear; clc;

mdl = "Project1_part2_opt2_combined";
load_system(mdl);
set_param(mdl, "StopTime", "25");

% Parameters
J1 = 100; b1 = 1;
J2 = 1;   b2 = 1;

Jeq = J1 + J2;
beq = b1 + b2;

% Initial conditions
omega0 = 0;
theta0 = 0;

% Input step amplitudes required
A_list = [1 100];

% Push parameters to base workspace (model should reference these)
assignin("base","J1",J1); assignin("base","b1",b1);
assignin("base","J2",J2); assignin("base","b2",b2);
assignin("base","Jeq",Jeq);
assignin("base","beq",beq);
assignin("base","omega0",omega0);
assignin("base","theta0",theta0);

% Solver configurations required
cfgs = struct( ...
    "solver",  {"ode1","ode1","ode4","ode4","ode45"}, ...
    "isFixed", { true,  true,  true,  true,  false}, ...
    "dt",      { 0.1,   1.0,   0.1,   1.0,   NaN } );

rows = struct([]);
r = 0;

% Warm-up run
assignin("base","A",A_list(1));
set_param(mdl, "SolverType","Fixed-step", "Solver","ode4", "FixedStep","0.1");
sim(mdl);

fprintf("Running Part 2 Option 2 sweep...\n");

for A = A_list
    assignin("base","A",A);

    for c = 1:numel(cfgs)
        solver  = cfgs(c).solver;
        isFixed = cfgs(c).isFixed;
        dt      = cfgs(c).dt;

        if isFixed
            set_param(mdl, ...
                "SolverType","Fixed-step", ...
                "Solver", solver, ...
                "FixedStep", num2str(dt));
        else
            set_param(mdl, ...
                "SolverType","Variable-step", ...
                "Solver", solver);
        end

        tStart = tic;
        simOut = sim(mdl, "ReturnWorkspaceOutputs","on");
        cpu_s  = toc(tStart);

        % These MUST match your To Workspace variable names:
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

% CPU time table
disp("CPU time summary (Option 2):");
disp(resultsTbl(:, {'A','solver','dt','cpu_s'}));

% Plots (one figure per A; compare all solvers on same plot)
for A = A_list
    figure("Name", sprintf("Option 2: omega vs time (A=%g)", A));
    hold on; grid on;

    m = resultsTbl.A==A;
    idx = find(m);

    for ii = idx'
        t = resultsTbl.omega(ii).Time;
        y = resultsTbl.omega(ii).Data;

        if isfinite(resultsTbl.dt(ii))
            lbl = sprintf("%s dt=%g", resultsTbl.solver(ii), resultsTbl.dt(ii));
        else
            lbl = sprintf("%s", resultsTbl.solver(ii));
        end

        plot(t, y, "DisplayName", lbl);
    end

    xlabel("Time (s)");
    ylabel("\omega (rad/s)");
    title(sprintf("Option 2 (combined), A=%g", A));
    legend("Location","best");
end

save("part2_option2_results.mat","resultsTbl");
fprintf("Saved results to part2_option2_results.mat\n");
%% Project 1 - Part 2 (Option 2: combined inertia/damping) runner
clear; clc;

mdl = "Project1_part2_opt2_combined";
load_system(mdl);
set_param(mdl, "StopTime", "25");

% Parameters
J1 = 100; b1 = 1;
J2 = 1;   b2 = 1;

Jeq = J1 + J2;
beq = b1 + b2;

% Initial conditions
omega0 = 0;
theta0 = 0;

% Input step amplitudes required
A_list = [1 100];

% Push parameters to base workspace (model should reference these)
assignin("base","J1",J1); assignin("base","b1",b1);
assignin("base","J2",J2); assignin("base","b2",b2);
assignin("base","Jeq",Jeq);
assignin("base","beq",beq);
assignin("base","omega0",omega0);
assignin("base","theta0",theta0);

% Solver configurations required
cfgs = struct( ...
    "solver",  {"ode1","ode1","ode4","ode4","ode45"}, ...
    "isFixed", { true,  true,  true,  true,  false}, ...
    "dt",      { 0.1,   1.0,   0.1,   1.0,   NaN } );

rows = struct([]);
r = 0;

% Warm-up run
assignin("base","A",A_list(1));
set_param(mdl, "SolverType","Fixed-step", "Solver","ode4", "FixedStep","0.1");
sim(mdl);

fprintf("Running Part 2 Option 2 sweep...\n");

for A = A_list
    assignin("base","A",A);

    for c = 1:numel(cfgs)
        solver  = cfgs(c).solver;
        isFixed = cfgs(c).isFixed;
        dt      = cfgs(c).dt;

        if isFixed
            set_param(mdl, ...
                "SolverType","Fixed-step", ...
                "Solver", solver, ...
                "FixedStep", num2str(dt));
        else
            set_param(mdl, ...
                "SolverType","Variable-step", ...
                "Solver", solver);
        end

        tStart = tic;
        simOut = sim(mdl, "ReturnWorkspaceOutputs","on");
        cpu_s  = toc(tStart);

        % These MUST match your To Workspace variable names:
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

% CPU time table
disp("CPU time summary (Option 2):");
disp(resultsTbl(:, {'A','solver','dt','cpu_s'}));

% Plots (one figure per A; compare all solvers on same plot)
for A = A_list
    figure("Name", sprintf("Option 2: omega vs time (A=%g)", A));
    hold on; grid on;

    m = resultsTbl.A==A;
    idx = find(m);

    for ii = idx'
        t = resultsTbl.omega(ii).Time;
        y = resultsTbl.omega(ii).Data;

        if isfinite(resultsTbl.dt(ii))
            lbl = sprintf("%s dt=%g", resultsTbl.solver(ii), resultsTbl.dt(ii));
        else
            lbl = sprintf("%s", resultsTbl.solver(ii));
        end

        plot(t, y, "DisplayName", lbl);
    end

    xlabel("Time (s)");
    ylabel("\omega (rad/s)");
    title(sprintf("Option 2 (combined), A=%g", A));
    legend("Location","best");
end

save("part2_option2_results.mat","resultsTbl");
fprintf("Saved results to part2_option2_results.mat\n");
