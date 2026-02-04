%% Project 1 - Part 2 (Option 1 only) : run + plots
clear; clc;

mdl = "Project1_part2_opt1_spring";
load_system(mdl);
set_param(mdl,"StopTime","25");

%% ---- Required parameters (per Part 2 prompt) ----
J1 = 100; b1 = 1;
J2 = 1;   b2 = 1;

A_list = [1 100];          % step input magnitudes
k_list = [10 100 1000];    % stiffness values

% Initial conditions (only if your Integrators reference these variables)
omega10 = 0; theta10 = 0;
omega20 = 0; theta20 = 0;

assignin("base","J1",J1); assignin("base","b1",b1);
assignin("base","J2",J2); assignin("base","b2",b2);
assignin("base","omega10",omega10); assignin("base","theta10",theta10);
assignin("base","omega20",omega20); assignin("base","theta20",theta20);

%% ---- Solver configurations required ----
cfgs = struct( ...
    "name",   {"ode1 dt=0.1","ode1 dt=1","ode4 dt=0.1","ode4 dt=1","ode45"}, ...
    "solver", {"ode1","ode1","ode4","ode4","ode45"}, ...
    "isFixed",{ true,  true,  true,  true,  false}, ...
    "dt",     { 0.1,   1.0,   0.1,   1.0,   NaN } );

%% ---- Helper to safely pull signals ----
getsig = @(simOut, nm) local_get(simOut, nm);

%% ---- Optional warm-up (reduces first-run overhead) ----
assignin("base","A",A_list(1));
assignin("base","k",k_list(1));
set_param(mdl,"SolverType","Fixed-step","Solver","ode4","FixedStep","0.1");
sim(mdl);

%% ---- Run sweep and store results ----
rows = struct([]);
r = 0;

fprintf("Running Option 1 sweep...\n");

for A = A_list
  for k = k_list

    % push parameters used by the model
    assignin("base","A",A);
    assignin("base","k",k);

    for c = 1:numel(cfgs)

      % set solver
      if cfgs(c).isFixed
        set_param(mdl, ...
          "SolverType","Fixed-step", ...
          "Solver",cfgs(c).solver, ...
          "FixedStep",num2str(cfgs(c).dt));
      else
        set_param(mdl, ...
          "SolverType","Variable-step", ...
          "Solver",cfgs(c).solver);
      end

      % run
      tStart = tic;
      simOut = sim(mdl,"ReturnWorkspaceOutputs","on");
      cpu_s  = toc(tStart);

      % REQUIRED outputs (must match your To Workspace names)
      w1 = getsig(simOut,"omega1_out");
      w2 = getsig(simOut,"omega2_out");

      % OPTIONAL diagnostics (only if you created To Workspace blocks)
      theta_rel = local_tryget(simOut,"theta_rel_out");
      tau_k     = local_tryget(simOut,"tau_k_out");
      tau_net2  = local_tryget(simOut,"tau_net2_out");

      r = r + 1;
      rows(r).A = A;
      rows(r).k = k;
      rows(r).solver = string(cfgs(c).solver);
      rows(r).cfgName = string(cfgs(c).name);
      rows(r).dt = cfgs(c).dt;
      rows(r).cpu_s = cpu_s;

      rows(r).w1 = w1;
      rows(r).w2 = w2;

      rows(r).theta_rel = theta_rel;
      rows(r).tau_k     = tau_k;
      rows(r).tau_net2  = tau_net2;
    end
  end
end

resultsTbl = struct2table(rows);
save("part2_option1_results.mat","resultsTbl");
disp("Saved: part2_option1_results.mat");

%% =========================
%  REQUIRED PLOTS (Option 1)
%% =========================

% 1) Shaft speed vs time, compare stiffnesses (k) for each A, for EACH solver config
for A = A_list
  for c = 1:numel(cfgs)

    figure("Name", sprintf("Option 1: omega vs time | A=%g | %s", A, cfgs(c).name));
    hold on; grid on;

    for k = k_list
      m = resultsTbl.A==A & resultsTbl.k==k & resultsTbl.cfgName==string(cfgs(c).name);
      i = find(m,1);
      if isempty(i), continue; end

      t = resultsTbl.w1(i).Time;

      plot(t, resultsTbl.w1(i).Data, ...
        "DisplayName", sprintf("\\omega_1 (k=%g)", k));
      plot(t, resultsTbl.w2(i).Data, "--", ...
        "DisplayName", sprintf("\\omega_2 (k=%g)", k));
    end

    xlabel("Time (s)");
    ylabel("\omega (rad/s)");
    title(sprintf("Option 1 (spring coupling): A=%g, %s", A, cfgs(c).name));
    legend("Location","best");
  end
end

% 2) CPU time table (this satisfies the "tabulate CPU time" part for Option 1)
cpuTable = resultsTbl(:,["A","k","cfgName","solver","dt","cpu_s"]);
disp("CPU time table (all runs):");
disp(cpuTable);

% A compact summary (median CPU for each config)
cpuSummary = groupsummary(cpuTable, ["A","k","cfgName"], "median", "cpu_s");
disp("CPU time summary (median cpu_s):");
disp(cpuSummary);

%% ===== Local functions =====
function ts = local_get(simOut, varName)
  % Works if To Workspace saved variable appears in simOut
  try
    ts = simOut.get(varName);
  catch
    error("Couldn't find '%s' in simOut. Check the To Workspace variable name.", varName);
  end
end

function ts = local_tryget(simOut, varName)
  try
    ts = simOut.get(varName);
  catch
    ts = []; % not present, that's ok
  end
end
