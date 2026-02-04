%% Project 1 - Part 1 runner (S1 only)
% Assumes To Workspace signals named:
%   out.omega_out, out.tau_net_out, out.tau_b_out   (Timeseries)
% Optional:
%   out.theta_out  (Timeseries)
%
% Produces:
%   - Results table in MATLAB workspace: resultsTbl
%   - Figures: error vs dt, cpu vs dt, error vs cpu, contour(eigs), contour(freq)
%
% Spec requirements: :contentReference[oaicite:3]{index=3} :contentReference[oaicite:4]{index=4}

clear; clc;


mdl = "Project1_s1";  
load_system(mdl);

%% ==== Part 1 parameter sweeps (from PDF) ====
omega0_list = [10.0, 0.0];        % rad/s
J_list      = [100, 0.01];        % kg-m^2
b_list      = [10, 0.1];          % N-m-s/rad
A_list      = [0.0, 100];         % N-m (constant or sine amplitude)
w_list      = [0.1, 100];         % rad/s (sine frequency)
Tstop       = 25;                 % sec

dt_list     = [0.001, 0.1, 1];    % sec fixed-step
fixedSolvers = ["ode1","ode4"];   % Euler, RK4
varSolvers   = ["ode45","ode23tb"];

%% ==== Helper: safe extraction of logged signals ====
get_ts = @(simOut, field) local_get_timeseries(simOut, field);

%% ==== Build a reference library for sinusoidal cases:
% Reference = ode4 with dt = 0.001, as required :contentReference[oaicite:5]{index=5}
refMap = containers.Map;  % key -> struct(t, omega)

fprintf("Building sine reference runs (ode4, dt=0.001)...\n");
for J = J_list
for b = b_list
for omega0 = omega0_list
for A = A_list
for w = w_list
    key = local_key(J,b,omega0,A,w);
    % Only meaningful sine if A>0 (still fine if A=0)
    simIn = Simulink.SimulationInput(mdl);
    simIn = simIn.setVariable("J",J).setVariable("b",b) ...
                 .setVariable("omega0",omega0).setVariable("A",A) ...
                 .setVariable("w",w).setVariable("useSine",1);
    simIn = simIn.setModelParameter( ...
        "StopTime", num2str(Tstop), ...
        "SolverType","Fixed-step", ...
        "Solver","ode4", ...
        "FixedStep","0.001" ...
    );

    try
        tStart = tic;
        simOut = sim(simIn);
        cpu = toc(tStart); %#ok<NASGU>
        omega_ts = get_ts(simOut, "omega_out");
        refMap(key) = struct("t", omega_ts.Time, "omega", omega_ts.Data);
    catch ME
        warning("Reference run failed for %s: %s", key, ME.message);
        refMap(key) = struct("t", [], "omega", []);
    end
end
end
end
end
end

%% ==== Main sweep: all methods, all cases ====
rows = struct([]);
fprintf("Running full sweep...\n");

% We'll run:
%   - Step/constant cases (useSine=0, A in [0 100], compare to analytic)
%   - Sinusoidal cases (useSine=1, compare to refMap)

case_id = 0;

for J = J_list
for b = b_list
for omega0 = omega0_list
for A = A_list

    %% --- CONSTANT TORQUE (step) ---
    case_id = case_id + 1;
    rows = [rows, local_run_case_set(mdl, Tstop, get_ts, ...
        case_id, J,b,omega0,A, NaN, 0, dt_list, fixedSolvers, varSolvers, refMap)]; %#ok<AGROW>

    %% --- SINUSOIDAL TORQUE (if you want both w values) ---
    for w = w_list
        case_id = case_id + 1;
        rows = [rows, local_run_case_set(mdl, Tstop, get_ts, ...
            case_id, J,b,omega0,A, w, 1, dt_list, fixedSolvers, varSolvers, refMap)]; %#ok<AGROW>
    end

end
end
end
end

resultsTbl = struct2table(rows);

%% ==== Plots required in Part 1 :contentReference[oaicite:6]{index=6} ====

% Separate step vs sine
isSine = resultsTbl.useSine == 1;
stepTbl = resultsTbl(~isSine,:);
sineTbl = resultsTbl(isSine,:);

%% 1) Max error vs dt (fixed-step only), for step inputs
figure("Name","Max Error vs dt (fixed-step, step inputs)");
hold on; grid on;

for s = fixedSolvers
    mask = stepTbl.isFixed & stepTbl.solver == s;

    [dts, errMax] = local_aggregate_by_dt(stepTbl(mask,:), dt_list);

    plot(dts, errMax, "-o", "DisplayName", solverLabel(s));
end

set(gca,"XScale","log","YScale","log");
xlabel("Fixed step dt (s)");
ylabel("Max |omega - omega_{theory}| (rad/s)");
legend("Location","best");
%% 2) CPU time vs dt (fixed-step only), step inputs
figure("Name","CPU Time vs dt (fixed-step, step inputs)");
hold on; grid on;

for s = fixedSolvers
    mask = stepTbl.isFixed & stepTbl.solver == s;

    [dts, cpuMean] = local_aggregate_cpu_by_dt(stepTbl(mask,:), dt_list);

    plot(dts, cpuMean, "-o", "DisplayName", solverLabel(s));
end

set(gca,"XScale","log","YScale","log");
xlabel("Fixed step dt (s)");
ylabel("CPU time (s)");
legend("Location","best");

%% 3) Max error vs CPU time (all solvers), step inputs
figure("Name","Max Error vs CPU (step inputs, all solvers)");
hold on; grid on;

solvers = unique(stepTbl.solver);
colors = lines(numel(solvers));

for i = 1:numel(solvers)
    s = solvers(i);
    mask = stepTbl.solver == s & isfinite(stepTbl.maxErr);

    scatter(stepTbl.cpu_s(mask), ...
            stepTbl.maxErr(mask), ...
            36, ...
            colors(i,:), ...
            "filled", ...
            "DisplayName", solverLabel(s));
end

set(gca,"XScale","log","YScale","log");
xlabel("CPU time (s)");
ylabel("Max error (rad/s)");
legend("Location","best");
title("Each point = one simulation run");


%% 4) Contour of constant eigenvalue lambda = -b/J over (CPU, Error) for step inputs
% (We grid scattered points using griddata)
figure("Name","Contour: eigenvalue over (CPU, Error) (step inputs)");
local_contour_cpu_err(stepTbl.cpu_s, stepTbl.maxErr, stepTbl.lambda, "Eigenvalue \lambda = -b/J");

%% 5) Contour of input frequency over (CPU, Error) for sine inputs
figure("Name","Contour: input frequency over (CPU, Error) (sine inputs)");
local_contour_cpu_err(sineTbl.cpu_s, sineTbl.maxErr, sineTbl.w, "Input frequency w (rad/s)");

fprintf("\nDONE. Results are in resultsTbl.\n");

%% ===================== LOCAL FUNCTIONS =====================

function rows = local_run_case_set(mdl, Tstop, get_ts, case_id, J, b, omega0, A, w, useSine, dt_list, fixedSolvers, varSolvers, refMap)

    key = local_key(J,b,omega0,A,w);

    % Build solver configs
    cfgs = struct([]);
    idx = 0;

    % Fixed-step
    for s = fixedSolvers
        for dt = dt_list
            idx = idx + 1;
            cfgs(idx).solver  = string(s);
            cfgs(idx).isFixed = true;
            cfgs(idx).dt      = dt;
        end
    end

    % Variable-step
    for s = varSolvers
        idx = idx + 1;
        cfgs(idx).solver  = string(s);
        cfgs(idx).isFixed = false;
        cfgs(idx).dt      = NaN;
    end

    % ---- Template struct (ALL fields defined once) ----
    template = struct( ...
        "case_id", 0, ...
        "J", 0, "b", 0, "omega0", 0, "A", 0, "w", NaN, "useSine", 0, ...
        "solver", "", "isFixed", false, "dt", NaN, ...
        "cpu_s", NaN, "maxErr", Inf, "lambda", NaN, ...
        "simOk", false ...
    );

    % Preallocate rows to same shape
    rows = repmat(template, 1, numel(cfgs));

    for k = 1:numel(cfgs)
        solver  = cfgs(k).solver;
        isFixed = cfgs(k).isFixed;
        dt      = cfgs(k).dt;

        simIn = Simulink.SimulationInput(mdl);
        simIn = simIn.setVariable("J",J).setVariable("b",b) ...
                     .setVariable("omega0",omega0).setVariable("A",A);

        if useSine == 1
            simIn = simIn.setVariable("w",w).setVariable("useSine",1);
        else
            simIn = simIn.setVariable("w",0.0).setVariable("useSine",0);
        end

        if isFixed
            simIn = simIn.setModelParameter( ...
                "StopTime", num2str(Tstop), ...
                "SolverType","Fixed-step", ...
                "Solver", char(solver), ...
                "FixedStep", num2str(dt) ...
            );
        else
            simIn = simIn.setModelParameter( ...
                "StopTime", num2str(Tstop), ...
                "SolverType","Variable-step", ...
                "Solver", char(solver) ...
            );
        end

        cpu_s = NaN;
        maxErr = Inf;
        simOk = true;

        try
            tStart = tic;
            simOut = sim(simIn);
            cpu_s = toc(tStart);

            omega_ts   = get_ts(simOut, "omega_out");
            tau_net_ts = get_ts(simOut, "tau_net_out");
            tau_b_ts   = get_ts(simOut, "tau_b_out"); %#ok<NASGU>

            t = omega_ts.Time;
            omega_sim = omega_ts.Data;

            if useSine == 0
                omega_th = local_omega_theory(t, omega0, A, b, J);
                maxErr = max(abs(omega_sim - omega_th));
            else
                ref = refMap(key);
                if isempty(ref.t)
                    maxErr = Inf;
                else
                    omega_ref = interp1(ref.t, ref.omega, t, "linear", "extrap");
                    maxErr = max(abs(omega_sim - omega_ref));
                end
            end

            if any(~isfinite(omega_sim)) || any(~isfinite(tau_net_ts.Data))
                maxErr = Inf;
                simOk = false;
            end

        catch
            simOk = false;
            maxErr = Inf;
        end

        lambda = -b / J;

        % Fill the preallocated row (fieldnames identical always)
        rows(k) = template;
        rows(k).case_id = case_id;

        rows(k).J = J; rows(k).b = b; rows(k).omega0 = omega0; rows(k).A = A;
        rows(k).w = w; rows(k).useSine = useSine;

        rows(k).solver = solver;
        rows(k).isFixed = isFixed;
        rows(k).dt = dt;

        rows(k).cpu_s = cpu_s;
        rows(k).maxErr = maxErr;
        rows(k).lambda = lambda;
        rows(k).simOk = simOk;
    end
end


function ts = local_get_timeseries(simOut, varName)
    % Works when To Workspace Variable name is omega_out, tau_net_out, etc.
    if isprop(simOut, varName)
        ts = simOut.(varName);
    else
        % Sometimes data is stored via get()
        try
            ts = simOut.get(varName);
        catch
            error("Could not find '%s' in SimulationOutput. Check To Workspace variable names.", varName);
        end
    end
end


function ts = local_try_get(get_ts, simOut, fieldName)
    try
        ts = get_ts(simOut, fieldName);
    catch
        ts = [];
    end
end

function key = local_key(J,b,omega0,A,w)
    if isnan(w); w = -999; end
    key = sprintf("J=%g|b=%g|w0=%g|A=%g|w=%g", J,b,omega0,A,w);
end

function omega_th = local_omega_theory(t, omega0, A, b, J)
    % For J*dw/dt + b*w = A (constant)
    if b == 0
        omega_th = omega0 + (A/J).*t;
    else
        omega_ss = A / b;
        omega_th = (omega0 - omega_ss).*exp(-(b/J).*t) + omega_ss;
    end
end

function [dts, errMax] = local_aggregate_by_dt(T, dt_list)
    dts = dt_list(:);
    errMax = nan(size(dts));
    for i = 1:numel(dts)
        mask = abs(T.dt - dts(i)) < 1e-12;
        errMax(i) = max(T.maxErr(mask));
    end
end

function [dts, cpuMean] = local_aggregate_cpu_by_dt(T, dt_list)
    dts = dt_list(:);
    cpuMean = nan(size(dts));
    for i = 1:numel(dts)
        mask = abs(T.dt - dts(i)) < 1e-12;
        cpuMean(i) = mean(T.cpu_s(mask), "omitnan");
    end
end

function local_contour_cpu_err(cpu, err, z, zlabelText)
    % Remove inf/nan rows for gridding
    good = isfinite(cpu) & isfinite(err) & isfinite(z) & err > 0 & cpu > 0;
    cpu = cpu(good); err = err(good); z = z(good);

    % Build grid in log space (better spread)
    x = log10(cpu); y = log10(err);
    xi = linspace(min(x), max(x), 50);
    yi = linspace(min(y), max(y), 50);
    [X,Y] = meshgrid(xi, yi);

    Z = griddata(x, y, z, X, Y, "natural");

    contourf(10.^X, 10.^Y, Z, 12);
    set(gca,"XScale","log","YScale","log");
    colorbar; ylabel(colorbar, zlabelText);
    grid on;
    xlabel("CPU time (s)");
    ylabel("Max error (rad/s)");
end

function lbl = solverLabel(s)
    switch char(s)
        case 'ode1'
            lbl = 'Euler (ode1)';
        case 'ode4'
            lbl = 'Runge–Kutta 4 (ode4)';
        case 'ode45'
            lbl = 'ode45 (adaptive RK)';
        case 'ode23tb'
            lbl = 'ode23tb (stiff)';
        otherwise
            lbl = char(s);
    end
end
