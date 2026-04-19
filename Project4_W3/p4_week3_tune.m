%% p4_week3_tune.m
% Sweep target speeds and report the fastest VALID Week 3 run

clc; clear; clf;
p4_init;

candidate_mph = carData.week3_speed_candidates_mph;
summary = nan(numel(candidate_mph), 6);

fprintf('Sweeping candidate target speeds for Week 3...\n\n');

for k = 1:numel(candidate_mph)
    mph = candidate_mph(k);
    mps = mph * 0.44704;

    % Rebuild workspace each run
    p4_init;
    carData.vxd_mph = mph;
    carData.vxd     = mps;
    assignin('base','carData',carData);

    fprintf('Trying %.1f mph...\n', mph);
    simout = sim('p4_model.slx', 'StopTime', num2str(sim_stop_time));

    X_ts = simout.X; Y_ts = simout.Y;
    if isa(X_ts,'timeseries')
        X = X_ts.Data; Y = Y_ts.Data;
    else
        X = X_ts; Y = Y_ts;
    end

    vx_ts  = simout.vx_out;
    soc_ts = simout.soc_out;
    if isa(vx_ts,'timeseries')
        vx  = vx_ts.Data; soc = soc_ts.Data;
    else
        vx  = vx_ts; soc = soc_ts;
    end

    results = evaluate_week3_run_local(X(:), Y(:), vx(:), soc(:), path);

    summary(k,:) = [mph, results.laps, min(soc)*100, max(soc)*100, ...
                    results.max_cte, results.week3_valid];

    fprintf('  laps = %.2f | min SOC = %.2f%% | max CTE = %.2f m | valid = %d\n\n', ...
        results.laps, min(soc)*100, results.max_cte, results.week3_valid);
end

T = array2table(summary, 'VariableNames', ...
    {'TargetSpeed_mph','Laps','MinSOC_percent','MaxSOC_percent','MaxCTE_m','Valid'});
disp(T);

valid_idx = find(T.Valid > 0.5);
if isempty(valid_idx)
    fprintf('No valid candidate found. Reduce speed or reduce lookahead distance.\n');
else
    [~, best_rel] = max(T.Laps(valid_idx));
    best_idx = valid_idx(best_rel);
    fprintf('\nBest valid target speed: %.1f mph\n', T.TargetSpeed_mph(best_idx));
    fprintf('Laps completed: %.2f\n', T.Laps(best_idx));
end

function results = evaluate_week3_run_local(X, Y, vx, soc, path)
    cx = path.xpath(:);
    cy = path.ypath(:);
    nv = length(X);

    ds = sqrt(diff(cx).^2 + diff(cy).^2);
    cumS = [0; cumsum(ds)];
    L = cumS(end);

    arc = zeros(nv,1);
    cte = zeros(nv,1);
    for i = 1:nv
        d2 = (cx - X(i)).^2 + (cy - Y(i)).^2;
        [d2min, idx] = min(d2);
        arc(i) = cumS(idx);
        cte(i) = sqrt(d2min);
    end

    delta_arc = diff(arc);
    delta_arc(delta_arc < -L/2) = delta_arc(delta_arc < -L/2) + L;
    delta_arc(delta_arc >  L/2) = delta_arc(delta_arc >  L/2) - L;
    laps = sum(max(delta_arc,0)) / L;

    results.laps        = laps;
    results.max_cte     = max(abs(cte));
    results.on_track    = all(abs(cte) <= path.width/2);
    results.soc_ok      = all(soc >= 0.10) && all(soc <= 0.95);
    results.week3_valid = results.on_track && results.soc_ok;
end