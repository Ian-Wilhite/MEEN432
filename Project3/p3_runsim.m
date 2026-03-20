if ~exist("cycleName", "var")
    cycleName = "drive_cycle";
end

simout = sim("p3_car.slx");

% Retrieve vehicle speed signal (m/s) regardless of logged variable name
if isprop(simout, "veh_speed") && ~isempty(simout.veh_speed)
    sim_vel = simout.veh_speed.Data;
elseif isprop(simout, "veh_speed_dyn") && ~isempty(simout.veh_speed_dyn)
    sim_vel = simout.veh_speed_dyn.Data;
elseif simout.has("veh_speed_dyn")
    sim_vel = simout.get("veh_speed_dyn");
elseif simout.has("veh_speed")
    sim_vel = simout.get("veh_speed");
else
    error("Vehicle speed signal not found in SimulationOutput.");
end

sim_time = simout.tout;          % s
mph2mps = 1609.344/3600;         % m/s per mph

sim_vel_mph = sim_vel * (1 / mph2mps);
drive_vel_mph = DriveData; % provided drive cycle is already mph

% Compare only over the overlapping time window to avoid extrapolation blow-up
t_end = min(sim_time(end), Time(end));
dc_mask = Time <= t_end;
sim_mask = sim_time <= t_end;

interp_sim = interp1(sim_time(sim_mask), sim_vel_mph(sim_mask), Time(dc_mask), "linear", "extrap");
speed_error = interp_sim - drive_vel_mph(dc_mask);
max_abs_err = max(abs(speed_error));
fprintf("%s: Max |speed error| = %.2f mph (target <= 3 mph)\n", cycleName, max_abs_err);

figure;
plot(sim_time, sim_vel_mph, 'b') % Remember, drive cycles are mph
hold on
plot(Time, drive_vel_mph, '--r')
plot(Time, drive_vel_mph + 3, '--k')
plot(Time, drive_vel_mph - 3, '--k')
xlabel("Time (s)")
ylabel("Velocity (mph)")
legend("Sim Velocity", "Drive Cycle Velocity", "3 mph Error Band")
title("Simulated Vehicle Velocity vs Time: " + cycleName)

% Save the plot for reporting (kept inside Project3/assets)
scriptDir = fileparts(mfilename("fullpath"));
assetsDir = fullfile(scriptDir, "assets");
if ~exist(assetsDir, "dir")
    mkdir(assetsDir);
end
exportgraphics(gcf, fullfile(assetsDir, cycleName + "_plot.png"));

% For seeing how large the errors are
% error = zeros(length(Time),1);
% for j = 1:length(Time)
%     time_dc = Time(j);
%     vel_dc = DriveData(j);
%     for i = 1:length(sim_time)
%         time_s = sim_time(i);
%         vel_s = sim_vel(i);
% 
%         if time_s == time_dc
%             err = vel_dc - vel_s;
%             error(j) = err;
%         else
%         end
%     end
% end
