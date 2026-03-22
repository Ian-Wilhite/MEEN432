if ~exist("cycleName", "var")
    cycleName = "drive_cycle";
end

% Use the folder containing this script
scriptDir = fileparts(mfilename("fullpath"));
cd(scriptDir)

% Load and simulate the model
load_system("p3_car.slx");
simout = sim("p3_car");

% Retrieve vehicle speed signal (m/s)
if isprop(simout, "veh_speed") && ~isempty(simout.veh_speed)
    sim_vel = simout.veh_speed.Data;
elseif isprop(simout, "veh_speed_dyn") && ~isempty(simout.veh_speed_dyn)
    sim_vel = simout.veh_speed_dyn.Data;
else
    error("Vehicle speed signal not found in SimulationOutput.");
end

sim_time = simout.tout;
mph2mps = 1609.344/3600;

% Retrieve motor electrical power
if isprop(simout, "motor_power_elec_W") && ~isempty(simout.motor_power_elec_W)
    motorPowerW = simout.motor_power_elec_W.Data;
else
    error("motor_power_elec_W signal not found in SimulationOutput.");
end

% Energy calculation
energy_J = trapz(sim_time, max(motorPowerW,0));
energy_kWh = energy_J / 3.6e6;

fprintf("%s: Energy consumed = %.4f kWh\n", cycleName, energy_kWh);

cumEnergy_kWh = cumtrapz(sim_time, max(motorPowerW,0)) / 3.6e6;

figure;
plot(sim_time, motorPowerW)
xlabel("Time (s)")
ylabel("Motor Electrical Power (W)")
title("Motor Electrical Power vs Time: " + cycleName)

figure;
plot(sim_time, cumEnergy_kWh)
xlabel("Time (s)")
ylabel("Cumulative Energy (kWh)")
title("Cumulative Energy Consumed: " + cycleName)

sim_vel_mph = sim_vel * (1 / mph2mps);
drive_vel_mph = DriveData;

% Compare only over overlapping time window
t_end = min(sim_time(end), Time(end));
dc_mask = Time <= t_end;
sim_mask = sim_time <= t_end;

interp_sim = interp1(sim_time(sim_mask), sim_vel_mph(sim_mask), Time(dc_mask), "linear", "extrap");
speed_error = interp_sim - drive_vel_mph(dc_mask);
max_abs_err = max(abs(speed_error));

fprintf("%s: Max |speed error| = %.2f mph (target <= 3 mph)\n", cycleName, max_abs_err);

figure;
plot(sim_time, sim_vel_mph, 'b')
hold on
plot(Time, drive_vel_mph, '--r')
plot(Time, drive_vel_mph + 3, '--k')
plot(Time, drive_vel_mph - 3, '--k')
xlabel("Time (s)")
ylabel("Velocity (mph)")
legend("Sim Velocity", "Drive Cycle Velocity", "3 mph Error Band")
title("Simulated Vehicle Velocity vs Time: " + cycleName)

% Save the plot
assetsDir = fullfile(scriptDir, "assets");
if ~exist(assetsDir, "dir")
    mkdir(assetsDir);
end
exportgraphics(gcf, fullfile(assetsDir, cycleName + "_plot.png"));