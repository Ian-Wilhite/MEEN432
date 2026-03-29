if ~exist("cycleName", "var")
    cycleName = "drive_cycle";
end

scriptDir = fileparts(mfilename("fullpath"));
cd(scriptDir)

load_system("p3_car.slx");
simout = sim("p3_car");

mph2mps = 1609.344/3600;

% ---------------- Retrieve vehicle speed ----------------
if isprop(simout, "veh_speed") && ~isempty(simout.veh_speed)
    sim_vel = simout.veh_speed.Data;
elseif isprop(simout, "veh_speed_dyn") && ~isempty(simout.veh_speed_dyn)
    sim_vel = simout.veh_speed_dyn.Data;
else
    error("Vehicle speed signal not found in SimulationOutput.");
end

sim_time = simout.tout;

% ---------------- Retrieve battery power ----------------
if isprop(simout, "battery_power_W") && ~isempty(simout.battery_power_W)
    battPowerW = simout.battery_power_W.Data;
else
    error("battery_power_W signal not found. Log this signal from the model.");
end

% ---------------- Retrieve SOC ----------------
if isprop(simout, "soc_out") && ~isempty(simout.soc_out)
    soc = simout.soc_out.Data;
elseif isprop(simout, "SOC") && ~isempty(simout.SOC)
    soc = simout.SOC.Data;
else
    error("SOC signal not found. Log SOC from the model.");
end

% ---------------- Energy calculations ----------------
tractionPowerW = max(battPowerW, 0);      % battery discharging
regenPowerW    = max(-battPowerW, 0);     % recovered charging power

energyUsed_J      = trapz(sim_time, tractionPowerW);
energyRecovered_J = trapz(sim_time, regenPowerW);
netEnergy_J       = trapz(sim_time, battPowerW);

energyUsed_kWh      = energyUsed_J / 3.6e6;
energyRecovered_kWh = energyRecovered_J / 3.6e6;
netEnergy_kWh       = netEnergy_J / 3.6e6;

fprintf("\n%s results:\n", cycleName);
fprintf("  Battery energy used      = %.4f kWh\n", energyUsed_kWh);
fprintf("  Regen energy recovered   = %.4f kWh\n", energyRecovered_kWh);
fprintf("  Net battery energy       = %.4f kWh\n", netEnergy_kWh);
fprintf("  Initial SOC              = %.2f %%\n", 100*soc(1));
fprintf("  Final SOC                = %.2f %%\n", 100*soc(end));
fprintf("  SOC drop                 = %.2f %%\n", 100*(soc(1)-soc(end)));

cumNetEnergy_kWh = cumtrapz(sim_time, battPowerW) / 3.6e6;

% ---------------- Speed tracking check ----------------
sim_vel_mph = sim_vel / mph2mps;
drive_vel_mph = DriveData;

t_end = min(sim_time(end), Time(end));
dc_mask = Time <= t_end;
sim_mask = sim_time <= t_end;

interp_sim = interp1(sim_time(sim_mask), sim_vel_mph(sim_mask), Time(dc_mask), "linear", "extrap");
speed_error = interp_sim - drive_vel_mph(dc_mask);
max_abs_err = max(abs(speed_error));

fprintf("  Max |speed error|        = %.2f mph (target <= 3 mph)\n", max_abs_err);

% ---------------- Plots ----------------
fig1 = figure;
plot(sim_time, battPowerW)
xlabel("Time (s)")
ylabel("Battery Power (W)")
title("Battery Power vs Time: " + cycleName)
grid on

fig2 = figure;
plot(sim_time, cumNetEnergy_kWh)
xlabel("Time (s)")
ylabel("Cumulative Net Battery Energy (kWh)")
title("Net Battery Energy vs Time: " + cycleName)
grid on

fig3 = figure;
plot(sim_time, soc)
xlabel("Time (s)")
ylabel("SOC")
title("Battery State of Charge: " + cycleName)
grid on

fig4 = figure;
plot(sim_time, sim_vel_mph, 'b')
hold on
plot(Time, drive_vel_mph, '--r')
plot(Time, drive_vel_mph + 3, '--k')
plot(Time, drive_vel_mph - 3, '--k')
xlabel("Time (s)")
ylabel("Velocity (mph)")
legend("Sim Velocity", "Drive Cycle Velocity", "±3 mph Band", "Location", "best")
title("Simulated Vehicle Velocity vs Time: " + cycleName)
grid on

% ---------------- Save plots ----------------
assetsDir = fullfile(scriptDir, "assets");
if ~exist(assetsDir, "dir")
    mkdir(assetsDir);
end

% Determine output directory and filename suffix based on drive cycle
if cycleName == "fueleconomy"
    plotDir = fullfile(scriptDir, "Week3Plots_Fuel_Economy");
    suf = "Fuel_Economy";
elseif cycleName == "urban_driving"
    plotDir = fullfile(scriptDir, "Week3Plots_Urban_Driving");
    suf = "Urbandriving";
else
    plotDir = fullfile(scriptDir, "Week3Plots_" + cycleName);
    suf = cycleName;
end
if ~exist(plotDir, "dir")
    mkdir(plotDir);
end

exportgraphics(fig1, fullfile(plotDir, "Battery_Power_vs_Time_" + suf + ".pdf"));
exportgraphics(fig2, fullfile(plotDir, "Net_Battery_Energy_vs_Time_" + suf + ".pdf"));
exportgraphics(fig3, fullfile(plotDir, "Battery_State_of_Charge_" + suf + ".pdf"));
exportgraphics(fig4, fullfile(plotDir, "Simulated_Vehicle_Velocity_vs_Time_" + suf + ".pdf"));

fprintf("  Plots saved to %s\n", plotDir);