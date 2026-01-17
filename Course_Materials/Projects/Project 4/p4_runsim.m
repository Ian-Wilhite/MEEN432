simout = sim("p4_car.slx");
sim_vel = simout.veh_speed.Data;
sim_time = simout.tout;

figure;
%plot(sim_time, sim_vel*(1/mph2mps), 'b') % Remember, drive cycles are mph
%hold on
plot(Time, DriveData, '--r') 
plot(Time, (DriveData)+3, '--k') 
plot(Time, (DriveData)-3, '--k') 
xlabel("Time (s)")
ylabel("Velocity (mph)") 
legend("Sim Velocity", "Drive Cycle Velocity", "3 mph Error Band")
title("Simulated Vehicle Velocity vs Time for a certain Drive Cycle")

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