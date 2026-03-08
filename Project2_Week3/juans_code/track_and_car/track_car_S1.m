clc;
clear;

% Calling Track and Car info (matches starter track geometry)
[scenario, egoVehicle] = track_car_info();

% figure set up
fig = figure('Name', 'Driving Scenario Animation', 'Position', [100, 100, 1200, 600]);

% 2D Plot of Scenario
ax2d = subplot(1, 2, 1);
plot(scenario, 'Parent', ax2d);
title('2D Top-Down View');
grid on;

% 3D View of Scenario
ax3d = subplot(1, 2, 2);
chasePlot(egoVehicle, 'Parent', ax3d);
title('3D Perspective View');
view(ax3d, 45, 30); % viewing angle
grid on;
axis(ax3d, 'equal');

while advance(scenario)
    % update scenario plots until time is done
    drawnow limitrate;
end
