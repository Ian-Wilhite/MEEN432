clc;
clear;

if ~isfile("Project_2_Kinematic_Model.slx")
    error("Simulink model Project_2_Kinematic_Model.slx not found in Project2 directory.");
end

init;      % vehicle parameters & initial conditions
gentrack;  % track definition (exports 'path' to base workspace)

try
    animate;  % runs Simulink model and animation
catch err
    fprintf(2,"Error during simulation/animation: %s\n", err.message);
end
