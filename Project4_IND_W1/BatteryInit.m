%% This file loads the Battery Data

datBat.SOC=[0, 0.01,   .1,  .2,  .3,  .4,  .5,  .6, .7,  .8,  .9,   1]; % soc
datBat.OCV=[0, 3.1,3.55,3.68,3.74,3.76,3.78,3.85,3.9,3.95,4.08,4.15];  % OCV per cell
datBat.Rint=0.1695; % internal resistance per cell
datBat.C = 150; % Amp hr total battery capacity
datBat.numSeries = 96;
datBat.numParallel = 74;

%{
figure
plot(datBat.SOC,datBat.OCV)
xlabel('Cell State of Charge ')
ylabel('Cell Open Circuit Volage (OCV) - Volts')
title('Lithium Ion Cell Characteristic - Project 3')
grid
%}
