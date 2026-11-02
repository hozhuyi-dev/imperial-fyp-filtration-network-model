% The script performs simulations with monodisperse base particles of size
% D3F and D5F; it is required that 3% and 5% of edges being clogged

clear;
clc;
close all;

% Add src/ to path
path(path, ['..', filesep, 'src']);

% Initialisation
simulator = Simulator();
simulator = simulator.createFilter(1.2);

% Determine D3F and D5F
D3F = simulator.filter.percentileEdgeDiameter(3);
D5F = simulator.filter.percentileEdgeDiameter(5);

% Create monodisperse base
simulator = simulator.simulateMonodisperse([D3F*ones(1,5),D5F*ones(1,5)]);

% Save
simulator.saveResults(FileLocation.simulatorValidationPercentageClogged);