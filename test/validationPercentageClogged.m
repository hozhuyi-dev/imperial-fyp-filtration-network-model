% The script performs simulations with monodisperse base particles of size
% D3F and D5F; it is required that 3% and 5% of edges being clogged

clear;
clc;
close all;

% Add src/ to path
path(path, ['..', filesep, 'src']);

% Load simulator
simulator = Simulator.loadResults(FileLocation.simulatorValidationPercentageClogged);

% Determine the percentage of edges clogged
p = simulator.percentageEdgeClogged(1);