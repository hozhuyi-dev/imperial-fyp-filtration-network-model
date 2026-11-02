% The scripts run simulations when varying network size

clear;
clc;
close all;

% Add src to path
path(path, ['..', filesep, 'src']);

% Base particle sizes to consider
d = 0.15:0.005:0.35;

% Number of time steps
maxFinesTimeStep = 5;
particleNum = 'less';

% Filters to consider
networkType = 'actual';
cuCoarse = 1.2;
triangulationMethod = 'Del';
density = 'dense';
assignmentMethod = 'actual';
sizeList = [1,1,1; ...
    1,1,2; ...
    1,1,3; ...
    2,2,1; ...
    2,2,2; ...
    2,2,3];

% Options
continuous = 0;
loadMethod = 'read';

for i = 1:size(sizeList,1)
    % Initialise simulator
    simulator = Simulator('maxFinesTimeStep',maxFinesTimeStep);
    
    % Simulations
    simulator = simulator.runnerMonodisperse(cuCoarse, d, 'density', density, 'size', sizeList(i,:), ...
        'triangulationMethod', triangulationMethod, 'networkType', networkType, 'assignmentMethod', assignmentMethod, ...
        'continuous', continuous, 'loadMethod', loadMethod, 'particleNum', particleNum);
end