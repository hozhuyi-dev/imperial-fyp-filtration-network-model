% Script runs simulations required for evaluation of network type on
% percentage of particles eroded

clear;
clc;
close all;

% Add src to path
path(path, ['..', filesep, 'src']);

% Base particle sizes to consider
d = 0.15:0.005:0.35;

% Filters to consider
transportLogicType = 'flowbiasedpos';
networkTypeList = {'actual'};%{'actual','squarenp','slantedsquarenp','cubicnp','slantedcubicnp','bccnp','fccnp'};
cuCoarse = 1.2;
triangulationMethod = 'Del';
density = 'dense';
assignmentMethod = 'actual';
size = [1,1,1];

% Options
continuous = 0;
loadMethod = 'read';

for i = 1:length(networkTypeList)
    % Initialise simulator
    simulator = Simulator('transportLogicType',transportLogicType);
    
    % Simulations
    simulator = simulator.runnerMonodisperse(cuCoarse, d, 'density', density, 'size', size, ...
        'triangulationMethod', triangulationMethod, 'networkType', networkTypeList{i}, 'assignmentMethod', assignmentMethod, ...
        'continuous', continuous, 'loadMethod', loadMethod);
end