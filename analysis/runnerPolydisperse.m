% The script performs polydisperse simulations

clear
clc
close all

% Add src to path
path(path, ['..', filesep, 'src'])

% Base to consider
r = 5;
cuFine = 1.5;

% Filter to consider
transportLogicType = 'flowbiased';
networkType = 'actual';
cuCoarseList = 1.2;
triangulationMethod = 'Del';
densityList = {'dense'};
assignmentMethod = 'actual';
size = [1,1,1];

% Option(s)
run = 1;
repeatNum = 10;
loadMethod = 'read';

for i = 1:length(cuCoarseList)
    for j = 1:length(densityList)
        % Initialise simulator
        simulator = Simulator('transportLogicType',transportLogicType);
        
        % Perform simulations
        simulator = simulator.repeatPolydisperse(cuCoarseList(i), r, cuFine, 'density', densityList{j}, 'size', size, ...
            'triangulationMethod', triangulationMethod, 'networkType', networkType, 'assignmentMethod', assignmentMethod, ...
            'loadMethod', loadMethod, 'partNum', repeatNum, 'run', 1);
    end
end