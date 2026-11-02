% The script determines the range of z-coordinates for inlets and outlets

clear;
clc;
close all;

% Add src to path
path(path, ['..', filesep, 'src']);

% List of filters
cuCoarseList = Constants.cuCoarse(1:4);
densityList = Constants.density;
size = [1,1,1];
triangulationMethod = 'Del';

% Initialise array(s)
minInlet = zeros(length(cuCoarseList),length(densityList));
maxInlet = minInlet;
minOutlet = minInlet;
maxOutlet = minInlet;

% Number of networks
n = length(cuCoarseList)*length(densityList);

for i = 1:length(cuCoarseList)
    for j = 1:length(densityList)
        % Load filter
        simulator = Simulator();
        simulator = simulator.createFilter(cuCoarseList(i),'density',densityList{j},'size',size,'triangulationMethod',triangulationMethod);
        
        % Store the values
        minInlet(i,j) = min(simulator.filter.nodeCoord(simulator.filter.nodeInlet==1,3));
        maxInlet(i,j) = max(simulator.filter.nodeCoord(simulator.filter.nodeInlet==1,3));
        minOutlet(i,j) = min(simulator.filter.nodeCoord(simulator.filter.nodeOutlet==1,3));
        maxOutlet(i,j) = max(simulator.filter.nodeCoord(simulator.filter.nodeOutlet==1,3));
    end
end

minInlet = minInlet';
minInlet = minInlet(:);
maxInlet = maxInlet';
maxInlet = maxInlet(:);
minOutlet = minOutlet';
minOutlet = minOutlet(:);
maxOutlet = maxOutlet';
maxOutlet = maxOutlet(:);

% Plot
figure;
set(gcf,'WindowState','maximized');
plot(1:n,minInlet,'LineStyle','none','Marker','x','LineWidth',5);
hold on;
plot(1:n,maxInlet,'LineStyle','none','Marker','x','LineWidth',5);
plot(1:n,minOutlet,'LineStyle','none','Marker','x','LineWidth',5);
plot(1:n,maxOutlet,'LineStyle','none','Marker','x','LineWidth',5);
xlim([0,n+1]);
xticks(0:n+1);