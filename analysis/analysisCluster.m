% The script plots the graphs showing the importance of spatial correlation
% in filtration simulations

clear;
clc;
close all;

% Add src to path
path(path, ['..', filesep, 'src']);

% Load actual
simulator = Simulator.loadResults(FileLocation.simulatorMonodisperse('flowbiased','actual','actual',1.2,'dense',[1,1,1],'Del'));
cActual = simulator.statistics.statistics18.spatialCorrelationFullNodes(0.235,simulator.filter,0.5,'distanceUpTo',5);
cVolumeActual = simulator.filter.spatialCorrelationNodeVolume('distanceUpTo',5);

% Load random
simulator = Simulator.loadResults(FileLocation.simulatorMonodisperse('flowbiased','actual','random',1.2,'dense',[1,1,1],'Del'));
cRandom = simulator.statistics.statistics18.spatialCorrelationFullNodes(0.235,simulator.filter,0.5,'distanceUpTo',5);
cVolumeRandom = simulator.filter.spatialCorrelationNodeVolume('distanceUpTo',5);

% Figure properties
colourList = Constants.colourList;
legendList = {'Type 3 - Clogged nodes','Type 2 - Clogged nodes', ...
    'Type 3 - Node volume','Type 2 - Node volume'};

% Plot
figure;
set(gcf,'WindowState','maximized');
plot(1:5,cActual,'LineWidth',5,'Color',colourList(1,:),'Marker','o','MarkerSize',30);
hold on;
plot(1:5,cRandom,'LineWidth',5,'Color',colourList(3,:),'Marker','x','MarkerSize',30);
plot(1:5,cVolumeActual,'LineWidth',5,'Color',colourList(1,:),'LineStyle',':','Marker','o','MarkerSize',30);
plot(1:5,cVolumeRandom,'LineWidth',5,'Color',colourList(3,:),'LineStyle',':','Marker','x','MarkerSize',30);
legend(legendList,'Location','northeast','FontName','Helvetica','FontSize',36);
xlim([1,5]);
xticks(1:5);
set(gca, 'FontName','Helvetica','FontSize',36);
xlabel('Number of edges apart','FontName','Helvetica','FontSize',36);
ylabel('Correlation coefficient','FontName','Helvetica','FontSize',36);