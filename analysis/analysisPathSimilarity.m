% The script evaluates the similarity of erodible paths for consecutive
% runs

clear
clc
close all

% Add src to path
path(path, ['..', filesep, 'src'])

% Networks to consider
transportLogicType = 'flowbiased';
networkType = 'actual';
assignmentMethod = 'actual';
cuCoarse = 1.2;
density = 'dense';
size = [1,1,1];
triangulationMethod = 'Del';
dList = 0.25:0.005:0.28;
partNum = 10;

% Number of base
baseNum = length(dList);

% Initialise S
S = zeros(baseNum,partNum-1);

for i = 1:baseNum
    % Load simulation data
    simulator = Simulator.loadResults(FileLocation.simulatorMonodisperseParts(transportLogicType,networkType, ...
        assignmentMethod,cuCoarse,density,size,triangulationMethod,dList(i),0));
    
    % For each consecutive runs
    for j = 1:partNum-1
        % Determine the indices of edges passed in run i
        v1 = simulator.statistics.(sprintf('statistics%i',j)).edgePathEroded>0;
        
        % Indices of edges passed in run i+1
        v2 = simulator.statistics.(sprintf('statistics%i',j+1)).edgePathEroded>0;
        
        % Determine similarities
        S(i,j) = Statistics.similarityIndex(v1,v2);
    end
end

%% 

for i = 1:baseNum
    fprintf('%.3f & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f\n',S(i,:));
end

%% Plot

% Figure properties
colourList = Constants.colourList;
legendList = {'Type 3 - Clogged nodes','Type 2 - Clogged nodes', ...
    'Type 3 - Node volume','Type 2 - Node volume'};

% Plot
figure;
set(gcf,'WindowState','maximized');
plot(1:5,cActual,'LineWidth',5,'Color',colourList(1,:),'Marker','o','MarkerSize',30);
hold on;
plot(1:5,cRandom,'LineWidth',5,'Color',colourList(2,:),'Marker','x','MarkerSize',30);
plot(1:5,cVolumeActual,'LineWidth',5,'Color',colourList(1,:),'LineStyle',':','Marker','o','MarkerSize',30);
plot(1:5,cVolumeRandom,'LineWidth',5,'Color',colourList(2,:),'LineStyle',':','Marker','x','MarkerSize',30);
legend(legendList,'Location','northeast','FontName','Helvetica','FontSize',36);
xlim([1,5]);
xticks(1:5);
set(gca, 'FontName','Helvetica','FontSize',36);
xlabel('Number of edges apart','FontName','Helvetica','FontSize',36);
ylabel('Correlation coefficient','FontName','Helvetica','FontSize',36);