% The script plots figures showing the validations when compared to Kenney
% et. al. (1985) and Roozbahani et. at. (2014)

clear;
clc;
close all;

% Add src to path
path(path, ['..', filesep, 'src']);

% Networks to consider
transportLogicType = 'flowbiased';
networkType = 'actual';
assignmentMethod = 'actual'; 
cuCoarseList = [1.2,3];
densityList = {'dense','loose'};
sizeList = [1,1,1; ...
    1,1,2; ...
    1,1,3; ...
    2,2,1; ...
    2,2,2; ...
    2,2,3];
triangulationMethod = 'Del';

% Number of networks
networkNum = 9;

% Initialise struct(s)
dNormal = struct;
pNormal = struct;
dLess = struct;
pLess = struct;

% Initialise legendList
legendList = cell(length(cuCoarseList)+networkNum,1);

% Load Simulator objects
legendIdx = 0;
k = 0;
for i = 1:length(cuCoarseList)
    % Current cuCoarse
    cuCoarse = cuCoarseList(i);
    density = densityList{i};
    sizeCurrent = sizeList(1,:);
    
    % Load results
    simulator = Simulator.loadResults(FileLocation.simulatorMonodisperse(transportLogicType, ...
        networkType, assignmentMethod, cuCoarse, density, sizeCurrent, triangulationMethod));
    
    % Size of base particles
    dNormal.(sprintf('d%i',i)) = simulator.monodisperseBaseDiameter();
    
    % Percentage of particles eroded
    pNormal.(sprintf('p%i',i)) = simulator.percentageErodedNum();
    
    % legendList
    legendIdx = legendIdx+1;
    legendList{legendIdx} = sprintf('F%i%s%i%i%i More',round(cuCoarse*10),density,sizeCurrent(1),sizeCurrent(2),sizeCurrent(3));
    
    if cuCoarse==3
        sizeNum = 3;
    else
        sizeNum = size(sizeList,1);
    end
    
    % For each size
    for j = 1:sizeNum
        % Current iteration
        k = k+1;
        sizeCurrent = sizeList(j,:);
        
        % Load results
        simulator = Simulator.loadResults(FileLocation.simulatorMonodisperseLess(transportLogicType, ...
            networkType, assignmentMethod, cuCoarse, density, sizeCurrent, triangulationMethod));
        
        % Size of base particles
        dLess.(sprintf('d%i',k)) = simulator.monodisperseBaseDiameter();
        
        % Percentage of particles eroded
        pLess.(sprintf('p%i',k)) = simulator.percentageErodedNum();
        
        % legendList
        legendIdx = legendIdx+1;
        legendList{legendIdx} = sprintf('F%i%s%i%i%i Less',round(cuCoarse*10),density,sizeCurrent(1),sizeCurrent(2),sizeCurrent(3));
    end
end

%%

% Figure properties
lineStyleList = {'-','--',':','none','none','none'};
markerList = {'none','none','none','x','o','+'};
colourList = Constants.colourList;

% Update legendList
% legendList = vertcat(legendList, ...
%     {'Roozbahani et. al. (2014) C_U = 1.0'; ...
%     'Kenney et. al. (1985) C_U = 1.2'; ...
%     'Kenney et. al. (1985) C_U = 3.0'});
legendList = {'Filter C_U = 1.2','Filter C_U = 3.0'};

%%

% Plot
figure;
set(gcf,'WindowState','maximized');
plot(dNormal.(sprintf('d%i',1)),pNormal.(sprintf('p%i',1)),'LineWidth',5,'LineStyle','-.','Color',colourList(1,:));
hold on;
% for i = 1:6
%     plot(dLess.(sprintf('d%i',i)),pLess.(sprintf('p%i',i)),'LineWidth',2, ...
%         'LineStyle',lineStyleList{(i==6)*6+mod(i,6)},'Color',colourList((i>6)+1,:), ...
%         'Marker',markerList{(i==6)*6+mod(i,6)},'MarkerSize',10);
% end
plot(dNormal.(sprintf('d%i',2)),pNormal.(sprintf('p%i',2)),'LineWidth',5,'LineStyle','-.','Color',colourList(2,:));
% for i = 7:9
%     plot(dLess.(sprintf('d%i',i)),pLess.(sprintf('p%i',i)),'LineWidth',2, ...
%         'LineStyle',lineStyleList{(i==6)*6+mod(i,6)},'Color',colourList((i>6)+1,:), ...
%         'Marker',markerList{(i==6)*6+mod(i,6)},'MarkerSize',10);
% end
% plot(0.27*ones(1,2),ylim,'LineWidth',3,'LineStyle',':','Color',colourList(1,:));
% plot(0.18*ones(1,2),ylim,'LineWidth',3,'LineStyle','--','Color',colourList(1,:));
% plot(0.25*ones(1,2),ylim,'LineWidth',3,'LineStyle','--','Color',colourList(2,:));
legend(legendList,'Location','northeast','FontName','Helvetica','FontSize',36);
xlim([0.15,0.45]);
set(gca, 'FontName','Helvetica','FontSize',36);
xlabel('Normalised diameter of monodisperse base particles [-]','FontName','Helvetica','FontSize',36);
ylabel('Percentage of base particles eroded','FontName','Helvetica','FontSize',36);