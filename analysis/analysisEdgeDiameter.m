% The script plots graphs on how cuCoarse and density affects percentage
% eroded

clear;
clc;
close all;

% Add src to path
path(path, ['..', filesep, 'src']);

%% Load

% Networks to consider
transportLogicType = 'flowbiased';
networkType = 'actual';
assignmentMethod = 'actual';
cuCoarseList = Constants.cuCoarse;
densityList = Constants.density;
size = [1,1,1];
triangulationMethod = 'Del';

% Initialise array(s)
p = struct;
d = struct;
r = struct;
nEroded = struct;
nRetained = struct;
pAffected = struct;

% Load data
k = 0;
for i = 1:length(cuCoarseList)
    for j = 1:length(densityList)
        k = k+1;
        
        % Load Simulator object
        simulator = Simulator.loadResults(FileLocation.simulatorMonodisperse(transportLogicType,networkType, ...
            assignmentMethod,cuCoarseList(i),densityList{j},size,triangulationMethod));
        
        % Determine the diameter of particles
        d.(sprintf('d%i',k)) = simulator.monodisperseBaseDiameter;
        
        % Determine D15F
        % D15F = simulator.filter.percentileParticleDiameter(15);
        % d.(sprintf('d%i',k)) = d.(sprintf('d%i',k))/D15F;
        
        % Determine the percentage of particles eroded
        p.(sprintf('p%i',k)) = simulator.percentageErodedNum;
        
        % Determine the percentage of nodes affected
        pAffected.(sprintf('p%i',k)) = simulator.percentageNodeAffected;
        
        % Determine the ratio of average distance travelled to shortest path
        r.(sprintf('r%i',k)) = simulator.directness();
        
        % Determine the average number of constrictions met
        nEroded.(sprintf('n%i',k)) = simulator.constrictionNum();
        nRetained.(sprintf('n%i',k)) = simulator.constrictionNum('mechanism','retained');
        
        fprintf('Loaded Cu = %.1f %6s %6.2f%%\n', cuCoarseList(i), densityList{j}, k/length(cuCoarseList)/length(densityList)*100);
    end
end
clear s;

%% Figures

% Figure properties
lineStyleList = {'-','--',':'};
colourList = Constants.colourList;
legendList = cell(length(cuCoarseList)*length(densityList),1);
k = 0;
for i = 1:length(cuCoarseList)
    for j = 1:length(densityList)
        k = k+1;
        legendList{k} = sprintf('F%i%s111', round(cuCoarseList(i)*10), densityList{j});
    end
end

% Plot percentage eroded by num against d
figure;
set(gcf,'WindowState','maximized');
k = 0;
for i = 1:length(cuCoarseList)
    for j = 1:length(densityList)
        k = k+1;
        plot(d.(sprintf('d%i',k)),p.(sprintf('p%i',k)),'LineWidth',3,'LineStyle',lineStyleList{j},'Color',colourList(i,:));
        hold on;
    end
end
plot(0.2*ones(1,2),ylim,'LineWidth',3,'LineStyle','--','Color','k');
text(0.195,2,'Kenney et. al. (1985)','FontName','Helvetica','FontSize',30,'HorizontalAlignment','left','Rotation',90)
% text(0.2,15,'D_c^* = 0.2 D_1_5_F','FontName','Helvetica','FontSize',30,'HorizontalAlignment','right')
xlim([0.15,0.45])
legend(legendList,'Location','northeast','FontName','Helvetica','FontSize',34);
set(gca, 'FontName','Helvetica','FontSize',36);
xlabel('Normalised diameter of monodisperse base particles by D_1_5_F [-]','FontName','Helvetica','FontSize',36);
ylabel('Percentage of base particles eroded','FontName','Helvetica','FontSize',36);

%%
% Plot percentage eroded by num against d
figure;
set(gcf,'WindowState','maximized');
k = 0;
for i = 1:length(cuCoarseList)
    for j = 1:length(densityList)
        k = k+1;
        plot(d.(sprintf('d%i',k)),pAffected.(sprintf('p%i',k)),'LineWidth',3,'LineStyle',lineStyleList{j},'Color',colourList(i,:));
        hold on;
    end
end
xlim([0.15,0.45])
ylim([0,100])
legend(legendList,'Location','northeast','FontName','Helvetica','FontSize',18);
set(gca, 'FontName','Helvetica','FontSize',18);
xlabel('Normalised diameter of monodisperse base particles [-]','FontName','Helvetica','FontSize',18);
ylabel('Percentage of base particles eroded','FontName','Helvetica','FontSize',18);

%%

% Plot directness ratio against d
figure;
set(gcf,'WindowState','maximized');
k = 0;
for i = 1:length(cuCoarseList)
    for j = 1:length(densityList)
        k = k+1;
        plot(d.(sprintf('d%i',k)),r.(sprintf('r%i',k)),'LineWidth',2,'LineStyle',lineStyleList{j},'Color',colourList(i,:));
        hold on;
    end
end
legend(legendList,'Location','northeast','FontName','Helvetica','FontSize',16);
set(gca, 'FontName','Helvetica','FontSize',16);
xlabel('Normalised diameter of monodisperse base particles [-]','FontName','Helvetica','FontSize',16);
ylabel('Ratio of average distance to shortestpath','FontName','Helvetica','FontSize',16);

% Plot number of constrictions against d
figure;
set(gcf,'WindowState','maximized');
k = 0;
for i = 1:length(cuCoarseList)
    for j = 1:length(densityList)
        k = k+1;
        plot(d.(sprintf('d%i',k)),nEroded.(sprintf('n%i',k)),'LineWidth',2,'LineStyle',lineStyleList{j},'Color',colourList(i,:));
        hold on;
    end
end
legend(legendList,'Location','northeast','FontName','Helvetica','FontSize',16);
set(gca, 'FontName','Helvetica','FontSize',16);
xlabel('Normalised diameter of monodisperse base particles [-]','FontName','Helvetica','FontSize',16);
ylabel('Average number of constrictions','FontName','Helvetica','FontSize',16);

% Plot number of constrictions against d
figure;
set(gcf,'WindowState','maximized');
k = 0;
for i = 1:length(cuCoarseList)
    for j = 1:length(densityList)
        k = k+1;
        plot(d.(sprintf('d%i',k)),nRetained.(sprintf('n%i',k)),'LineWidth',2,'LineStyle',lineStyleList{j},'Color',colourList(i,:));
        hold on;
    end
end
legend(legendList,'Location','northeast','FontName','Helvetica','FontSize',16);
set(gca, 'FontName','Helvetica','FontSize',16);
xlabel('Normalised diameter of monodisperse base particles [-]','FontName','Helvetica','FontSize',16);
ylabel('Average number of constrictions','FontName','Helvetica','FontSize',16);

%%

% Percentages to consider
percentage = [95,80,60,40,20,5];

% Initialise array(s)
ratioDensity = zeros(length(cuCoarseList),length(percentage));

for i = 1:length(cuCoarseList)
    % Loose and dense indices
    idxLoose = (i-1)*3+1;
    idxDense = (i*3);
    
    % Current p and d
    pDense = p.(sprintf('p%i',idxDense));
    pLoose = p.(sprintf('p%i',idxLoose));
    dDense = d.(sprintf('d%i',idxDense));
    dLoose = d.(sprintf('d%i',idxLoose));
    
    % Unique points
    dDense = dDense(pDense~=0);
    dLoose = dLoose(pLoose~=0);
    pDense = pDense(pDense~=0);
    pLoose = pLoose(pLoose~=0);
    
    for j = 1:length(percentage)
        ratioDensity(i,j) = interp1(pDense,dDense,percentage(j))/interp1(pLoose,dLoose,percentage(j));
    end
end