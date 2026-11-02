% The script plots figures showing how the effectiveness of filters affect
% certain parameters (Monodisperse simulations

clear;
clc;
close all;

% Add src to path
path(path, ['..', filesep, 'src']);

% Networks to consider
transportLogicType = 'flowbiased';
networkType = 'actual';
assignmentMethod = 'actual';
cuCoarse = 1.2;
density = 'dense';
size = [1,1,1];
triangulationMethod = 'Del';
dList = 0.25:0.005:0.3;
% dList = 0.27;
partNum = 10;

% Number of base
baseNum = length(dList);

% Initialise struct
d = struct;
p = struct;

for i = 1:baseNum
    % Load simulation data
    simulator = Simulator.loadResults(FileLocation.simulatorMonodisperseParts(transportLogicType,networkType, ...
        assignmentMethod,cuCoarse,density,size,triangulationMethod,dList(i),0));
    
    % Determine the diameter of particles
    d.(sprintf('d%i',i)) = simulator.monodisperseBaseDiameter();
    
    % Determine percentage eroded
    p.(sprintf('p%i',i)) = simulator.percentageErodedNum('continuous',1);
    p.(sprintf('p%i',i))(isnan(p.(sprintf('p%i',i)))) = 0;
end

%% Plot

% Figure properties
markerList = {'x','o','^','s','d'};
colourList = Constants.colourList;
legendList = cell(partNum,1);
for i = 1:partNum
    legendList{i} = sprintf('Run %i',i);
end

% Plot
figure;
set(gcf,'WindowState','maximized');
for i = 1:baseNum
    for j = 1:partNum
        plot(d.(sprintf('d%i',i))(j),p.(sprintf('p%i',i))(j),'LineWidth',3,'LineStyle','none','Marker',markerList{ceil(j/2)},'MarkerSize',18,'Color',colourList(j,:));
        hold on;
    end
end
plot(0.2825*ones(1,2),ylim,'LineWidth',3,'LineStyle','--','Color','k');
plot(0.2925*ones(1,2),ylim,'LineWidth',3,'LineStyle','--','Color','k');
text(0.282,50,'Partial or no seal','FontName','Helvetica','FontSize',18,'HorizontalAlignment','right');
text(0.282,47,'with continuous erosion','FontName','Helvetica','FontSize',18,'HorizontalAlignment','right');
text(0.2875,48.5,'Seal with some erosion','FontName','Helvetica','FontSize',18,'HorizontalAlignment','center');
text(0.293,50,'Seal with','FontName','Helvetica','FontSize',18,'HorizontalAlignment','left');
text(0.293,47,'no erosion','FontName','Helvetica','FontSize',18,'HorizontalAlignment','left');
legend(legendList,'Location','northeast','FontName','Helvetica','FontSize',18);
xlim([min(dList),max(dList)]);
set(gca, 'FontName','Helvetica','FontSize',18);
xlabel('Normalised diameter of monodisperse base particles [-]','FontName','Helvetica','FontSize',18);
ylabel('Percentage of base particles eroded','FontName','Helvetica','FontSize',18);