% The script plots graphs related to how network size affects percentage
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
cuCoarse = 1.2;
densityList = {'dense'};
sizeList = [1,1,1; ...
    1,1,2; ...
    1,1,3; ...
    2,2,1; ...
    2,2,2; ...
    2,2,3];
triangulationMethod = 'Del';
networkNum = size(sizeList,1);

% Initialise array(s)
p = struct;
d = struct;
dRetained = struct;
r = struct;
nEroded = struct;
nRetained = struct;

% Load data
j = 0;
for h = 1:length(densityList)
    for i = 1:networkNum
        j = j+1;
        
        % Load Simulator object
        simulator = Simulator.loadResults(FileLocation.simulatorMonodisperseLess(transportLogicType,networkType, ...
            assignmentMethod,cuCoarse,densityList{h},sizeList(i,:),triangulationMethod));
        
        % Determine the diameter of particles
        d.(sprintf('d%i',j)) = simulator.monodisperseBaseDiameter();
        
        % Determine the percentage of particles eroded
        p.(sprintf('p%i',j)) = simulator.percentageErodedNum();
        
        % Determine the ratio of average distance travelled to shortest path
        r.(sprintf('r%i',j)) = simulator.directness();
        
        % Determine the average number of constrictions met
        nEroded.(sprintf('n%i',j)) = simulator.constrictionNum();
        nRetained.(sprintf('n%i',j)) = simulator.constrictionNum('mechanism','retained');
        
        % Distance travelled
        dRetained.(sprintf('d%i',j)) = simulator.distanceTravelled('mechanism','retained');
        if i<=1
            [~,shortest] = simulator.filter.shortestPath;
        end
        rRetained.(sprintf('r%i',j)) = dRetained.(sprintf('d%i',j))/shortest;
        
        if i==3
            % Divide filters into 10 layers per cube
            layerNum = simulator.filter.filterParameters.size(3)*10;
            
            % Number of particles retained in each layer
            nParticleRetained113 = simulator.layerRetainedNum(layerNum);
        end
        
        fprintf('Loaded %6.2f%%\n', i/networkNum*100);
    end
end
clear s;

%% Figures

% Figure properties
lineStyleList = {'-','-','-','none','none','none'};
markerList = {'none','none','none','x','o','+'};
colourList = Constants.colourList;
% legendList = cell(networkNum+length(densityList),1);
% j = 0;
% for i = 1:length(densityList)
%     j = j+1;
%     legendList{j} = sprintf('Density: %s',densityList{i});
% end
% for i = 1:networkNum
%     j = j+1;
%     legendList{j} = sprintf('Size [%i,%i,%i]',sizeList(i,1),sizeList(i,2),sizeList(i,3));
% end
legendList = cell(networkNum,1);
for i = 1:networkNum
    legendList{i} = sprintf('F12dense%i%i%i',sizeList(i,1),sizeList(i,2),sizeList(i,3));
end

%%

% Plot percentage eroded by num against d
figure
set(gcf,'WindowState','maximized')
for i = 1:length(densityList)
    plot(0,0,'LineWidth',3,'LineStyle',lineStyleList{1},'Marker','none','MarkerSize',20,'Color',colourList(i,:))
    hold on
end
for i = 1:networkNum
    plot(0,0,'LineWidth',3,'LineStyle',lineStyleList{i},'Marker',markerList{i},'MarkerSize',20,'Color','k')
end
j = 0;
for h = 1:length(densityList)
    for i = 1:networkNum
        j = j+1;
        plot(d.(sprintf('d%i',j)),p.(sprintf('p%i',j)),'LineWidth',5,'LineStyle',lineStyleList{i},'Marker',markerList{i},'MarkerSize',20,'Color',colourList(h,:));
        hold on;
    end
end
% plot((d.d1),((p.p1/100).^2)*100,'LineWidth',3,'LineStyle',lineStyleList{(i>2)+1});
% legend(legendList,'Location','northeast','FontName','Helvetica','FontSize',36);
xlim([min(d.d1),max(d.d1)]);
set(gca, 'FontName','Helvetica','FontSize',36);
xlabel('Normalised diameter of monodisperse base particles [-]','FontName','Helvetica','FontSize',36);
ylabel('Percentage of base particles eroded','FontName','Helvetica','FontSize',36);

%%
% Plot directness ratio against d
figure;
set(gcf,'WindowState','maximized');
for i = 1:networkNum
    plot(d.(sprintf('d%i',i)),rRetained.(sprintf('r%i',i)),'LineWidth',5,'LineStyle',lineStyleList{i},'Marker',markerList{i},'MarkerSize',20);
    hold on;
end
plot(0.1975*ones(1,2),ylim,'LineWidth',5,'Color','k','LineStyle','--');
plot(0.2175*ones(1,2),ylim,'LineWidth',5,'Color','k','LineStyle','--');
legend(legendList,'Location','northeast','FontName','Helvetica','FontSize',36);
xlim([0.15,0.35]);
set(gca, 'FontName','Helvetica','FontSize',36);
xlabel('Normalised diameter of monodisperse base particles [-]','FontName','Helvetica','FontSize',36);
ylabel('Normalised average distance','FontName','Helvetica','FontSize',36);

%%
% Plot number of constrictions against d
figure;
set(gcf,'WindowState','maximized');
for i = 1:networkNum
    plot(d.(sprintf('d%i',i)),nEroded.(sprintf('n%i',i)),'LineWidth',3,'LineStyle',lineStyleList{(i>2)+1});
    hold on;
end
legend(legendList,'Location','northwest','FontName','Helvetica','FontSize',18);
xlim([min(d.d1),max(d.d1)]);
set(gca, 'FontName','Helvetica','FontSize',18);
xlabel('Normalised diameter of monodisperse base particles [-]','FontName','Helvetica','FontSize',18);
ylabel('Average number of constrictions','FontName','Helvetica','FontSize',18);

% Plot number of constrictions against d
figure;
set(gcf,'WindowState','maximized');
for i = 1:networkNum
    plot(d.(sprintf('d%i',i)),nRetained.(sprintf('n%i',i)),'LineWidth',3,'LineStyle',lineStyleList{i},'Marker',markerList{i},'MarkerSize',18);
    hold on;
end
legend(legendList,'Location','northeast','FontName','Helvetica','FontSize',18);
xlim([min(d.d1),max(d.d1)]);
set(gca, 'FontName','Helvetica','FontSize',18);
xlabel('Normalised diameter of monodisperse base particles [-]','FontName','Helvetica','FontSize',18);
ylabel('Average number of constrictions','FontName','Helvetica','FontSize',18);

%%
% Plot number of particles retained in each layer
figure;
set(gcf,'WindowState','maximized');
for i = [8,12,18]
    % Plot
    plot(1:layerNum,nParticleRetained113(:,i),'LineWidth',5);
    hold on;
end
plot(10.5*ones(1,2),ylim,'LineWidth',5,'Color','k','LineStyle','--');
plot(20.5*ones(1,2),ylim,'LineWidth',5,'Color','k','LineStyle','--');
text(5.5,2000,'Within [1,1,1]','FontName','Helvetica','FontSize',36,'HorizontalAlignment','center')
text(15.5,2000,'Within [1,1,2]','FontName','Helvetica','FontSize',36,'HorizontalAlignment','center')
text(25.5,2000,'Within [1,1,3]','FontName','Helvetica','FontSize',36,'HorizontalAlignment','center')
legend({'d^n = 0.185','d^n = 0.205','d^n = 0.235'},'Location','northeast','FontName','Helvetica','FontSize',36);
xlim([0.5,30.5]);
xticks(0:2:30);
set(gca,'FontName','Helvetica','FontSize',36);
xlabel('Layer','FontName','Helvetica','FontSize',36);
ylabel('Number of retained particles','FontName','Helvetica','FontSize',36);

% Save figure
% print(fig,[FileLocation.dataFolder, 'figures', filesep, 'Monodisperse12dense111'],'-dpng','-r300')