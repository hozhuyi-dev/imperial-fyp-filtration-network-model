% The script analyses the effect of network size on the percolation
% threshold for regular cubic networks

clear;
clc;
close all;

% Add src to path
path(path, ['..', filesep, 'src']);

% Colourlist
colourList = Constants.colourList;

%%

% Network sizes
networkSizeList = [10,10,10; ...
    20,20,20; ...
    30,30,30; ...
    40,40,40; ...
    50,50,50; ...
    60,60,60; ...
    10,10,30; ...
    20,20,30; ...
    30,30,30; ...
    40,40,30; ...
    50,50,30; ...
    60,60,30; ...
    30,30,10; ...
    30,30,20; ...
    30,30,30; ...
    30,30,40; ...
    30,30,50; ...
    30,30,60; ...
    30,30,30; ...
    30,30,60; ...
    30,30,90; ...
    30,30,120];

% Initialisation
periodicNum = size(networkSizeList,1);
n = 150;
p = NaN(periodicNum,n);

% Load file
for i = 1:periodicNum
    data = load(FileLocation.percolationThreshold('cubic',[],[],[],networkSizeList(i,:),1));
    p(i,:) = data.p;
end

% Create label list
labelList = cell(periodicNum,1);
for i = 1:periodicNum
    labelList{i} = sprintf('(%i,%i,%i)',networkSizeList(i,1),networkSizeList(i,2),networkSizeList(i,3));
end

% Plot
figure;
set(gcf,'WindowState','maximized');
boxplot(p','Labels',labelList,'Whisker',1000);
set(findobj(gca,'type','line'),'linew',3);
hold on;
plot(xlim,[0.2488,0.2488],'LineWidth',3,'LineStyle','--','Color','k');
ylim([0.14,0.32]);
plot(6.5*ones(2,1),ylim,'LineWidth',1,'Color','k');
plot(12.5*ones(2,1),ylim,'LineWidth',1,'Color','k');
plot(18.5*ones(2,1),ylim,'LineWidth',1,'Color','k');
set(gca, 'FontName','Helvetica','FontSize',30,'XTickLabelRotation',90);
xlabel('Size of periodic cubic network','FontName','Helvetica','FontSize',Constants.fontSize);
ylabel('Percolation threshold','FontName','Helvetica','FontSize',Constants.fontSize);
text(3.5,0.31,'Increasing overall size','FontName','Helvetica','FontSize',26,'HorizontalAlignment','center')
text(9.5,0.31,'Increasing width','FontName','Helvetica','FontSize',26,'HorizontalAlignment','center')
text(15.5,0.31,'Increasing length','FontName','Helvetica','FontSize',26,'HorizontalAlignment','center')
text(20.5,0.31,'Stacking network','FontName','Helvetica','FontSize',26,'HorizontalAlignment','center')
yticks(0.14:0.02:0.32);

%%

% Non-periodic networks
networkSizeNPList = [30,30,30; ...
    30,30,60; ...
    30,30,90; ...
    30,30,120];

% Number of regular networks non-periodic
nonPeriodicNum = size(networkSizeNPList,1);

% Actual networks
cuCoarseList = [12,15,20,30];
triangulationMethodList = Constants.triangulationMethod;

% Number of actual networks
actualNum = length(cuCoarseList)*length(triangulationMethodList);
actualActualNum = 3;

% Initialisation
p = NaN(periodicNum+nonPeriodicNum+actualNum*2,n);
nodeNum = NaN(periodicNum+nonPeriodicNum+actualNum*2,1);
edgeNum = nodeNum;
ratio = ones(periodicNum+nonPeriodicNum+actualNum*2,1);

% Determine the number of edges
nodeNum(1:periodicNum) = networkSizeList(:,1).*networkSizeList(:,2).*networkSizeList(:,3);
nodeNum(periodicNum+1:periodicNum+nonPeriodicNum) = networkSizeNPList(:,1).*networkSizeNPList(:,2).*networkSizeNPList(:,3);
nodeNum(periodicNum+nonPeriodicNum+1:periodicNum+nonPeriodicNum+actualNum) = [21030;9822;23483;10693;30091;13937;55503;25231];
edgeNum(1:periodicNum) = nodeNum(1:periodicNum)*3-(networkSizeList(:,1).*networkSizeList(:,2));
edgeNum(periodicNum+1:periodicNum+nonPeriodicNum) = nodeNum(periodicNum+1:periodicNum+nonPeriodicNum)*3-(networkSizeNPList(:,1).*networkSizeNPList(:,2))- ...
    (networkSizeNPList(:,1).*networkSizeNPList(:,3))-(networkSizeNPList(:,2).*networkSizeNPList(:,3));
edgeNum(periodicNum+nonPeriodicNum+1:periodicNum+nonPeriodicNum+actualNum) = [40542;26666;45305;29378;58154;37790;107922;69301];

% Determine the ratio of length to width
ratio(1:periodicNum) = networkSizeList(:,3)./networkSizeList(:,1);
ratio(periodicNum+1:periodicNum+nonPeriodicNum) = networkSizeNPList(:,3)./networkSizeNPList(:,1);

% Load data for regular networks
for i = 1:periodicNum
    data = load(FileLocation.percolationThreshold('cubic',[],[],[],networkSizeList(i,:),1));
    p(i,:) = data.p;
end

% Load data for regular networks without periodic boundary
for i = 1:nonPeriodicNum
    data = load(FileLocation.percolationThreshold('cubic',[],[],[],networkSizeNPList(i,:),0));
    p(periodicNum+i,:) = data.p;
end

% Load data for actual networks (randomly assigned properties)
k = 0;
for i = 1:length(cuCoarseList)
    for j = 1:length(triangulationMethodList)
        k = k+1;
        data = load(FileLocation.percolationThreshold('actual','random',cuCoarseList(i),triangulationMethodList{j},[],[]));
        p(periodicNum+nonPeriodicNum+k,:) = data.p;
    end
end

% Load data for actual networks
k = 0;
for i = 1:length(cuCoarseList)
    for j = 1:length(triangulationMethodList)
        k = k+1;
        data = load(FileLocation.percolationThreshold('actual','actual',cuCoarseList(i),triangulationMethodList{j},[],[]));
        p(periodicNum+nonPeriodicNum+actualNum+k,1:actualActualNum) = data.p;
    end
end

% Determine the mean and the std
pMean = mean(p,2);
pStd = std(p,[],2);
pMean(periodicNum+nonPeriodicNum+actualNum+1:end) = mean(p(periodicNum+nonPeriodicNum+actualNum+1:end,1:actualActualNum),2);
pStd(periodicNum+nonPeriodicNum+actualNum+1:end) = std(p(periodicNum+nonPeriodicNum+actualNum+1:end,1:actualActualNum),[],2);

% Determine average node connectivity
meanNodeConnectivity = zeros(periodicNum+nonPeriodicNum+actualNum,1);
for i = 1:periodicNum
    filter = Filter([],[],[],[],[],'cubic',[],[networkSizeList(i,1),networkSizeList(i,2),networkSizeList(i,3)],[],'undirected');
    filter = filter.generateNetworkTopologyRegular();
    meanNodeConnectivity(i) = mean(filter.nodeConnectivity);
end
for i = 1:nonPeriodicNum
    filter = Filter([],[],[],[],[],'cubicnp',[],[networkSizeNPList(i,1),networkSizeNPList(i,2),networkSizeNPList(i,3)],[],'undirected');
    filter = filter.generateNetworkTopologyRegular();
    meanNodeConnectivity(periodicNum+i) = mean(filter.nodeConnectivity);
end
k = 0;
for i = 1:length(cuCoarseList)
    for j = 1:length(triangulationMethodList)
        k = k+1;
        filter = Filter(1,cuCoarseList(i)/10,'dense',[0,0,0],triangulationMethodList{j},'actual','actual',[],'xpos','undirected');
        filter = filter.loadDTM('read');
        filter = filter.generateNetwork();
        meanNodeConnectivity(periodicNum+nonPeriodicNum+k) = mean(filter.nodeConnectivity);
    end
end

% Figure properties
legendList = {'Periodic cubic network: Increasing overall size','Periodic cubic network: Increasing width', ...
    'Periodic cubic network: Increasing length','Periodic cubic network: Stacking network', ...
    'Non-periodic cubic network: Stacking network','Actual network with size [1,1,1] (wDTM)','Actual network with size [1,1,1] (mDTM)'};

% Plot standatd deviation against no of nodes
figure;
set(gcf,'WindowState','maximized');
plot(nodeNum(1:6),pStd(1:6),'LineWidth',Constants.lineWidth,'Color',colourList(1,:));
hold on;
plot(nodeNum(7:12),pStd(7:12),'LineWidth',Constants.lineWidth,'Color',colourList(2,:));
plot(nodeNum(13:18),pStd(13:18),'LineWidth',Constants.lineWidth,'Color',colourList(3,:));
plot(nodeNum(19:22),pStd(19:22),'LineWidth',Constants.lineWidth,'Color',colourList(4,:));
plot(nodeNum(periodicNum+1:periodicNum+nonPeriodicNum),pStd(periodicNum+1:periodicNum+nonPeriodicNum),'LineWidth',Constants.lineWidth,'Color',colourList(5,:));

% Plot points
% plot(nodeNum(periodicNum+nonPeriodicNum+1:2:periodicNum+nonPeriodicNum+actualNum),pStd(periodicNum+nonPeriodicNum+1:2:periodicNum+nonPeriodicNum+actualNum),'LineStyle','none','Marker','^','MarkerSize',20,'LineWidth',0.1,'Color',colourList(6,:),'MarkerFaceColor',colourList(6,:));
% plot(nodeNum(periodicNum+nonPeriodicNum+2:2:periodicNum+nonPeriodicNum+actualNum),pStd(periodicNum+nonPeriodicNum+2:2:periodicNum+nonPeriodicNum+actualNum),'LineStyle','none','Marker','v','MarkerSize',20,'LineWidth',0.1,'Color',colourList(7,:),'MarkerFaceColor',colourList(7,:));
plot(nodeNum(1:periodicNum+nonPeriodicNum),pStd(1:periodicNum+nonPeriodicNum),'LineStyle','none','Marker','o','MarkerSize',10,'LineWidth',0.1,'Color','k','MarkerFaceColor','k');

% Plot limit
plot(20000*ones(1,2),ylim,'LineWidth',2,'LineStyle','--','Color','k');

xlim([0,2.6e5]);
xticks([0,0.2,0.5:0.5:4]*1e5)
set(gca, 'FontName','Helvetica','FontSize',Constants.fontSize);
xlabel('Number of nodes','FontName','Helvetica','FontSize',Constants.fontSize);
ylabel('Std of percolation threshold','FontName','Helvetica','FontSize',Constants.fontSize);
legend(legendList{1:5}, 'Location','northeast','FontName','Helvetica','FontSize',Constants.fontSize);

% Texts
text(nodeNum(1),pStd(1),'---(10,10,10)','FontName','Helvetica','FontSize',30,'HorizontalAlignment','left')
text(nodeNum(3),pStd(3),'---------------(30,30,30)','FontName','Helvetica','FontSize',30,'HorizontalAlignment','left')
text(nodeNum(6),pStd(6),'---(60,60,60)','FontName','Helvetica','FontSize',30,'HorizontalAlignment','left')

%%
% Legend list
legendList = {'Periodic cubic network: Increasing overall size','Periodic cubic network: Increasing width', ...
    'Periodic cubic network: Increasing length','Periodic cubic network: Stacking network', ...
    'Non-periodic cubic network: Stacking network'};

% Plot mean against the ratio of length to width
figure;
set(gcf,'WindowState','maximized');
plot(ratio(1:6),pMean(1:6),'LineWidth',Constants.lineWidth,'Color',colourList(1,:));
hold on;
plot(ratio(7:12),pMean(7:12),'LineWidth',Constants.lineWidth,'Color',colourList(2,:));
plot(ratio(13:18),pMean(13:18),'LineWidth',Constants.lineWidth,'Color',colourList(3,:));
plot(ratio(19:22),pMean(19:22),'LineWidth',Constants.lineWidth,'Color',colourList(4,:));
plot(ratio(periodicNum+1:periodicNum+nonPeriodicNum),pMean(periodicNum+1:periodicNum+nonPeriodicNum),'LineWidth',Constants.lineWidth,'Color',colourList(5,:));

% Plot points
% plot(ratio(periodicNum+nonPeriodicNum+1:2:periodicNum+nonPeriodicNum+actualNum),pMean(periodicNum+nonPeriodicNum+1:2:periodicNum+nonPeriodicNum+actualNum),'LineStyle','none','Marker','^','MarkerSize',10,'LineWidth',0.1,'Color',colourList(6,:),'MarkerFaceColor',colourList(6,:));
% plot(ratio(periodicNum+nonPeriodicNum+2:2:periodicNum+nonPeriodicNum+actualNum),pMean(periodicNum+nonPeriodicNum+2:2:periodicNum+nonPeriodicNum+actualNum),'LineStyle','none','Marker','v','MarkerSize',10,'LineWidth',0.1,'Color',colourList(7,:),'MarkerFaceColor',colourList(7,:));
plot(ratio(1:periodicNum+nonPeriodicNum),pMean(1:periodicNum+nonPeriodicNum),'LineStyle','none','Marker','o','MarkerSize',10,'LineWidth',0.1,'Color','k','MarkerFaceColor','k');
xlim([-0.4,5.4]);
ylim([0.2,0.28]);
plot(xlim,[0.2488,0.2488],'LineWidth',Constants.lineWidth,'LineStyle','--','Color','k');
set(gca, 'FontName','Helvetica','FontSize',Constants.fontSize);
xlabel('Ratio of network length to network width','FontName','Helvetica','FontSize',Constants.fontSize);
ylabel('Mean of percolation threshold','FontName','Helvetica','FontSize',Constants.fontSize);
legend(legendList, 'Location','southeast','FontName','Helvetica','FontSize',Constants.fontSize);

% Texts
text(ratio(3),pMean(3),'(30,30,30)--------','FontName','Helvetica','FontSize',30,'HorizontalAlignment','right')
text(ratio(7),pMean(7),'---(10,10,30)','FontName','Helvetica','FontSize',30,'HorizontalAlignment','left')
text(ratio(12),pMean(12),'(60,60,30)--','FontName','Helvetica','FontSize',30,'HorizontalAlignment','right')
text(ratio(13),pMean(13),'---(30,30,10)','FontName','Helvetica','FontSize',30,'HorizontalAlignment','left')
text(ratio(22),pMean(22),'---(30,30,120)','FontName','Helvetica','FontSize',30,'HorizontalAlignment','left')
text(ratio(23),pMean(23),'(30,30,30)(NP)---','FontName','Helvetica','FontSize',30,'HorizontalAlignment','right')
text(ratio(26),pMean(26),'---(30,30,120)(NP)','FontName','Helvetica','FontSize',30,'HorizontalAlignment','left')
text(3.5,0.2515,' Infinitely large cubic network','FontName','Helvetica','FontSize',30,'HorizontalAlignment','center')

%%
% Legend
legendList = {'Periodic cubic network','Non-periodic cubic network', ...
    'Actual network with wDTM randomly distributed','Actual network with mDTM randomly distributed', ...
    'Actual network with wDTM','Actual network with mDTM', ...
    'Infinitely large square network','Infinitely large diamond network','Infinitely large cubic network'};

% Plot mean against mean of connectivity
figure;
set(gcf,'WindowState','maximized');
plot(meanNodeConnectivity(1:periodicNum),pMean(1:periodicNum),'LineStyle','none','Marker','o','MarkerSize',20,'LineWidth',0.1,'Color',colourList(1,:),'MarkerFaceColor',colourList(1,:));
hold on;
plot(meanNodeConnectivity(periodicNum+1:periodicNum+nonPeriodicNum),pMean(periodicNum+1:periodicNum+nonPeriodicNum),'LineStyle','none','Marker','s','MarkerSize',20,'LineWidth',0.1,'Color',colourList(5,:),'MarkerFaceColor',colourList(5,:));
plot(meanNodeConnectivity(periodicNum+nonPeriodicNum+1:2:periodicNum+nonPeriodicNum+actualNum),pMean(periodicNum+nonPeriodicNum+1:2:periodicNum+nonPeriodicNum+actualNum),'LineStyle','none','Marker','^','MarkerSize',20,'LineWidth',0.1,'Color',colourList(6,:),'MarkerFaceColor',colourList(6,:));
plot(meanNodeConnectivity(periodicNum+nonPeriodicNum+2:2:periodicNum+nonPeriodicNum+actualNum),pMean(periodicNum+nonPeriodicNum+2:2:periodicNum+nonPeriodicNum+actualNum),'LineStyle','none','Marker','v','MarkerSize',20,'LineWidth',0.1,'Color',colourList(7,:),'MarkerFaceColor',colourList(7,:));
plot(meanNodeConnectivity(periodicNum+nonPeriodicNum+1:2:periodicNum+nonPeriodicNum+actualNum),pMean(periodicNum+nonPeriodicNum+actualNum+1:2:periodicNum+nonPeriodicNum+actualNum*2),'LineStyle','none','Marker','^','MarkerSize',20,'LineWidth',0.1,'Color',colourList(2,:),'MarkerFaceColor',colourList(2,:));
plot(meanNodeConnectivity(periodicNum+nonPeriodicNum+2:2:periodicNum+nonPeriodicNum+actualNum),pMean(periodicNum+nonPeriodicNum+actualNum+2:2:periodicNum+nonPeriodicNum+actualNum*2),'LineStyle','none','Marker','v','MarkerSize',20,'LineWidth',0.1,'Color',colourList(3,:),'MarkerFaceColor',colourList(3,:));
xlim([3.8,6.05]);
ylim([0.18,0.65]);
plot(4,0.5,'LineWidth',Constants.lineWidth,'Marker','x','MarkerSize',30,'LineStyle','none','Color',colourList(1,:));
plot(4,0.388,'LineWidth',Constants.lineWidth,'Marker','x','MarkerSize',30,'LineStyle','none','Color',colourList(2,:));
plot(6,0.2488,'LineWidth',Constants.lineWidth,'Marker','x','MarkerSize',30,'LineStyle','none','Color',colourList(3,:));
set(gca, 'FontName','Helvetica','FontSize',Constants.fontSize);
xlabel('Mean of node connectivity','FontName','Helvetica','FontSize',Constants.fontSize);
ylabel('Mean of percolation threshold','FontName','Helvetica','FontSize',Constants.fontSize);
legend(legendList, 'Location','northeast','FontName','Helvetica','FontSize',Constants.fontSize);