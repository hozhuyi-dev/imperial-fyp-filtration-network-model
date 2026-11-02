% The script is used to evaluate the effects of network type and assignment
% method on the largest particle eroded and percentage of particles eroded

clear;
clc;
close all;

% Add src to path
path(path, ['..', filesep, 'src']);

%% Load data

% Networks to consider
networkList = {'flowbiased','actual','actual',1.2,'dense',[1,1,1],'Del'; ...
    'flowbiased','actual','actual',1.2,'dense',[0,0,0],'Mod'; ...
    'flowbiased','actual','random',1.2,'dense',[1,1,1],'Del'; ...
    'flowbiased','squarenp','random',1.2,'dense',[1,1,1],'Del'; ...
    'flowbiased','slantedsquarenp','random',1.2,'dense',[1,1,1],'Del'; ...
    'flowbiased','cubicnp','random',1.2,'dense',[1,1,1],'Del'; ...
    'flowbiased','slantedcubicnp','random',1.2,'dense',[1,1,1],'Del'; ...
    'flowbiased','bccnp','random',1.2,'dense',[1,1,1],'Del'; ...
    'flowbiased','fccnp','random',1.2,'dense',[1,1,1],'Del'};
% networkList = {'flowbiased','actual','actual',1.2,'dense',[1,1,1],'Del'; ...
%     'flowbiasedpos','actual','actual',1.2,'dense',[1,1,1],'Del'};
% networkList = {'flowbiased','actual','actual',3,'loose',[1,1,1],'Del'; ...
%     'flowbiased','actual','random',3,'loose',[1,1,1],'Del'; ...
%     'flowbiased','squarenp','random',3,'loose',[1,1,1],'Del'; ...
%     'flowbiased','slantedsquarenp','random',3,'loose',[1,1,1],'Del'; ...
%     'flowbiased','cubicnp','random',3,'loose',[1,1,1],'Del'; ...
%     'flowbiased','slantedcubicnp','random',3,'loose',[1,1,1],'Del'; ...
%     'flowbiased','bccnp','random',3,'loose',[1,1,1],'Del'; ...
%     'flowbiased','fccnp','random',3,'loose',[1,1,1],'Del'};
networkNum = size(networkList,1);

% Initialise array(s)
p = struct;
d = struct;
r = struct;
nEroded = struct;
nRetained = struct;
pAffected = struct;
mClusterSize = struct;
mSmaller = struct;

% Load data
for i = 1:networkNum
    if i==1
        endWith = '';
    else
        endWith = '';
    end
    
    % Load Simulator object
    simulator = Simulator.loadResults(FileLocation.simulatorMonodisperse(networkList{i,1},networkList{i,2}, ...
        networkList{i,3},networkList{i,4},networkList{i,5},networkList{i,6},networkList{i,7},endWith));
    
    % Determine the diameter of particles
    d.(sprintf('d%i',i)) = simulator.monodisperseBaseDiameter;
    
    % Determine the percentage of particles eroded
    p.(sprintf('p%i',i)) = simulator.percentageErodedNum;
    
    % Determine the ratio of average distance travelled to shortest path
    if i==1
        [~,shortestD] = simulator.filter.shortestPath;
    end
    if i==2
        r.(sprintf('r%i',i)) = simulator.directness('d',shortestD);
    else
        r.(sprintf('r%i',i)) = simulator.directness;
    end
    
    % Determine the average number of constrictions met
    nEroded.(sprintf('n%i',i)) = simulator.constrictionNum;
    nRetained.(sprintf('n%i',i)) = simulator.constrictionNum('mechanism','retained');
    
    % Determine the percentage of nodes affected
    pAffected.(sprintf('p%i',i)) = simulator.percentageNodeAffected;
    
    % Determine the mean cluster size
    mClusterSize.(sprintf('m%i',i)) = simulator.meanClusterSize;
    mSmaller.(sprintf('m%i',i)) = simulator.meanClusterEdgeDiameter;
    
    % Save shortest path
    % simulator.saveResults(simulator.fileLocation);
    
    fprintf('Loaded %6.2f%%\n', i/networkNum*100);
end
clear s;

%% Figures

% Figure properties
colourList = Constants.colourList;
lineStyleList = {'-','--'};
legendList = {'Type 3 wDTM (22437)','Type 3 mDTM (9822)','Type 2 (22437)', ...
    'Type 1 Simple Square (24389)','Type 1 Slanted Square (24389)', ...
    'Type 1 Simple Cubic (24389)','Type 1 Slanted Cubic (26100)', ...
    'Type 1 BCC (25578)','Type 1 FCC (24389)'};
% legendList = {'Flow-biased probability','Flow-biased probability and positive z-direction only'};

%%
% Plot percentage eroded by num against d
figure;
set(gcf,'WindowState','maximized');
for i = 1:networkNum
    if i<=7
        colourIdx = i;
    else
        colourIdx = i-7;
    end
    plot(d.(sprintf('d%i',i)),p.(sprintf('p%i',i)),'LineWidth',5,'LineStyle',lineStyleList{(i>3)+1},'Color',colourList(colourIdx,:),'Marker','none');
    hold on;
end
hold off;
% plot((d.d1),((p.p1/100).^2)*100,'LineWidth',5,'LineStyle',lineStyleList{(i>2)+1});
legend(legendList,'Location','northeast','FontName','Helvetica','FontSize',26);
xlim([min(d.d1),max(d.d1)]);
set(gca, 'FontName','Helvetica','FontSize',36);
xlabel('Normalised diameter of monodisperse base particles [-]','FontName','Helvetica','FontSize',36);
ylabel('Percentage of base particles eroded','FontName','Helvetica','FontSize',36);

%%
% yyaxis right
% for i = 1:networkNum
%     if i<=7
%         colourIdx = i;
%     else
%         colourIdx = i-7;
%     end
%     plot(d.(sprintf('d%i',i)),pAffected.(sprintf('p%i',i)),'LineWidth',5,'LineStyle',lineStyleList{(i>2)+1},'Color',colourList(colourIdx,:),'Marker','none');
%     hold on;
% end
% hold off;
% ylabel('Percentage of nodes affected','FontName','Helvetica','FontSize',36);

% Plot directness ratio against d
figure;
set(gcf,'WindowState','maximized');
for i = 1
    plot(d.(sprintf('d%i',i)),r.(sprintf('r%i',i)),'LineWidth',5,'LineStyle',lineStyleList{(i>3)+1});
    hold on;
end
% legend(legendList,'Location','northwest','FontName','Helvetica','FontSize',36);
xlim([min(d.d1),max(d.d1)]);
ylim([1.2,1.6])
set(gca, 'FontName','Helvetica','FontSize',36);
xlabel('Normalised diameter of monodisperse base particles [-]','FontName','Helvetica','FontSize',36);
ylabel('Normalised average distance travelled','FontName','Helvetica','FontSize',36);

%%

% Plot directness ratio against d
figure
set(gcf,'WindowState','maximized')
plot(1:5,[1.422,1.495,1.517,1.53,1.537],'LineWidth',5,'LineStyle',lineStyleList{(i>3)+1})
% legend(legendList,'Location','northwest','FontName','Helvetica','FontSize',36);
xlim([0,6])
ylim([1.2,1.6])
set(gca, 'FontName','Helvetica','FontSize',36)
xlabel('Simulation time','FontName','Helvetica','FontSize',36)
ylabel('Normalised average distance travelled','FontName','Helvetica','FontSize',36)

%%
% Plot number of constrictions met by eroded particles against d
figure;
set(gcf,'WindowState','maximized');
for i = 1:networkNum
    plot(d.(sprintf('d%i',i)),nEroded.(sprintf('n%i',i)),'LineWidth',5,'LineStyle',lineStyleList{(i>3)+1});
    hold on;
end
% legend(legendList,'Location','northwest','FontName','Helvetica','FontSize',36);
xlim([min(d.d1),max(d.d1)]);
set(gca, 'FontName','Helvetica','FontSize',36);
xlabel('Normalised diameter of monodisperse base particles [-]','FontName','Helvetica','FontSize',36);
ylabel('Average number of constrictions','FontName','Helvetica','FontSize',36);

% Plot number of constrictions met by retained particles against d
figure;
set(gcf,'WindowState','maximized');
for i = 1:networkNum
    plot(d.(sprintf('d%i',i)),nRetained.(sprintf('n%i',i)),'LineWidth',5,'LineStyle',lineStyleList{(i>3)+1});
    hold on;
end
% legend(legendList,'Location','northwest','FontName','Helvetica','FontSize',36);
xlim([min(d.d1),max(d.d1)]);
set(gca, 'FontName','Helvetica','FontSize',36);
xlabel('Normalised diameter of monodisperse base particles [-]','FontName','Helvetica','FontSize',36);
ylabel('Average number of constrictions','FontName','Helvetica','FontSize',36);

% Plot percentage of nodes affected by node volume against d
figure;
set(gcf,'WindowState','maximized');
for i = 1:networkNum
    plot(d.(sprintf('d%i',i)),pAffected.(sprintf('p%i',i)),'LineWidth',5,'LineStyle',lineStyleList{(i>3)+1});
    hold on;
end
legend(legendList,'Location','northwest','FontName','Helvetica','FontSize',36);
xlim([min(d.d1),max(d.d1)]);
set(gca, 'FontName','Helvetica','FontSize',36);
xlabel('Normalised diameter of monodisperse base particles [-]','FontName','Helvetica','FontSize',36);
ylabel('Percentage of nodes affected','FontName','Helvetica','FontSize',36);

% Plot mean cluster size of full nodes against d
figure;
set(gcf,'WindowState','maximized');
for i = [1,3]
    plot(d.(sprintf('d%i',i)),mClusterSize.(sprintf('m%i',i)),'LineWidth',5,'LineStyle',lineStyleList{(i>3)+1},'Color',colourList(i,:));
    hold on;
end
legendList = {'Type 3 wDTM (22437)','Type 2 (22437)'};
legend(legendList,'Location','northwest','FontName','Helvetica','FontSize',36);
xlim([min(d.d1),max(d.d1)]);
set(gca, 'FontName','Helvetica','FontSize',36);
xlabel('Normalised diameter of monodisperse base particles [-]','FontName','Helvetica','FontSize',36);
ylabel('Mean cluster size','FontName','Helvetica','FontSize',36);

% Plot mean cluster size of small nodes against d
figure;
set(gcf,'WindowState','maximized');
for i = 1:networkNum
    plot(d.(sprintf('d%i',i)),mSmaller.(sprintf('m%i',i)),'LineWidth',5,'LineStyle',lineStyleList{(i>3)+1});
    hold on;
end
% legend(legendList,'Location','northwest','FontName','Helvetica','FontSize',36);
xlim([min(d.d1),max(d.d1)]);
set(gca, 'FontName','Helvetica','FontSize',36);
xlabel('Normalised diameter of monodisperse base particles [-]','FontName','Helvetica','FontSize',36);
ylabel('Mean cluster size','FontName','Helvetica','FontSize',36);

% Save figure
% print(fig,[FileLocation.dataFolder, 'figures', filesep, 'Monodisperse12dense111'],'-dpng','-r300')