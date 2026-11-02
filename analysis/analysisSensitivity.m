% The script plots graphs relevant to show the sensitivity of results to
% gravity and node volume

clear
clc
close all

% Add src to path
path(path, ['..', filesep, 'src'])

% Load the default
simulator = Simulator.loadResults(FileLocation.simulatorMonodisperse('flowbiased', ...
    'actual','actual',1.2,'dense',[1,1,1],'Del'));
d.d1 = simulator.monodisperseBaseDiameter;
p.p1 = simulator.percentageErodedNum;
r.r1 = simulator.directness;
pAffected.p1 = simulator.percentageNodeAffected;

% Sensitivity to gravity
simulator = Simulator.loadResults(FileLocation.simulatorMonodisperse('flowbiasedpos', ...
    'actual','actual',1.2,'dense',[1,1,1],'Del'));
d.d2 = simulator.monodisperseBaseDiameter;
p.p2 = simulator.percentageErodedNum;
r.r2 = simulator.directness;
pAffected.p2 = simulator.percentageNodeAffected;

% Sensitivity to node volume
simulator = Simulator.loadResults(FileLocation.simulatorMonodisperse('flowbiasednode', ...
    'actual','actual',1.2,'dense',[1,1,1],'Del'));
d.d3 = simulator.monodisperseBaseDiameter;
p.p3 = simulator.percentageErodedNum;
r.r3 = simulator.directness;
pAffected.p3 = simulator.percentageNodeAffected;

% Combination of both
simulator = Simulator.loadResults(FileLocation.simulatorMonodisperse('flowbiasedextreme', ...
    'actual','actual',1.2,'dense',[1,1,1],'Del'));
d.d4 = simulator.monodisperseBaseDiameter;
p.p4 = simulator.percentageErodedNum;
r.r4 = simulator.directness;
pAffected.p4 = simulator.percentageNodeAffected;

%%

% Figure properties
colourList = Constants.colourList;
lineStyleList = {'-','--',':','-.'};
legendList = {'Default parameters','(1) Restrict to move in positive z-direction only','(2) Accelerated node clogging','(3) Combination of both'};

% Plot percentage eroded by num against d
figure;
set(gcf,'WindowState','maximized');
for i = 1:4
    plot(d.(sprintf('d%i',i)),p.(sprintf('p%i',i)),'LineWidth',4,'LineStyle',lineStyleList{i});
    hold on;
end
hold off;
% plot((d.d1),((p.p1/100).^2)*100,'LineWidth',5,'LineStyle',lineStyleList{(i>2)+1});
% legend(legendList,'Location','southwest','FontName','Helvetica','FontSize',36);
xlim([min(d.d1),max(d.d1)]);
set(gca, 'FontName','Helvetica','FontSize',36);
xlabel('Normalised diameter of monodisperse base particles [-]','FontName','Helvetica','FontSize',36);
ylabel('Percentage of base particles eroded','FontName','Helvetica','FontSize',36);

% Plot percentage of nodes affected by node volume against d
figure;
set(gcf,'WindowState','maximized');
for i = 1:4
    plot(d.(sprintf('d%i',i)),pAffected.(sprintf('p%i',i)),'LineWidth',4,'LineStyle',lineStyleList{i});
    hold on;
end
legend(legendList,'Location','northeast','FontName','Helvetica','FontSize',36);
xlim([min(d.d1),max(d.d1)]);
ylim([0,100])
set(gca, 'FontName','Helvetica','FontSize',36);
xlabel('Normalised diameter of monodisperse base particles [-]','FontName','Helvetica','FontSize',36);
ylabel('Percentage of nodes affected','FontName','Helvetica','FontSize',36);