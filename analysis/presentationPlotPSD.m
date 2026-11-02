% The script plots the PSD of a base with d0B of 0.25 with cu = 1.2

clear
clc
close all

path(path, '../src')

simulator = Simulator;
simulator = simulator.createBase(0.3);

figure
set(gcf,'WindowState','maximized')
semilogx(simulator.base.base1.psdFineVol(:,2),simulator.base.base1.psdFineVol(:,1),'LineWidth',5,'Color','k')
xlim([0.2,0.5])
set(gca, 'FontName','Helvetica','FontSize',36);
xlabel('Diameter of base particles','FontName','Helvetica','FontSize',36);
ylabel('Percentage smaller','FontName','Helvetica','FontSize',36);