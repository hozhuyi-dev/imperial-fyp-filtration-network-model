% The script plots the psd by volume and psd by surface area

clear
clc
close all

path(path, '../src')

simulator = Simulator;
simulator = simulator.createBase(0.3,'cuFine',4.5);

psdFineSA = Statistics.convertPSDVolToArea(simulator.base.base1.psdFineVol);

figure
set(gcf,'WindowState','maximized')
semilogx(simulator.base.base1.psdFineVol(:,2),simulator.base.base1.psdFineVol(:,1),'LineWidth',5,'Color','k')
hold on
semilogx(psdFineSA(:,2),psdFineSA(:,1),'LineWidth',5,'Color','k')
set(gca, 'FontName','Helvetica','FontSize',36);
xlabel('Diameter of base particles','FontName','Helvetica','FontSize',36);
ylabel('Percentage smaller','FontName','Helvetica','FontSize',36);