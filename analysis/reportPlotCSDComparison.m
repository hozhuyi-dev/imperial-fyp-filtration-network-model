% The script compares the CSD obtained from the wDTM and from analytical
% method

clear
clc
close all

% Add src to path
path(path, ['..',  filesep, 'src'])

% Filters to consider
cuCoarseList = 1.2;
densityList = Constants.density;
RdList = [1,0.75,0.5,0.25,0];

% Initialise
simulator = Simulator;

% Plot properties
colourList = Constants.colourList;

% Legend list
legendList = {'DEM-DTM(wDTM) Loose','DEM-DTM(wDTM) Medium','DEM-DTM(wDTM) Dense', ...
    'DEM-DTM(mDTM) Dense', ...
    'ICM (3 Particles) R_d = 100%','ICM (3 Particles) R_d = 75%', ...
    'ICM (3 Particles) R_d = 50%','ICM (3 Particles) R_d = 25%', ...
    'ICM (3 Particles) R_d = 0%','ICM (4 Particles)'};

% For each cuCoarse
for i = 1:length(cuCoarseList)
    % Initialise figure
    fig = figure;
    set(gcf,'WindowState','maximized')
    
    % For each density
    for j = 1:length(densityList)
        % Load filter
        simulator = simulator.createFilter(cuCoarseList(i),'density',densityList{j});
        
        % CSD from DEM-wDTM
        csd = simulator.filter.csdNum;
        
        % Plot
        semilogx(csd(:,2),csd(:,1),'LineWidth',3)
        hold on
    end
    for j = 3
        % Load filter
        simulator = simulator.createFilter(cuCoarseList(i),'density',densityList{j},'triangulationMethod','Mod','size',[0,0,0]);
        
        % CSD from DEM-wDTM
        csd = simulator.filter.csdNum;
        
        % Plot
        semilogx(csd(:,2),csd(:,1),'LineWidth',3,'LineStyle',':')
    end
    
    % For each Rd
    for j = 1:length(RdList)
        csd = simulator.filter.analyticalCSDCombined(3,RdList(j));
        semilogx(csd(:,2),csd(:,1),'LineWidth',3,'LineStyle','--','Color',colourList(j,:))
    end
    
    % Loosest analytical csd
    csd = simulator.filter.analyticalCSDCombined(4,0);
    semilogx(csd(:,2),csd(:,1),'LineWidth',3,'LineStyle','--','Color',colourList(7,:))
    
    % Figure labels
    plot(xlim,35*ones(1,2),'LineStyle','--','LineWidth',3,'Color','k')
    legend(legendList,'Location','southeast','FontName','Helvetica','FontSize',18)
    set(gca,'FontName','Helvetica','FontSize',18)
    xlabel('Constriction size normalised by D_0_F [-]','FontName','Helvetica','FontSize',18)
    ylabel('Percentage smaller','FontName','Helvetica','FontSize',18)
    
    print(fig,[FileLocation.dataFolder, 'figures', filesep, 'CSDComparison'],'-dpng','-r300')
end