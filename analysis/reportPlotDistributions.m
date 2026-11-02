% The script generates plots of geometric properties distributions

clear;
clc;
close all;

% Add src/ to path
path(path, ['..', filesep, 'src']);

% Data available
cuCoarseList = Constants.cuCoarse;
densityList = Constants.density;
size = [1,1,1];

% Initialise struct
filterDel = struct;
filterMod = struct;

% Load data for Del
k = 0;
for i = 1:length(cuCoarseList)
    for j = 1:length(densityList)
        % Current iteration
        k = k+1;
        
        % Create filter object
        filterDel.(sprintf('filter%i',k)) = Filter(1,cuCoarseList(i),densityList{j},size,'Del','actual','actual',[],'zpos','undirected');
        filterDel.(sprintf('filter%i',k)) = filterDel.(sprintf('filter%i',k)).loadDTM('read');
    end
end

% Load data for Mod
k = 0;
for i = 1:length(cuCoarseList)
    % Current iteration
    k = k+1;
    
    % Create filter object
    simulator = Simulator;
    simulator = simulator.createFilter(cuCoarseList(i),'density','dense','size',[0,0,0],'triangulationMethod','Mod');
    filterMod.(sprintf('filter%i',k)) = simulator.filter;
end

% Load Shire's data
% csdNum = struct;
% for i = 1:length(cuCoarseList)
%     [csdNum.(sprintf('csdNum%i',i)), ~, ~, ~, ~, ~, ~] = Validate.generateInputShire(1.2, 0.5, ...
%         sprintf('Cu%ipt%i',floor(cuCoarseList(i)),round((cuCoarseList(i)-floor(cuCoarseList(i)))*10)), 'dense');
% end

% Plot properties
colourList = Constants.colourList;
lineStyleList = {'-','--',':'};
triangulationMethodLegendList = {'wDTM','mDTM'};

%% Filter PSD

% Legend list
legendList = cell(length(cuCoarseList),1);
for i = 1:length(cuCoarseList)
    legendList{i} = sprintf('C_u = %.1f', cuCoarseList(i));
end

% Plot
figure;
set(gcf,'WindowState','maximized');
for i = 1:length(cuCoarseList)
    idx = i*length(densityList)-1;
    
    % Create CDF
    CDF = Statistics.convertDataToPSDVol(filterDel.(sprintf('filter%i',idx)).particleRadius*2);
    
    % Plot
    semilogx(CDF(:,2),CDF(:,1),'LineWidth',Constants.lineWidth);
    hold on;
end
legend(legendList,'Location','southeast','FontName','Helvetica','FontSize',Constants.fontSize);
set(gca,'FontName','Helvetica','FontSize',Constants.fontSize);
xlabel('Normalised filter particle diameter (D_F/D_0_F)','FontName','Helvetica','FontSize',Constants.fontSize);
ylabel('Percentage smaller [%]','FontName','Helvetica','FontSize',Constants.fontSize);

%% Other properties

% Inputs
property = 'edgeDiameter';

% Legend list
% legendList = cell(length(cuCoarseList)*length(densityList),1);
% k = 0;
% for i = 1:length(cuCoarseList)
%     for j = 1:length(densityList)
%         k = k+1;
%         legendList{k} = sprintf('C_u = %.1f Density = %s', cuCoarseList(i), densityList{j});
%     end
% end
legendList = {'C_u = 1.2 Density = dense (mDTM)','C_u = 1.5 Density = dense (mDTM)','C_u = 2.0 Density = dense (mDTM)','C_u = 3.0 Density = dense (mDTM)', ...
    'C_u = 1.2 Density = dense (wDTM)','C_u = 1.5 Density = dense (wDTM)','C_u = 2.0 Density = dense (wDTM)','C_u = 3.0 Density = dense (wDTM)'};

% Plot
figure;
set(gcf,'WindowState','maximized');
% k = 0;
% for i = 1:length(cuCoarseList)
%     k = k+1;
%     
%     % Create CDF
%     CDF = Statistics.createCDF(filterMod.(sprintf('filter%i',k)).(sprintf('%s',property)),1e3);
%     
%     % Plot
%     semilogx(CDF(:,2),CDF(:,1),'LineWidth',Constants.lineWidth,'LineStyle','-.','Color',colourList(i,:));
%     hold on
% end
k = 0;
for i = 1%:length(cuCoarseList)
    for j = 3%1:length(densityList)
        k = k+1;
        
        if j==3
            % Create CDF
            CDF = Statistics.createCDF(filterDel.(sprintf('filter%i',k)).(sprintf('%s',property)),1e3);
            
            % Plot
            semilogx(CDF(:,2),CDF(:,1),'LineWidth',Constants.lineWidth,'LineStyle','-','Color','k')%colourList(i,:));
        end
    end
end
xlim([0.158,1])
ylim([0,100])
xticks([0.05,0.1,0.5,1,5,10])
% xticks(0.1:0.1:0.8)
% legend(legendList,'Location','southeast','FontName','Helvetica','FontSize',26);
set(gca,'FontName','Helvetica','FontSize',Constants.fontSize);
xlabel('Normalised constriction size (D_c/D_0_F)','FontName','Helvetica','FontSize',Constants.fontSize);
ylabel('Percentage smaller','FontName','Helvetica','FontSize',Constants.fontSize);

%% Edge diameter normalised by D5F or D15F

% Legend list
legendList = cell(length(cuCoarseList)*length(densityList),1);
k = 0;
for i = 1:length(cuCoarseList)
    for j = 1:length(densityList)
        k = k+1;
        legendList{k} = sprintf('C_u = %.1f %s', cuCoarseList(i), densityList{j});
    end
end

% Plot 1 (D5F)
figure;
set(gcf,'WindowState','maximized');
k = 0;
for i = 1:length(cuCoarseList)
    for j = 1:length(densityList)
        k = k+1;
        
        % Normalising factor
        CDF = Statistics.convertDataToPSDVol(filterDel.(sprintf('filter%i',k)).particleRadius*2);
        D5F = interp1(CDF(:,1),CDF(:,2),5);
        
        % Create CDF
        CDF = Statistics.createCDF(filterDel.(sprintf('filter%i',k)).edgeDiameter/D5F,1e3);
        
        % Plot
        semilogx(CDF(:,2),CDF(:,1),'LineWidth',Constants.lineWidth,'LineStyle',lineStyleList{j},'Color',colourList(i,:));
        hold on;
    end
end
semilogx([0.25,0.25],[0,100],'Color','k','LineStyle','--','LineWidth',1.5);
legend(legendList,'Location','southeast','FontName','Helvetica','FontSize',Constants.fontSize);
set(gca,'FontName','Helvetica','FontSize',Constants.fontSize);
xlabel('Normalised constriction diameter (D_c/D_5_F)','FontName','Helvetica','FontSize',Constants.fontSize);
ylabel('Percentage smaller','FontName','Helvetica','FontSize',Constants.fontSize);
xticks([0.25,0.5,1]);
xticklabels({'0.25','0.50','1.00'});
text(0.173,80,'Kenney et. al. (1985):','FontName','Helvetica','FontSize',14)
text(0.173,76.5,'D_c^* = 0.25 D_5_F','FontName','Helvetica','FontSize',14)

% Plot 2 (D15F)
figure;
set(gcf,'WindowState','maximized');
k = 0;
for i = 1:length(cuCoarseList)
    for j = 1:length(densityList)
        k = k+1;
        
        % Normalising factor
        CDF = Statistics.convertDataToPSDVol(filterDel.(sprintf('filter%i',k)).particleRadius*2);
        D15F = interp1(CDF(:,1),CDF(:,2),15);
        
        % Create CDF
        CDF = Statistics.createCDF(filterDel.(sprintf('filter%i',k)).edgeDiameter/D15F,1e3);
        
        % Plot
        semilogx(CDF(:,2),CDF(:,1),'LineWidth',Constants.lineWidth,'LineStyle',lineStyleList{j},'Color',colourList(i,:));
        hold on;
    end
end
semilogx([0.2,0.2],[0,100],'Color','k','LineStyle','--','LineWidth',1.5);
legend(legendList,'Location','southeast','FontName','Helvetica','FontSize',Constants.fontSize);
set(gca,'FontName','Helvetica','FontSize',Constants.fontSize);
xlabel('Normalised constriction diameter (D_c/D_15_F)','FontName','Helvetica','FontSize',Constants.fontSize);
ylabel('Percentage smaller','FontName','Helvetica','FontSize',Constants.fontSize);
xticks([0.2,0.25,0.5,1]);
xticklabels({'0.20','0.25','0.50','1.00'});
text(0.14,80,'Kenney et. al. (1985):','FontName','Helvetica','FontSize',14)
text(0.14,76.5,'D_c^* = 0.2 D_1_5_F','FontName','Helvetica','FontSize',14)