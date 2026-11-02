% The script analyses polydisperse simulations

clear
clc
close all

% Add src to path
path(path, ['..', filesep, 'src'])

% Simulation to consider
transportLogicType = 'flowbiased';
networkType = 'actual';
assignmentMethod = 'actual';
cuCoarseList = Constants.cuCoarse;
densityList = Constants.density;
siz = [1,1,1];
triangulationMethod = 'Del';
cuFineList = [1.5,3,4.5];
rList = 0.5:0.5:7;

% Initialise arrays
arrayRow = length(cuFineList)*length(rList);
arrayColumn = length(cuCoarseList)*length(densityList);
markerList = cell(arrayRow,arrayColumn);
ratio1 = zeros(arrayRow,arrayColumn);
ratio2 = ratio1;
ratio3 = ratio1;
ratio4 = ratio1;
ratio5 = ratio1;
ratio6 = ratio1;
normalisedVolumeEroded = ratio1;
category = ratio1;

% Plot properties
markers = {'o','x'};

fprintf('Loading data ')
idxColumn = 0;
for i = 1:length(cuCoarseList)
    for j = 1:length(densityList)
        idxColumn = idxColumn+1;
        
        idxRow = 0;
        for k = 1:length(cuFineList)
            for l = 1:length(rList)
                idxRow = idxRow+1;
                
                % Load simulator
                simulator = Simulator.loadResults(FileLocation.simulatorPolydisperseRepeat(transportLogicType, ...
                    networkType, assignmentMethod, cuCoarseList(i), densityList{j}, siz, triangulationMethod, ...
                    cuFineList(k), rList(l), 'continuous'));
                
                % D15F/d85B
                r1 = simulator.terzaghiRatio;
                r2 = simulator.fosterFellRatio;
                r3 = simulator.testRatio1;
                r4 = simulator.testRatio2;
                r5 = simulator.testRatio3;
                r6 = simulator.indraratnaRautRatio;
                ratio1(idxRow,idxColumn) = r1(1);
                ratio2(idxRow,idxColumn) = r2(1);
                ratio3(idxRow,idxColumn) = r3(1);
                ratio4(idxRow,idxColumn) = r4(1);
                ratio5(idxRow,idxColumn) = r5(1);
                ratio6(idxRow,idxColumn) = r6(1);
                
                % Volume of particles eroded normalised by total volume of
                % void space
                normalisedVolumeEroded(idxRow,idxColumn) = simulator.totalErodedVolArea;
                
                % Category of effectiveness
                bool = simulator.filterSeal;
                category(idxRow,idxColumn) = bool(end)+1;
                
                % Plot properties
                markerList{idxRow,idxColumn} = markers{category(idxRow,idxColumn)};
            end
        end
        
        % Display progress
        if idxColumn~=1
            fprintf(repmat('\b', 1, 7));
        end
        fprintf('%6.2f%%', idxColumn/arrayColumn*100);
    end
end
fprintf('\n')

% Finalise
normalisedVolumeEroded(normalisedVolumeEroded==0) = 1e-14;

%%

close all

% Legend list
legendList = {'Filter C_u = 1.2','Filter C_u = 1.5','Filter C_u = 2.0','Filter C_u = 3.0', ...
    'Base C_u = 1.5','Base C_u = 3.0','Base C_u = 4.5', ...
    'Seal','No seal'};
% legendList = {'Density = loose','Density = medium','Density = dense', ...
%     'Base C_u = 1.5','Base C_u = 3.0','Base C_u = 4.5', ...
%     'Seal','No seal'};

% Options
yl = [1e-3,1e2];
boundaryY = 0.1;
boundaryX1 = 4;
boundaryX2 = 1;
yString = 'Normalised total volume eroded';

% Plot properties
lineStyleList = {'-','--',':'};
colourList = Constants.colourList;

% Plot
figure
set(gcf,'WindowState','maximized')
idxColumn = 0;
for i = 1:length(cuCoarseList)
    semilogy(1e-14*ones(1,2),1e-14*ones(1,2),'LineWidth',5,'LineStyle','-','Color',colourList(i,:))
    hold on
end
for i = 1:length(cuFineList)
    semilogy(1e-14*ones(1,2),1e-14*ones(1,2),'LineWidth',5,'LineStyle',lineStyleList{i},'Color','k')
end
for i = [2,1]
    semilogy(1e-14,1e-14,'LineWidth',3,'LineStyle','none','Color','k', ...
        'Marker',markers{i},'MarkerFaceColor','k','MarkerSize',20)
end
for i = 1:length(cuCoarseList)
    for j = 1:length(densityList)
        idxColumn = idxColumn+1;
        
        for k = 1:length(cuFineList)
            idxLower = (k-1)*length(rList)+1;
            idxUpper = k*length(rList);
            semilogy(ratio1(idxLower:idxUpper,idxColumn),normalisedVolumeEroded(idxLower:idxUpper,idxColumn),'LineWidth',5,'LineStyle',lineStyleList{k},'Color',colourList(idxColumn,:))
        end
    end
end
idxColumn = 0;
for i = 1:length(cuCoarseList)
    for j = 1:length(densityList)
        idxColumn = idxColumn+1;
        
        idxRow = 0;
        for k = 1:length(cuFineList)
            for l = 1:length(rList)
                idxRow = idxRow+1;
                semilogy(ratio1(idxRow,idxColumn),normalisedVolumeEroded(idxRow,idxColumn),'LineWidth',3,'LineStyle','none','Color',colourList(idxColumn,:), ...
                    'Marker',markerList{idxRow,idxColumn},'MarkerFaceColor',colourList(idxColumn,:),'MarkerSize',20)
            end
        end
    end
end
% plot(xlim,0.15*ones(1,2),'LineWidth',5,'LineStyle','--','Color','k')
% xlim([1.5,7.5])
ylim(yl)
yticks([0.001,0.01,0.1,1,10,100])
semilogy(boundaryX1*ones(1,2),ylim,'LineWidth',5,'LineStyle','--','Color','k')
semilogy(xlim,boundaryY*ones(1,2),'LineWidth',5,'LineStyle','--','Color','k')
legend(legendList,'Location','northwest','FontName','Helvetica','FontSize',20)
set(gca,'FontName','Helvetica','FontSize',36)
xlabel('Ratio of D_1_5_F to d_8_5_B [-]','FontName','Helvetica','FontSize',36)
ylabel(yString,'FontName','Helvetica','FontSize',36)

figure
set(gcf,'WindowState','maximized')
idxColumn = 0;
for i = 1:length(cuCoarseList)
    for j = 1:length(densityList)
        idxColumn = idxColumn+1;
        
        for k = 1:length(cuFineList)
            idxLower = (k-1)*length(rList)+1;
            idxUpper = k*length(rList);
            semilogy(ratio2(idxLower:idxUpper,idxColumn),normalisedVolumeEroded(idxLower:idxUpper,idxColumn),'LineWidth',5,'LineStyle',lineStyleList{k},'Color',colourList(idxColumn,:))
            hold on
        end
    end
end
idxColumn = 0;
for i = 1:length(cuCoarseList)
    for j = 1:length(densityList)
        idxColumn = idxColumn+1;
        
        idxRow = 0;
        for k = 1:length(cuFineList)
            for l = 1:length(rList)
                idxRow = idxRow+1;
                semilogy(ratio2(idxRow,idxColumn),normalisedVolumeEroded(idxRow,idxColumn),'LineWidth',3,'LineStyle','none','Color',colourList(idxColumn,:), ...
                    'Marker',markerList{idxRow,idxColumn},'MarkerFaceColor',colourList(idxColumn,:),'MarkerSize',20)
            end
        end
    end
end
% xlim([1.5,7.5])
ylim(yl)
yticks([0.001,0.01,0.1,1,10,100])
plot(boundaryX1*ones(1,2),ylim,'LineWidth',5,'LineStyle','--','Color','k')
set(gca,'FontName','Helvetica','FontSize',36)
xlabel('Ratio of D_1_5_F to d_9_5_B [-]','FontName','Helvetica','FontSize',36)
ylabel(yString,'FontName','Helvetica','FontSize',36)

figure
set(gcf,'WindowState','maximized')
idxColumn = 0;
for i = 1:length(cuCoarseList)
    for j = 1:length(densityList)
        idxColumn = idxColumn+1;
        
        for k = 1:length(cuFineList)
            idxLower = (k-1)*length(rList)+1;
            idxUpper = k*length(rList);
            semilogy(ratio3(idxLower:idxUpper,idxColumn),normalisedVolumeEroded(idxLower:idxUpper,idxColumn),'LineWidth',5,'LineStyle',lineStyleList{k},'Color',colourList(idxColumn,:))
            hold on
        end
    end
end
idxColumn = 0;
for i = 1:length(cuCoarseList)
    for j = 1:length(densityList)
        idxColumn = idxColumn+1;
        
        idxRow = 0;
        for k = 1:length(cuFineList)
            for l = 1:length(rList)
                idxRow = idxRow+1;
                semilogy(ratio3(idxRow,idxColumn),normalisedVolumeEroded(idxRow,idxColumn),'LineWidth',3,'LineStyle','none','Color',colourList(idxColumn,:), ...
                    'Marker',markerList{idxRow,idxColumn},'MarkerFaceColor',colourList(idxColumn,:),'MarkerSize',20)
            end
        end
    end
end
% xlim([1.5,7.5])
ylim(yl)
yticks([0.001,0.01,0.1,1,10,100])
plot(boundaryX1*ones(1,2),ylim,'LineWidth',5,'LineStyle','--','Color','k')
plot(xlim,boundaryY*ones(1,2),'LineWidth',5,'LineStyle','--','Color','k')
set(gca,'FontName','Helvetica','FontSize',36)
xlabel('Ratio of D_1_5_F to d_8_5_B_(_S_A_) [-]','FontName','Helvetica','FontSize',36)
ylabel(yString,'FontName','Helvetica','FontSize',36)

figure
set(gcf,'WindowState','maximized')
idxColumn = 0;
for i = 1:length(cuCoarseList)
    for j = 1:length(densityList)
        idxColumn = idxColumn+1;
        
        for k = 1:length(cuFineList)
            idxLower = (k-1)*length(rList)+1;
            idxUpper = k*length(rList);
            semilogy(ratio4(idxLower:idxUpper,idxColumn),normalisedVolumeEroded(idxLower:idxUpper,idxColumn),'LineWidth',5,'LineStyle',lineStyleList{k},'Color',colourList(idxColumn,:))
            hold on
        end
    end
end
idxColumn = 0;
for i = 1:length(cuCoarseList)
    for j = 1:length(densityList)
        idxColumn = idxColumn+1;
        
        idxRow = 0;
        for k = 1:length(cuFineList)
            for l = 1:length(rList)
                idxRow = idxRow+1;
                semilogy(ratio4(idxRow,idxColumn),normalisedVolumeEroded(idxRow,idxColumn),'LineWidth',3,'LineStyle','none','Color',colourList(idxColumn,:), ...
                    'Marker',markerList{idxRow,idxColumn},'MarkerFaceColor',colourList(idxColumn,:),'MarkerSize',20)
            end
        end
    end
end
plot(boundaryX2*ones(1,2),ylim,'LineWidth',5,'LineStyle','--','Color','k')
plot(xlim,boundaryY*ones(1,2),'LineWidth',5,'LineStyle','--','Color','k')
% xlim([1.5,7.5])
ylim(yl)
yticks([0.001,0.01,0.1,1,10,100])
set(gca,'FontName','Helvetica','FontSize',36)
xlabel('Ratio of D_3_5_c to d_8_5_B [-]','FontName','Helvetica','FontSize',36)
ylabel(yString,'FontName','Helvetica','FontSize',36)

figure
set(gcf,'WindowState','maximized')
idxColumn = 0;
for i = 1:length(cuCoarseList)
    for j = 1:length(densityList)
        idxColumn = idxColumn+1;
        
        for k = 1:length(cuFineList)
            idxLower = (k-1)*length(rList)+1;
            idxUpper = k*length(rList);
            semilogy(ratio5(idxLower:idxUpper,idxColumn),normalisedVolumeEroded(idxLower:idxUpper,idxColumn),'LineWidth',5,'LineStyle',lineStyleList{k},'Color',colourList(idxColumn,:))
            hold on
        end
    end
end
idxColumn = 0;
for i = 1:length(cuCoarseList)
    for j = 1:length(densityList)
        idxColumn = idxColumn+1;
        
        idxRow = 0;
        for k = 1:length(cuFineList)
            for l = 1:length(rList)
                idxRow = idxRow+1;
                semilogy(ratio5(idxRow,idxColumn),normalisedVolumeEroded(idxRow,idxColumn),'LineWidth',3,'LineStyle','none','Color',colourList(idxColumn,:), ...
                    'Marker',markerList{idxRow,idxColumn},'MarkerFaceColor',colourList(idxColumn,:),'MarkerSize',20)
            end
        end
    end
end
plot(boundaryX2*ones(1,2),ylim,'LineWidth',5,'LineStyle','--','Color','k')
% xlim([1.5,7.5])
ylim(yl)
yticks([0.001,0.01,0.1,1,10,100])
set(gca,'FontName','Helvetica','FontSize',36)
xlabel('Ratio of D_3_5_c to d_9_5_B [-]','FontName','Helvetica','FontSize',36)
ylabel(yString,'FontName','Helvetica','FontSize',36)

figure
set(gcf,'WindowState','maximized')
idxColumn = 0;
for i = 1:length(cuCoarseList)
    for j = 1:length(densityList)
        idxColumn = idxColumn+1;
        
        for k = 1:length(cuFineList)
            idxLower = (k-1)*length(rList)+1;
            idxUpper = k*length(rList);
            semilogy(ratio6(idxLower:idxUpper,idxColumn),normalisedVolumeEroded(idxLower:idxUpper,idxColumn),'LineWidth',5,'LineStyle',lineStyleList{k},'Color',colourList(idxColumn,:))
            hold on
        end
    end
end
idxColumn = 0;
for i = 1:length(cuCoarseList)
    for j = 1:length(densityList)
        idxColumn = idxColumn+1;
        
        idxRow = 0;
        for k = 1:length(cuFineList)
            for l = 1:length(rList)
                idxRow = idxRow+1;
                semilogy(ratio6(idxRow,idxColumn),normalisedVolumeEroded(idxRow,idxColumn),'LineWidth',3,'LineStyle','none','Color',colourList(idxColumn,:), ...
                    'Marker',markerList{idxRow,idxColumn},'MarkerFaceColor',colourList(idxColumn,:),'MarkerSize',20)
            end
        end
    end
end
plot(boundaryX2*ones(1,2),ylim,'LineWidth',5,'LineStyle','--','Color','k')
plot(xlim,boundaryY*ones(1,2),'LineWidth',5,'LineStyle','--','Color','k')
% xlim([1.5,7.5])
ylim(yl)
yticks([0.001,0.01,0.1,1,10,100])
set(gca,'FontName','Helvetica','FontSize',36)
xlabel('Ratio of D_3_5_c to d_8_5_B_(_S_A_) [-]','FontName','Helvetica','FontSize',36)
ylabel(yString,'FontName','Helvetica','FontSize',36)

%% Overall plot

close all

% Legend list
legendList = {'Base C_u = 1.5','Base C_u = 3.0','Base C_u = 4.5', ...
    'Seal','No seal'};

% Options
yl = [1e-3,1e2];
boundaryY = 0.1;
boundaryX1 = 4;
boundaryX2 = 1;
yString = 'Normalised total volume eroded';

% Plot properties
lineStyleList = {'-','--',':'};
colourList = Constants.colourList;

% Plot
figure
set(gcf,'WindowState','maximized')
for i = 1:length(cuFineList)
    semilogy(1e-14,1e-14,'LineWidth',3,'LineStyle','none','Color',colourList(i,:), ...
        'Marker',markers{1},'MarkerFaceColor',colourList(i,:),'MarkerSize',20)
    hold on
end
for i = [2,1]
    semilogy(1e-14,1e-14,'LineWidth',3,'LineStyle','none','Color','k', ...
        'Marker',markers{i},'MarkerFaceColor','k','MarkerSize',20)
end
idxColumn = 0;
for i = 1:length(cuCoarseList)
    for j = 1:length(densityList)
        idxColumn = idxColumn+1;
        
        idxRow = 0;
        for k = 1:length(cuFineList)
            for l = 1:length(rList)
                idxRow = idxRow+1;
                semilogy(ratio5(idxRow,idxColumn),normalisedVolumeEroded(idxRow,idxColumn),'LineWidth',3,'LineStyle','none','Color',colourList(k,:), ...
                    'Marker',markerList{idxRow,idxColumn},'MarkerFaceColor',colourList(k,:),'MarkerSize',20)
            end
        end
    end
end
plot(boundaryX2*ones(1,2),ylim,'LineWidth',5,'LineStyle','--','Color','k')
% xlim([1.5,7.5])
legend(legendList,'Location','northwest','FontName','Helvetica','FontSize',20)
ylim(yl)
yticks([0.001,0.01,0.1,1,10,100])
set(gca,'FontName','Helvetica','FontSize',36)
xlabel('Ratio of D_3_5_c to d_9_5_B [-]','FontName','Helvetica','FontSize',36)
ylabel(yString,'FontName','Helvetica','FontSize',36)

figure
set(gcf,'WindowState','maximized')
idxColumn = 0;
for i = 1:length(cuCoarseList)
    for j = 1:length(densityList)
        idxColumn = idxColumn+1;
        
        idxRow = 0;
        for k = 1:length(cuFineList)
            for l = 1:length(rList)
                idxRow = idxRow+1;
                semilogy(ratio6(idxRow,idxColumn),normalisedVolumeEroded(idxRow,idxColumn),'LineWidth',3,'LineStyle','none','Color',colourList(k,:), ...
                    'Marker',markerList{idxRow,idxColumn},'MarkerFaceColor',colourList(k,:),'MarkerSize',20)
                hold on
            end
        end
    end
end
plot(boundaryX2*ones(1,2),ylim,'LineWidth',5,'LineStyle','--','Color','k')
plot(xlim,boundaryY*ones(1,2),'LineWidth',5,'LineStyle','--','Color','k')
% xlim([1.5,7.5])
ylim(yl)
yticks([0.001,0.01,0.1,1,10,100])
set(gca,'FontName','Helvetica','FontSize',36)
xlabel('Ratio of D_3_5_c to d_8_5_B_(_S_A_) [-]','FontName','Helvetica','FontSize',36)
ylabel(yString,'FontName','Helvetica','FontSize',36)