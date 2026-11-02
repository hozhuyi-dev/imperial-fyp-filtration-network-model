% The script analyses how the network model compares to Terzaghi's rule

clear
clc
close all

% Add src to path
path(path, ['..', filesep, 'src'])

% Simulations to consider
transportLogicType = 'flowbiased'; 
networkType = 'actual';
assignmentMethod = 'actual';
cuCoarseList = Constants.cuCoarse;
densityList = Constants.density;
siz = [1,1,1];
triangulationMethod = 'Del';
cuFineList = 2;

for i = 1:length(cuCoarseList)
    for j = 1:length(densityList)
        for k = 1:length(cuFineList)
            % Load simulator
            simulator = Simulator.loadResults(FileLocation.simulatorPolydisperse(transportLogicType, ...
                networkType, assignmentMethod, cuCoarseList(i), densityList{j}, siz, triangulationMethod, cuFineList(k)));
            
            % Determine Terzaghi ratio
            r.(sprintf('r%i%i%i',i,j,k)) = simulator.terzaghiRatio;
            
            % Determine percentage eroded by volume
            p.(sprintf('p%i%i%i',i,j,k)) = simulator.percentageErodedVol;
            
            fprintf('%i%i%i',i,j,k)
        end
        fprintf('\n')
    end
end

%% Plot

legendList = cell(length(cuCoarseList)*length(densityList),1);
k = 0;
for i = 1:length(cuCoarseList)
    for j = 1:length(densityList)
        k = k+1;
        legendList{k} = sprintf('cuCoarse %.1f, %s',cuCoarseList(i),densityList{j});
    end
end
markerList = {'x','o','^'};
lineStyleList = {'-','--',':'};
colourList = Constants.colourList;

figure
set(gcf,'WindowState','maximized')
for i = 1:length(cuCoarseList)
    for j = 1:length(densityList)
        for k = 1:length(cuFineList)
            plot(r.(sprintf('r%i%i%i',i,j,k)),p.(sprintf('p%i%i%i',i,j,k)), ...
                'Color',colourList(i,:),'LineWidth',3,'LineStyle',lineStyleList{j},'Marker','none','MarkerSize',18)
            hold on
        end
    end
end
% plot((d.d1),((p.p1/100).^2)*100,'LineWidth',3,'LineStyle',lineStyleList{(i>2)+1})
legend(legendList,'Location','southeast','FontName','Helvetica','FontSize',18)
% xlim([min(d.d1),max(d.d1)])
set(gca, 'FontName','Helvetica','FontSize',18)
xlabel('Normalised diameter of monodisperse base particles [-]','FontName','Helvetica','FontSize',18)
ylabel('Percentage of base particles eroded','FontName','Helvetica','FontSize',18)