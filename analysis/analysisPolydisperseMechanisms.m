% The script analyses how the behaviours observed in monodisperse
% simulations are manifested in combinatory 

clear
clc
close all

% Add src to path
path(path, ['..', filesep, 'src'])

%% Largest particle eroded

% Percolation thresholds based on directed networks
pDirected = [0.5809,0.5985,0.5879; ...
    0.5864,0.5380,0.5844; ...
    0.5714,0.5787,0.5875; ...
    0.5780,0.5865,0.5672];
% pDirected = 0.5809;

% Simulations to consider
transportLogicType = 'flowbiased';
networkType = 'actual';
assignmentMethod = 'actual';
cuCoarseList = Constants.cuCoarse;
densityList = Constants.density;
siz = [1,1,1];
triangulationMethod = 'Del';
cuFineList = [1.5,3,4.5];
rList = 0.5:0.5:7;

% Plot properties
markerList = {'o','s','^','d'};
colourList = Constants.colourList;
lineStyleList = {'-','--',':'};

%%

% Initialise arrays
arraySize = length(cuCoarseList)*length(densityList)*length(cuFineList)*length(rList);
dCorresponding = zeros(length(cuCoarseList),length(densityList));
d = zeros(length(cuCoarseList),length(densityList),length(cuFineList),length(rList));
d95 = d;
dMax = d;
dMin = d;
D35c = dCorresponding;

fprintf('Loading largest particle eroded ')
idx = 0;
for i = 1:length(cuCoarseList)
    for j = 1:length(densityList)
        for k = 1:length(cuFineList)
            for l = 1:length(rList)
                idx = idx+1;
                
                % Load simulator
                simulator = Simulator.loadResults(FileLocation.simulatorPolydisperseRepeat(transportLogicType, ...
                    networkType, assignmentMethod, cuCoarseList(i), densityList{j}, siz, triangulationMethod, ...
                    cuFineList(k), rList(l), 'continuous'));
                
                % Determine the largest particle eroded
                d(i,j,k,l) = simulator.largestParticleEroded;
                d95(i,j,k,l) = simulator.largestParticleEroded95;
                
                % Largest simulated base particles
                dMax(i,j,k,l) = simulator.base.base1.maxFine;
                dMin(i,j,k,l) = simulator.base.base1.minFine;
                
                if idx~=1
                    fprintf(repmat('\b', 1, 7));
                end
                fprintf('%6.2f%%', idx/arraySize*100);
            end
        end
        
        % Use of 35%
        D35c(i,j) = simulator.filter.percentileEdgeDiameter(35);
        
        % Corresponding largest particle eroded
        dCorresponding(i,j) = simulator.filter.percentileEdgeDiameter((1-pDirected(i,j))*100);
    end
end
fprintf('\n')

% Finalise
% d(d==0) = NaN;
% d95(d95==0) = NaN;

%% Save

save(FileLocation.analysisPolydisperseMechanisms, 'd', 'd95', 'dMax', 'dMin', 'D35c', 'dCorresponding');

%% Load

data = load(FileLocation.analysisPolydisperseMechanisms);
d = data.d;
d95 = data.d95;
dMax = data.dMax;
dMin = data.dMin;
D35c = data.D35c;
dCorresponding = data.dCorresponding;

%% Plot

% legendList
legendList = cell(length(cuFineList)+3,1);
for i = 1:length(cuFineList)
    legendList{i} = sprintf('Base C_u = %.1f',cuFineList(i));
end
legendList{4} = 'Smallest and largest base simulated';
legendList{5} = sprintf('Predicted largest particle eroded = %.3f',dCorresponding(1,3));
legendList{6} = sprintf('D_{35c}^{*} = %.3f',D35c(1,3));

for i = 1
    figure
    set(gcf,'WindowState','maximized')
    for j = 3
        for k = 1:length(cuFineList)
            dCurrent = d(i,j,k,:);
            dCurrent = dCurrent(:);
            
            plot(rList,dCurrent,'LineWidth',5,'Color',colourList(k,:),'LineStyle',lineStyleList{1})
            hold on
        end
    end
    plot([0,0],[0,0],'LineWidth',2,'Color','k','LineStyle','--');
    for j = 3
        plot(xlim,dCorresponding(i,j)*ones(1,2),'LineWidth',3,'Color','k','LineStyle',':');
        plot(xlim,D35c(i,j)*ones(1,2),'LineWidth',3,'Color','k','LineStyle','-.');
    end
    for j = 3
        for k = 1:length(cuFineList)
            % Largest particle simulated
            dCurrent = dMax(i,j,k,:);
            dCurrent = dCurrent(:);
            plot(rList,dCurrent,'LineWidth',3,'Color',colourList(k,:),'LineStyle','--');
            
            % Smallest particle simulated
            dCurrent = dMin(i,j,k,:);
            dCurrent = dCurrent(:);
            plot(rList,dCurrent,'LineWidth',3,'Color',colourList(k,:),'LineStyle','--');
        end
    end
    xlim([0.5,7])
    ylim([0,0.6])
    xticks(0.5:0.5:7)
    legend(legendList,'Location','northeast','FontName','Helvetica','FontSize',26)
    set(gca,'FontName','Helvetica','FontSize',36)
    xlabel('Ratio of D_1_5_F to d_8_5_B','FontName','Helvetica','FontSize',36)
    ylabel('Normalised largest erodible particle','FontName','Helvetica','FontSize',36)
end

%% Plot correlation between largest particle eroded 

% Calculate correlation coefficient
data = zeros(length(cuCoarseList)*length(densityList)*length(cuFineList),1);
idx = 0;
for i = 1:length(cuCoarseList)
    for j = 1:length(densityList)
        % Current D35c
        D35cCurrent = D35c(i,j);
        
        % For each cuFine
        for k = 1:length(cuFineList)
            idx = idx+1;
            % Current vector of d
            dCurrent = d(i,j,k,:);
            dCurrent = dCurrent(:);
            dCurrent = max(dCurrent);
            
            data(idx,1) = D35cCurrent;
            data(idx,2) = dCurrent;
        end
    end
end
R = corrcoef(data(:,1),data(:,2));
P = R(1,2);

legendList = {'Base C_{u} = 1.5','Base C_{u} = 3.0','Base C_{u} = 4.5', ...
    'Filter C_{u} = 1.2','Filter C_{u} = 1.5','Filter C_{u} = 2.0','Filter C_{u} = 3.0', ...
    'x = y'};

fig = figure;
set(gcf,'WindowState','maximized')
for i = 1:length(cuFineList)
    plot(0,0,'LineStyle','none','Marker','o','MarkerSize',20,'Color',colourList(i,:),'MarkerFaceColor',colourList(i,:))
    hold on
end
for i = 1:length(cuCoarseList)
    plot(0,0,'LineStyle','none','Marker',markerList{i},'MarkerSize',20,'Color','k','MarkerFaceColor','k')
end
plot([0,0.5],[0,0.5],'LineStyle','--','Color','k','LineWidth',3)
for i = 1:length(cuCoarseList)
    for j = 1:length(densityList)
        % Current D35c
        D35cCurrent = D35c(i,j);
        
        % For each cuFine
        for k = 1:length(cuFineList)
            % Current vector of d
            dCurrent = d(i,j,k,:);
            dCurrent = dCurrent(:);
            dCurrent = max(dCurrent);
            
            % Plot
            plot(D35cCurrent,dCurrent,'LineStyle','none','Marker',markerList{i},'MarkerSize',20,'Color',colourList(k,:),'MarkerFaceColor',colourList(k,:))
        end
    end
end
text(0.495,0.22,sprintf('Correlation coefficient = %.4f',P),'FontName','Helvetica','FontSize',26,'HorizontalAlignment','right')
xlim([0.2,0.5])
ylim([0.2,0.5])
legend(legendList,'Location','northwest','FontName','Helvetica','FontSize',26)
set(gca,'FontName','Helvetica','FontSize',36)
xlabel('Controlling constriction size D_{35c}^{*}','FontName','Helvetica','FontSize',36)
ylabel('Normalised largest erodible particle','FontName','Helvetica','FontSize',36)

%% 

% Simulation to consider
transportLogicType = 'flowbiased';
networkType = 'actual';
assignmentMethod = 'actual';
cuCoarse = 1.2;
density = 'dense';
siz = [1,1,1];
triangulationMethod = 'Del';
cuFine = 4.5;
r = 4;

% Load simulator
simulator = Simulator.loadResults(FileLocation.simulatorPolydisperseRepeat(transportLogicType, ...
    networkType, assignmentMethod, cuCoarse, density, siz, triangulationMethod, ...
    cuFine, r, 'continuous'));

%% PSDs of particles retained with time

simulator.statistics.statistics1.plotTimeRetainedPSD(simulator.base.base1.psdFineBucket, ...
    'idxPlot',[2,5,50,100,277],'csdNum',simulator.filter.csdNum,'psdVol',simulator.base.base1.psdFineVol, ...
    'xl',[0.02,1]);

%% PSDs of particles eroded with time

simulator.statistics.statistics1.plotTimeErodedPSD(simulator.base.base1.psdFineBucket, ...
    'idxPlot',[47,50,60,277],'csdNum',simulator.filter.csdNum,'psdVol',simulator.base.base1.psdFineVol, ...
    'xl',[0.02,1]);

%% Volume of particles retained in each filter layer

simulator.statistics.statistics1.plotLayerRetainedVol(simulator.base.base1.psdFineBucket, ...
    'idxRange',[1,68,86,101]);