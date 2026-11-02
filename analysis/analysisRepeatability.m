% The script plots p vs d to show repeatability

clear;
clc;
close all;

% Add src to path
path(path, ['..', filesep, 'src']);

% Networks to consider
transportLogicType = 'flowbiased';
networkType = 'actual';
assignmentMethod = 'actual';
cuCoarse = 1.2;
density = 'dense';
sizee = [1,1,1];
triangulationMethod = 'Del';
repeatNum = 10;

% Load
for i = 1:repeatNum
    % Load Simulator object
    simulator = Simulator.loadResults(FileLocation.simulatorMonodisperse(transportLogicType,networkType, ...
        assignmentMethod,cuCoarse,density,sizee,triangulationMethod,sprintf('%i',i)));
    
    % Determine the diameter of particles
    d.(sprintf('d%i',i)) = simulator.monodisperseBaseDiameter();
    
    % Determine the percentage of particles eroded
    p.(sprintf('p%i',i)) = simulator.percentageErodedNum('continuous',0);

    fprintf('Loaded %6.2f%%\n', i/repeatNum*100);
end
%%
% Determine the maximum, minimum and mean
d = d.d1;
pMean = zeros(size(d));
pMin = pMean;
pMax = pMean;
for i = 1:length(d)
    % Initialise pData
    pData = zeros(repeatNum,1);
    
    for j = 1:repeatNum
        pData(j) = p.(sprintf('p%i',j))(i);
    end
    
    % Mean 
    pMean(i) = mean(pData);
    
    % Min
    pMin(i) = min(pData);
    
    % Max
    pMax(i) = max(pData);
end

% Figure properties
legendList = {'Mean percentage','Minimum percentage','Maximum percentage'};

% Plot
figure;
set(gcf,'WindowState','maximized');
plot(d,pMean,'LineWidth',4);
hold on;
plot(d,pMin,'LineWidth',4,'LineStyle',':');
plot(d,pMax,'LineWidth',4,'LineStyle',':');
legend(legendList,'Location','northeast','FontName','Helvetica','FontSize',18);
xlim([0.15,0.35]);
set(gca, 'FontName','Helvetica','FontSize',18);
xlabel('Normalised diameter of monodisperse base particles [-]','FontName','Helvetica','FontSize',18);
ylabel('Percentage of base particles eroded','FontName','Helvetica','FontSize',18);