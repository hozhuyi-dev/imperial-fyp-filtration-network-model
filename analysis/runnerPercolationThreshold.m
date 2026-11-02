% The script evaluates how network size affects percolation threshold

clear;
clc;
close all;

% Add src to path
path(path, ['..', filesep, 'src']);

%% For cubic network

% Network size(s) to consider
networkType = 'cubicnp';
isPeriodic = 0;
networkSizeList = [30,30,30; ...
    30,30,60; ...
    30,30,90; ...
    30,30,120];
minPList = [0.22;0.25;0.252;0.255];
m = size(networkSizeList,1);

% Determine p n times
n = 150;

% Initialise array(s)
pData = NaN(m,n);

for i = 1:m
    fprintf('Network %i\n', i);
    % Generate regular network
    filter = Filter([],[],[],[],[],networkType,[],networkSizeList(i,:),[],'undirected');
    filter = filter.generateNetworkTopologyRegular();
    
    for j = 1:n
        fprintf('Run %i\n', j);
        % Determine percolation threshold
        pData(i,j) = filter.percolationThreshold('minP',minPList(i));
    end
    
    % Save file
    p = pData(i,:);
    save(FileLocation.percolationThreshold(networkType,'random',[],[],networkSizeList(i,:),isPeriodic), 'p');
end

%% For actual network

% Data to consider
networkType = 'actual';
isPeriodic = 0; % No option for periodic
cuCoarseList = 1.2;
triangulationMethodList = Constants.triangulationMethod;
minPList = [0.5,0.25];
m = length(cuCoarseList)*length(triangulationMethodList);
flowDirectionList = {'xpos','ypos','zpos'};
assignmentMethod = 'actual';

% Determine p n times (n/3 times for each flow direction)
n = 150;

% Initialise array(s)
pData = NaN(m,n);
idxX = 0;

for i = 1:length(cuCoarseList)
    for j = 1:length(triangulationMethodList)
        % Current iteration
        idxX = idxX+1;
        fprintf('Network %i\n', idxX);
        idxY = 0;
        
        for k = 1:length(flowDirectionList)
            % Generate regular network
            filter = Filter(1,cuCoarseList(i),'dense',[0,0,0],triangulationMethodList{j},'actual','actual',[],flowDirectionList{k},'undirected');
            filter = filter.loadDTM('read');
            filter = filter.generateNetwork();
            
            for l = 1:n/3
                idxY = idxY+1;
                fprintf('Run %i\n', idxY);
                
                % Determine percolation threshold
                if strcmp(assignmentMethod, 'actual')
                    pData(idxX,idxY) = filter.percolationThresholdAssigned('minP',minPList(idxX));
                else
                    pData(idxX,idxY) = filter.percolationThreshold('minP',minPList(idxX));
                end
            end
        end
        
        % Save data
        p = pData(i,:);
        save(FileLocation.percolationThreshold(networkType,assignmentMethod,cuCoarseList(i),triangulationMethodList{j},[],[]), 'p');
    end
end