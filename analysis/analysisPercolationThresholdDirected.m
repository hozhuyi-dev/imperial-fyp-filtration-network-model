% The script calculates the percolation threshold and corresponding
% d(largest), mean and standard deviation

clear
clc
close all

% Add src to path
path(path, ['..', filesep, 'src'])

% Networks to consider
cuCoarseList = Constants.cuCoarse;
densityList = Constants.density;
sizeList = Constants.size(2:end,:);

% Initialise array(s)
p = zeros(length(cuCoarseList)*length(densityList),size(sizeList,1));
distance = p;
shortestPath = p;
d = p;
r = p;

% Previously obtained percolation threshold
minP = [0.5809,0.6103,0.6103,0.5324,0.5495,0.5497; ...
    0.5985,0.6166,0.6168,0.5228,0.5677,0.5678; ...
    0.5879,0.5968,0.6041,0.5754,0.5953,0.5951; ...
    0.5864,0.6148,0.6232,0.5602,0.5710,0.5708; ...
    0.5380,0.6111,0.6157,0.5392,0.5581,0.6078; ...
    0.5844,0.5998,0.6067,0.5685,0.5923,0.5921; ...
    0.5714,0.5773,0.5900,0.5615,0.5753,0.5825; ...
    0.5787,0.5802,0.5804,0.5625,0.5815,0.5819; ...
    0.5875,0.6064,0.6067,0.5465,0.6037,0.6064; ...
    0.5780,0.6054,0.6051,0,0,0; ...
    0.5865,0.5898,0.5902,0,0,0; ...
    0.5672,0,0,0,0,0];

i = 0;
for i1 = 1:length(cuCoarseList)
    for i2 = 1:length(densityList)
        i = i+1;
        if i1==4
            if i2==3
                endWith = 1;
            else
                endWith = 3;
            end
        else
            endWith = size(sizeList,1);
        end
        for j = 1:endWith
            % Create filter
            simulator = Simulator;
            simulator = simulator.createFilter(cuCoarseList(i1),'density',densityList{i2},'size',sizeList(j,:));
            
            % Determine percolation threshold
            [p(i,j),distance(i,j)] = simulator.percolationThresholdAssignedDirected('minP',minP(i,j));
            
            % Determine the shortest path
            [~,shortestPath(i,j)] = simulator.filter.shortestPath;
            
            % Determine the path tortuousity
            r(i,j) = distance(i,j)/shortestPath(i,j);
            
            % Corresponding d
            d(i,j) = simulator.filter.percentileEdgeDiameter(100*(1-p(i,j)));
            
            fprintf('%i%i\n',i,j)
        end
    end
end

%% 

% For report
% i = 0;
% for i1 = 1:length(cuCoarseList)
%     for i2 = 1:length(densityList)
%         i = i+1;
%         fprintf('(%.3f) & (%.3f) & (%.3f) & (%.3f) & (%.3f) & (%.3f)\n',d(i,:))
%     end
% end

% Plot percolation threshold against length to width ratio

% Plot path tortuousity against length to width ratio