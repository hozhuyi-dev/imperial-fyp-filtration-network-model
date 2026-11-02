% Update all networks with nodeEdgeDirection

clear;
clc;
close all;

% Add src to path
path(path, ['..', filesep, 'src']);

% Data to consider
cuCoarse = 3;
densityList = Constants.density;
sizeList = Constants.size;
sizeList = sizeList(2:4,:);
% sizeList = [1,1,1];
triangulationMethod = 'Del';
flowDirectionList = Constants.flowDirection;

% For each characteristic
for i = 1:length(densityList)
    for j = 1:size(sizeList,1)
        % For each flow direction
        for k = 1:length(flowDirectionList)
            % Determine the indices
            switch flowDirectionList{k}(1)
                case 'x'
                    idxCoord = 1;
                case 'y'
                    idxCoord = 2;
                case 'z'
                    idxCoord = 3;
            end
            idxOdd = idxCoord*2-1;
            switch flowDirectionList{k}(2:4)
                case 'pos'
                    idxBoundary = [idxOdd, idxOdd+1];
                    coef = 1;
                case 'neg'
                    idxBoundary = [idxOdd+1, idxOdd];
                    coef = -1;
            end
            
            % Load filter
            simulator = Simulator();
            simulator = simulator.createFilter(cuCoarse,'density',densityList{i},'size',sizeList(j,:), ...
                'triangulationMethod',triangulationMethod,'flowDirection',flowDirectionList{k});
            
            if isempty(simulator.filter.nodeEdgeToPos)
                % Determine nodeEdgeToPos
                simulator.filter.nodeEdgeToPos = zeros(size(simulator.filter.nodeEdgeTo));
                for i1 = 1:simulator.filter.nodeNum
                    idxNodes = simulator.filter.nodeNodeTo(i1,1:simulator.filter.nodeConnectivity(i1));
                    for i2 = 1:length(idxNodes)
                        if (coef*(simulator.filter.nodeCoord(idxNodes(i2),idxCoord)-simulator.filter.nodeCoord(i1,idxCoord)))>0
                            simulator.filter.nodeEdgeToPos(i1,i2) = 1;
                        end
                    end
                end
                
                % Save nodeEdgeDirection in network data file
                simulator.filter.saveState(FileLocation.filterNetwork(simulator.filter), 'variables', {'nodeInlet', 'nodeLattice', 'nodeOutlet', ...
                    'nodeConnectivity', 'nodeEdgeTo', 'nodeNodeTo', 'nodeEdgeToPos', 'nodeEdgeDirection'});
            end
        end
    end
end