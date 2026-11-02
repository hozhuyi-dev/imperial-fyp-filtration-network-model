% The script determines the number of nodes and edges, as well as
% connectivity

clear;
clc;
close all;

path(path, ['..', filesep, 'src']);

cuCoarse = 3;
triangulationMethod = 'Del';
densityList = Constants.density;
sizeList = [1,1,1; 1,1,4];

nodeNum = zeros(length(densityList)*size(sizeList,1),1);
edgeNum = nodeNum;
nodeConnectivity = nodeNum;

k = 0;
for i = 1:length(densityList)
    for j = 1:size(sizeList,1)
        k = k+1;
        
        simulator = Simulator();
        simulator = simulator.createFilter(cuCoarse,triangulationMethod,'density',densityList{i},'size',sizeList(j,:));
        
        nodeNum(k) = simulator.filter.nodeNum;
        edgeNum(k) = simulator.filter.edgeNum;
        nodeConnectivity(k) = mean(simulator.filter.nodeConnectivity);
    end
end