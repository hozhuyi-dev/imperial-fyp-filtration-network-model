% The script add the boundary of the filter

clear;
clc;
close all;

% Data available
% cuCoarseList = Constants.cuCoarse;
densityList = Constants.density;
% sizeList = Constants.size;
cuCoarseList = 3;
% densityList = {'dense'};
sizeList = [1,1,1; 1,1,4];
triangulationMethod = 'Del';

for i = 1:length(cuCoarseList)
    for j = 1:length(densityList)
        for k = 1:size(sizeList,1)
            % Create filter
            filter = Filter(1,cuCoarseList(i),densityList{j},sizeList(k,:),triangulationMethod,'actual','actual',[],'zpos','undirected');
            
            % Load dtmCleaned
            filter = filter.loadState(FileLocation.filterDTMCleaned(filter), ...
                'variables', {'nodeNum', 'nodeCoord', 'nodeBoundary', 'nodeVolume', 'nodeSurfaceArea', ...
                'edgeNum', 'edgeArea', 'edgeCoord', 'edgeDiameter', 'edgeLength', 'edgeList', 'edgeNormal', 'edgeParticles', 'edgeReduced', ...
                'particleNum', 'particleCoord', 'particleRadius'});
            
            % Read the boundary coordinates
            s = fileread(FileLocation.filterBoundarySize(filter));
            idxStart = regexp(s,']');
            idxEnd = regexp(s,'N');
            idxEnd = idxEnd(end);
            s = s(idxStart+1:idxEnd-1);
            s = regexp(s,'\s','split');
            boundary = zeros(6,1);
            i2 = 0;
            for i1 = 1:length(s)
                if str2double(s{i1})>0
                    i2 = i2+1;
                    boundary(i2) = str2double(s{i1});
                end
            end
            
            % Load particle
            particleData = load(FileLocation.filterDTMRaw(filter,'particle'));
            particleData = particleData.(sprintf('%s',FileLocation.filterDTMRawVariable(filter,'particle')));
            particleRadius = particleData(:,size(particleData,2));
            
            % Normalize
            boundary = boundary*1/(2*min(particleRadius));
            
            % All nodes must be within boundary
            assert(min(filter.nodeCoord(:,1))>boundary(1));
            assert(max(filter.nodeCoord(:,1))<boundary(2));
            assert(min(filter.nodeCoord(:,2))>boundary(3));
            assert(max(filter.nodeCoord(:,2))<boundary(4));
            assert(min(filter.nodeCoord(:,3))>boundary(5));
            assert(max(filter.nodeCoord(:,3))<boundary(6));
            
            % Store
            filter.boundary = boundary;
            
            % Save dtmCleaned
            filter.saveState(FileLocation.filterDTMCleaned(filter), ...
                'variables', {'nodeNum', 'nodeCoord', 'nodeBoundary', 'nodeVolume', 'nodeSurfaceArea', ...
                'edgeNum', 'edgeArea', 'edgeCoord', 'edgeDiameter', 'edgeLength', 'edgeList', 'edgeNormal', 'edgeParticles', 'edgeReduced', ...
                'particleNum', 'particleCoord', 'particleRadius', ...
                'boundary'});
        end
    end
end