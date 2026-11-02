classdef Constants
    % The class stores constants and lists that are frequently used
    
    properties (Constant)
        % Data available
        cuCoarse = [1.2,1.5,2,3];
        density = {'loose','medium','dense'};
        availableSize = [0,0,0; ...
            1,1,4; ...
            2,2,4];
        size = [0,0,0; ...
            1,1,1; ...
            1,1,2; ...
            1,1,3; ...
            2,2,1; ...
            2,2,2; ...
            2,2,3];
        triangulationMethod = {'Del','Mod'};
        networkType = {'square', 'slantedsquare', 'cubic', 'slantedcubic', 'bcc', 'fcc', ...
            'squarenp', 'slantedsquarenp', 'cubicnp', 'slantedcubicnp', 'bccnp', 'fccnp', 'actual'};
        assignmentMethod = {'actual', 'random'};
        flowDirection = {'xpos', 'xneg', 'ypos', 'yneg', 'zpos', 'zneg'};
        directedOrUndirected = {'directed', 'undirected'};
        
        % Options for Simulator class
        loadMethod = {'read','calculate'};
        introductionMethod = {'refillInlets','fixedParticleNum','fixedParticleNum2'};
        transportLogicType = {'flowbiased', 'flowbiasedpos', 'flowbiasednode', 'flowbiasedextreme', 'noclogging', 'noflow', 'noflowrandom'};
        storeOption = {'normal','storePath','more','storeEdgePath','storeNodeClogging'};
        
        % Constants
        dynamicViscosity = 1.002e-3;
        shapeFactor = 0.5;
        
        % Plot properties
        colourList = [0, 0.4470, 0.7410; ...
            0.8500, 0.3250, 0.0980; ...
            0.9290, 0.6940, 0.1250; ...
            0.4940, 0.1840, 0.5560; ...
            0.4660, 0.6740, 0.1880; ...
            0.3010, 0.7450, 0.9330; ...
            0.6350, 0.0780, 0.1840; ...
            0.75, 0.75, 0; ...
            0, 0.5, 0; ...
            0.25, 0.25, 0.25];
        fontSize = 36;
        lineWidth = 5;
    end
    
    methods (Static)
        % Output the original size of the filter to be cut from
        function [size] = cutFrom(currentSize)
            if any(all(currentSize==Constants.availableSize,2))
                size = currentSize;
            elseif any(all(currentSize==[1,1,1; 1,1,2; 1,1,3],2))
                size = [1,1,4];
            elseif any(all(currentSize==[2,2,1; 2,2,2; 2,2,3],2))
                size = [2,2,4];
            end
        end
        
        % Determine the size of regular networks such that the number of
        % nodes is comparable (compared to the actual network wDTM)
        function [nxyz] = comparableNxyz(nodeNum, networkType, size)
            % Initialise nxyz
            nxyz = zeros(1,3);
            
            % Determine the corresponding nxyz
            assert(size(1)==size(2));
            ratio = size(3)/size(1);
            x = ceil((nodeNum/ratio)^(1/3));
            if ismember(networkType, {'square','squarenp','slantedsquare','slantedsquarenp'})
                nxyz(1) = x^2;
                nxyz(2) = 1;
                nxyz(3) = ceil(x*ratio);
            elseif ismember(networkType, {'cubic','cubicnp','bcc','bccnp'})
                nxyz(1:2) = x;
                nxyz(3) = ceil(x*ratio);
            elseif ismember(networkType, {'slantedcubic','slantedcubicnp'})
                nxyz(1:2) = ceil(x/2);
                nxyz(3) = ceil(x*ratio);
            elseif ismember(networkType, {'fcc','fccnp'})
                nxyz(1:2) = ceil(((x^2)/2)^(1/2));
                nxyz(3) = ceil(x*ratio);
            end
        end
        
        % The void ratio of the filter as given by Shire (2014)
        function [e] = voidRatio(cuCoarseList,densityList)
            % Gives the void ratios of the samples by Shire (2014)
            eList = [0.714,0.651,0.558; ...
                0.658,0.608,0.531; ...
                0.555,0.523,0.467; ...
                0.455,0.426,0.382]; 
            
            e = zeros(length(cuCoarseList),length(densityList));
            for i = 1:length(cuCoarseList)
                for j = 1:length(densityList)
                    [~,x] = ismember(cuCoarseList(i),Constants.cuCoarse);
                    [~,y] = ismember(densityList{j},Constants.density);
                    
                    e(i,j) = eList(x,y);
                end
            end
        end
    end
end