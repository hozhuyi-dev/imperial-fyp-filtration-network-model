classdef Statistics
    % The output of each simulation are stored as this object. The class
    % also includes methods for post-processing
    
    properties
        % Network state
        % The matrix nodeRetained stores:
        % 1: Total number of particles retained in the node
        % 2: Total surface area of particles retained in the node
        % 3: Total volume of particles retained in the node
        % 4: Total number of particles that clog through Mechanism 1
        % 5: Total number of particles that clog through Mechanism 3
        % 6: If the node is clogged (for 'flowbiasedextreme')
        nodeRetained
        % Only for 'storeOption' = 'more' (LARGE DATA SET, NOT RECOMMENDED)
        % The struct timeNodeRetained stores (for each time step): 
        % 1: Number of particles retained in the node
        % 2: Number of open entrances
        % 3: Number of open exits
        % 4: Number of clogged entrances
        % 5: Number of clogged exits
        timeNodeRetained
        % Stores the time step at which a node becomes full, and how each
        % exit is clogged (either through particle straining or node
        % clogging)
        % 1: Time step at which it is clogged
        % 2: Number of exits
        % 3: Number of exits clogged due to particle straining
        % 4: Number of exits clogged due to clogging of adjacent nodes
        nodeClogging
        
        % The matrix edgeFlowing stores:
        % 1: Whether the edge is clogged
        % 2: The mechanism through which it is clogged (1:straining, 3:nospace)
        % 3: The node from which it is clogged
        % 4: The diameter of the particle
        % 5: Time step when it was clogged
        edgeFlowing
        % The array edgePath stores the number of particles passing through
        % each edge
        edgePath % All
        edgePathEroded
        edgePathRetained
        
        % The array particlePath stores the number of constrictions met by
        % each particle before being retained or eroded
        particlePathRetained
        particlePathEroded
        % The array particleDistance stores the distance travelled by each
        % particle before being retained or eroded
        particleDistanceRetained
        particleDistanceEroded
        % The array particleRoute stores the nodes that each particle
        % passed through (Only for 'storeOption' = 'storePath')
        particleRoute
        
        % The matrix layerCaptured stores the PSD of particles captured in
        % each layer
        layerRetained
        
        % The matrices timeEroded, timeCaptured, timeIntroduced store:
        % 1: Total number of particles eroded/captured/introduced
        % 2: Total volume of particles eroded/captured/introduced
        % 3: Largest particle eroded/captured/introduced
        % 4:end: Cumulative PSD of particles eroded/captured/introduced
        timeEroded
        timeIntroduced
        timeRetained
        timeOutside % Only for introductionMethod = fixedParticleNum2 (NOT USED)
        timeTotalFlow
    end
    
    methods
        % ===== Initialisation ============================================
        function [obj] = Statistics(nodeNum, edgeNum, nodeLattice, psdFineBucket, maxTimeStep)
            % Initialise statistics
            n = size(psdFineBucket, 1);
            
            % Initialise array(s)
            obj.nodeRetained = zeros(nodeNum, 6);
            obj.timeNodeRetained = struct;
            obj.nodeClogging = NaN(nodeNum, 4);
            
            obj.edgeFlowing = zeros(edgeNum, 5);
            obj.edgeFlowing(:,1) = 1;
            obj.edgePath = zeros(edgeNum, 1);
            obj.edgePathEroded = obj.edgePath;
            obj.edgePathRetained = obj.edgePath;
            
            obj.particlePathRetained = zeros(max(nodeLattice(:,end)), 3);
            obj.particlePathEroded = zeros(max(nodeLattice(:,end)), 1);
            obj.particleDistanceRetained = struct('mechanism1',[],'mechanism2',[],'mechanism3',[]);
            obj.particleDistanceEroded = [];
            
            obj.particleRoute = struct;
            
            obj.layerRetained = zeros(max(nodeLattice(:,end)), n);
            
            obj.timeEroded = zeros(maxTimeStep, n+3);
            obj.timeIntroduced = obj.timeEroded;
            obj.timeRetained = obj.timeEroded;
            obj.timeOutside = obj.timeEroded;
            
            obj.timeTotalFlow = zeros(maxTimeStep,1);
        end
        
        % ===== nodeRetained ==============================================
        % Determine the volume of particles retained in each layer
        % expressed as a percentage of the total node volume in that layer
        function [p] = percentageLayerVolume(obj, nodeLattice, nodeVolume, nodeInlet)
            % Number of layers
            layerNum = max(nodeLattice(:,end));
            
            % Determine the volume of particles retained in each layer and
            % the total volume in each layer
            particleVolTotal = zeros(layerNum,1);
            layerVolTotal = particleVolTotal;
            for i = 1:layerNum
                particleVolTotal(i) = sum(obj.nodeRetained(and(nodeLattice(:,end)==i, nodeInlet==0), 3));
                layerVolTotal(i) = sum(nodeVolume(and(nodeLattice(:,end)==i, nodeInlet==0)));
            end
            
            % Determine the percentages
            p = particleVolTotal./layerVolTotal;
        end
        
        % Plot the number of particles retained in each layer (number of
        % layers fixed at 20)
        function [] = plotNodeRetainedNum(obj, nodeLattice)
            % Number of layers
            layerNum = max(nodeLattice(:,end));
            
            % Determine the numbers of particles retained in each layer
            particleNumTotal = zeros(layerNum,1);
            particleNum1 = particleNumTotal;
            particleNum2 = particleNumTotal;
            particleNum3 = particleNumTotal;
            for i = 1:layerNum
                particleNumTotal(i) = sum(obj.nodeRetained(nodeLattice(:,end)==i,1));
                particleNum1(i) = sum(obj.nodeRetained(nodeLattice(:,end)==i,4));
                particleNum3(i) = sum(obj.nodeRetained(nodeLattice(:,end)==i,5));
                particleNum2(i) = particleNumTotal(i)-particleNum1(i)-particleNum3(i);
            end
            
            % Plot
            figure;
            set(gcf,'WindowState','maximized');
            plot(1:layerNum,particleNumTotal,'LineWidth',2);
            hold on;
            plot(1:layerNum,particleNum1,'LineWidth',2);
            plot(1:layerNum,particleNum2,'LineWidth',2);
            plot(1:layerNum,particleNum3,'LineWidth',2);
            legend({'Total number of particles retained','Particles retained due to particle straining', ...
                'Particles retained due to clogged exits','Particles retained due to filled node'}, ...
                'Location','northeast','FontName','Helvetica','FontSize',16);
            xlim([1,layerNum]);
            % ylim([0,1]);
            set(gca, 'FontName','Helvetica','FontSize',16);
            xlabel('Distance from inlets','FontName','Helvetica','FontSize',16);
            ylabel('Number of particles retained','FontName','Helvetica','FontSize',16);
            xticks(1:layerNum);
        end
        
        % Number of particles retained in each layer (variable number of
        % layers)
        function [n] = layerRetainedNum(obj, filter, layerNum, direction)
            % Initialise n
            n = zeros(layerNum,1);
            
            % Determine the nodes in each layer
            nodeLayer = filter.layerNode(layerNum, direction);
            
            % For each layer
            for i = 1:layerNum
                % Index of nodes in that layer
                idxNodes = find(nodeLayer==i);
                nodeNum = length(idxNodes);
                
                % For each node in layer i
                for j = 1:nodeNum
                    % Current node
                    node = idxNodes(j);
                    
                    % Number of particles retained in node j
                    n(i) = n(i)+obj.nodeRetained(node,1);
                end
            end
        end
        
        % Output array of boolean values representing whether each node is 
        % full
        function [nodeFull] = idxNodeFull(obj, d50, filter, densityLimit)
            nodeFull = filter.nodeVolume*densityLimit-obj.nodeRetained(:,3)<(pi/6)*(d50^3);
        end
        
        % Boolean value representing whether filer is seal
        function [bool] = inletFull(obj, filter, densityLimit)
            nodeFull = filter.nodeVolume*densityLimit<obj.nodeRetained(:,3);
            
            % If all inlets are full
            bool = all(nodeFull(filter.nodeInlet==1));
        end
        
        % Percentage of nodes full in layers (variable number of layers)
        function [p] = percentageLayerFull(obj, d50, filter, densityLimit, layerNum)
            % Determine the boundaries of layers
            b = filter.layerBoundary(layerNum, filter.networkParameters.flowDirection);
            
            % Determine the indices of nodes that are not inlet or outlet
            idxInternal = ~or(filter.nodeInlet, filter.nodeOutlet);
            
            % Initialise array(s)
            p = zeros(layerNum,1);
            
            % Direction considered
            switch filter.networkParameters.flowDirection(1)
                case 'x'
                    idx = 1;
                case 'y' 
                    idx = 2;
                case 'z'
                    idx = 3;
            end
            switch filter.networkParameters.flowDirection(2:4)
                case 'pos'
                    coef = 1;
                case 'neg'
                    coef = -1;
            end
            
            % Index of nodes that are full
            nodeFull = obj.idxNodeFull(d50, filter, densityLimit);
            
            % For each layer
            for i = 1:layerNum
                % Determine the index of nodes that are within layer i
                idxNode1 = coef*filter.nodeCoord(:,idx)>coef*b(i);
                idxNode2 = coef*filter.nodeCoord(:,idx)<=coef*b(i+1);
                idxNodeWithin = and(idxNode1, idxNode2);
                idxNodeWithin = and(idxNodeWithin, idxInternal); % Ignore inlets and outlets
                
                % Number of nodes
                nodeTotal = sum(idxNodeWithin);
                
                % Number of full nodes
                nodeFullTemp = and(nodeFull, idxNodeWithin);
                nodeFullTemp = sum(nodeFullTemp);
                
                % Percentage of full nodes
                p(i) = nodeFullTemp/nodeTotal*100;
            end
            assert(all(p<=100));
            assert(all(p>=0));
        end
        
        % Clusters of nodes that are full, as well as methods to plot
        % plotNode* = plot the number of nodes that are full within a
        %             region (not as good)
        % plotEdge* = plot edges connecting the nodes within the clusters
        function [bin, binsize] = monoNodeFullClusterSize(obj, d50, filter, densityLimit)
            % Number of nodes
            nodeNum = filter.nodeNum;
            
            % Determine the indices of nodes that are not inlet or outlet
            idxInternal = ~or(filter.nodeInlet, filter.nodeOutlet);
            
            % Determine the indices of nodes that are full
            boolFullNodes = and(idxInternal, filter.nodeVolume*densityLimit-obj.nodeRetained(:,3)<(pi/6)*(d50^3));
            idxFullNodes = find(boolFullNodes==1);
            
            % Determine the indices of edges that connect two full nodes
            idxFullEdges1 = ismember(filter.edgeList(:,1), idxFullNodes);
            idxFullEdges2 = ismember(filter.edgeList(:,2), idxFullNodes);
            idxFullEdges = and(idxFullEdges1, idxFullEdges2);
            
            % Create a Graph object
            G = graph(filter.edgeList(idxFullEdges==1,1), filter.edgeList(idxFullEdges==1,2));
            
            % Determine the clusters and the average size
            [bin, binsize] = conncomp(G);
            bin(~boolFullNodes) = 0;
            bin(length(bin):nodeNum) = 0;
            
            % Unique list of bin
            uniqueBin = unique(bin);
            uniqueBin = uniqueBin(2:end);
            uniqueNum = length(uniqueBin);
            
            % Sort bin and binsize
            for i = 1:uniqueNum
                idx = uniqueBin(i);
                
                bin(bin==idx) = i;
            end
            binsize = binsize(uniqueBin);
        end
        function [] = plotNodeFullCluster(obj, d50, filter, densityLimit)
            % Determine the indices of nodes that are not inlet or outlet
            idxInternal = ~or(filter.nodeInlet, filter.nodeOutlet);
            
            % Determine the indices of nodes that are full
            idxFullNodes = filter.nodeVolume*densityLimit-obj.nodeRetained(:,3)<(pi/6)*(d50^3);
            idxFullNodes = idxFullNodes(idxInternal==1);
            
            % Plot
            filter.plotSumNodeContour(idxFullNodes);
        end
        function [] = plotEdgeFullCluster(obj, d50, filter, densityLimit, varargin)
            % Input parser
            parser = inputParser;
            addParameter(parser, 'plotFirstFew', 3)
            parse(parser, varargin{:});
            plotFirstFew = parser.Results.plotFirstFew;
            
            % Determine the clusters
            [bin, binsize] = obj.monoNodeFullClusterSize(d50, filter, densityLimit);
            
            % Initialise array(s)
            idxEdge = zeros(filter.edgeNum,1);
            idxColour = idxEdge;
            
            % Plot the first few largest clusters
            for i = 1:plotFirstFew
                % Determine the index of the largest cluster
                maxBinSize = max(binsize);
                idxMax = find(binsize==maxBinSize,1);
                binsize(idxMax) = 0;
                
                % Nodes to be plotted
                idxNode = find(bin==idxMax);
                
                % Determine the edges to be plotted
                idxEdge1 = ismember(filter.edgeList(:,1),idxNode);
                idxEdge2 = ismember(filter.edgeList(:,2),idxNode);
                idxEdgeTemp = and(idxEdge1, idxEdge2);
                
                % Edges to be plotted
                idxEdge = or(idxEdge, idxEdgeTemp);
                idxColour(idxEdgeTemp) = i;
            end
            
            % Plot
            filter.plotEdges(idxEdge, idxColour);
        end
        
        % Spatial correlation between boolean values representing whether a
        % node is full
        function [c] = spatialCorrelationFullNodes(obj, d50, filter, densityLimit, varargin)
            % Input parser
            parser = inputParser;
            addParameter(parser, 'distanceUpTo', 3)
            parse(parser, varargin{:});
            distanceUpTo = parser.Results.distanceUpTo;
            
            % Determine the nodes that are full
            [idxFullNodes] = obj.idxNodeFull(d50, filter, densityLimit);
            idxFullNodes = double(idxFullNodes);
            
            % Initialise c
            c = zeros(distanceUpTo,1);
            
            % Determine the spatial correlation
            for i = 1:distanceUpTo
                c(i) = filter.spatialCorrelation(idxFullNodes,i);
            end
        end
        
        % ===== nodeClogging ==============================================
        % NOT FULLY DEVELOPED
        % Plotting the number of nodes that are clogged, categorised as: 
        % 1: any of the exits clogged due to particle straining
        % 2: all exits clogged due to node clogging
        function [nodeCloggedNum] = plotTimeNodeClogging(obj)
            % Number of time steps introducing particles
            introducingT = size(obj.timeIntroduced,1);
            
            % Number of time steps 
            overallT = size(obj.timeEroded,1);
            
            % Number of node
            nodeNum = size(obj.nodeClogging,1);
            
            % Initialise array
            % Column 1: number of open nodes
            % Column 2: number of nodes clogged due to 1
            % Column 3: number of nodes clogged due to 2
            nodeCloggedNum = zeros(overallT+1,3);
            nodeCloggedNum(1,1) = nodeNum;
            
            % For each time step
            for i = 1:overallT
                nodeCloggedNum(i+1,:) = nodeCloggedNum(i,:);
                
                % Index of nodes that are clogged at this time step
                idxCloggedCurrent = find(obj.nodeClogging(:,1)==i);
                cloggedCurrentNum = length(idxCloggedCurrent);
                
                if cloggedCurrentNum>0
                    % Number of nodes clogged categorised as (2)
                    nodeCloggedNum(i+1,3) = nodeCloggedNum(i+1,3)+sum(obj.nodeClogging(idxCloggedCurrent,3)==0);
                    
                    % Number of nodes clogged categorised as (1)
                    nodeCloggedNum(i+1,2) = nodeCloggedNum(i+1,2)+cloggedCurrentNum-sum(obj.nodeClogging(idxCloggedCurrent,3)==0);
                    
                    % Number of open nodes
                    nodeCloggedNum(i+1,1) = nodeNum-sum(nodeCloggedNum(i+1,2:3));
                end
            end
            
            % Scaled by total number of nodes
            nodeCloggedNum = nodeCloggedNum/nodeNum*100;
            
            % Plot
            figure
            set(gcf,'WindowState','maximized')
            for i = 1:3
                plot(0:overallT,nodeCloggedNum(:,i),'LineWidth',2)
                hold on
            end
            plot((introducingT+0.5)*ones(1,2),ylim,'LineWidth',2,'LineStyle','--','Color','k')
            legend({'Open nodes','Clogged nodes with due to particle straining', ...
                'Clogged nodes due to node clogging'}, ...
                'Location','northeast','FontName','Helvetica','FontSize',16);
            xlim([0,overallT]);
            set(gca, 'FontName','Helvetica','FontSize',16)
            xlabel('Time step t','FontName','Helvetica','FontSize',16)
            ylabel('Percentage of nodes','FontName','Helvetica','FontSize',16)
        end
        
        % ===== edgeFlowing ===============================================
        % Number of edges that are clogged at the end of simulation
        function [n] = edgeCloggedNum(obj)
            % Percentage of edges clogged
            n = sum(obj.edgeFlowing(:,1)==0);
        end
        
        % Percentage of edges that are clogged, options available for each
        % retaining mechanism
        function [p] = percentageEdgeClogged(obj, mechanism)
            % Number of edges
            edgeNum = size(obj.edgeFlowing,1);
            
            % Percentage of edges clogged
            p = sum(and(obj.edgeFlowing(:,1)==0, obj.edgeFlowing(:,2)==mechanism))/edgeNum*100;
        end
        
        % Plot the percentages of edge that are clogged (total, as well as
        % for each mechanism)
        function [] = plotTimeEdgeClogged(obj)
            % Number of edges
            edgeNum = size(obj.edgeFlowing,1);
            
            % Number of time steps
            n = size(obj.timeTotalFlow,1);
            
            % Initialise array(s)
            timeEdgeClogged = zeros(n,1);
            timeEdgeClogged1 = timeEdgeClogged;
            timeEdgeClogged3 = timeEdgeClogged;
            
            % Determine timeEdgeClogged
            timeEdgeClogged(1) = sum(obj.edgeFlowing(:,5)==1);
            timeEdgeClogged1(1) = sum(and(obj.edgeFlowing(:,5)==1,obj.edgeFlowing(:,2)==1));
            timeEdgeClogged3(1) = sum(and(obj.edgeFlowing(:,5)==1,obj.edgeFlowing(:,2)==3));
            for i = 2:n
                timeEdgeClogged(i) = sum(obj.edgeFlowing(:,5)==i)+timeEdgeClogged(i-1);
                timeEdgeClogged1(i) = sum(and(obj.edgeFlowing(:,5)==i,obj.edgeFlowing(:,2)==1))+timeEdgeClogged1(i-1);
                timeEdgeClogged3(i) = sum(and(obj.edgeFlowing(:,5)==i,obj.edgeFlowing(:,2)==3))+timeEdgeClogged3(i-1);
            end
            
            % Plot
            figure;
            set(gcf,'WindowState','maximized');
            plot(0:n,[0;timeEdgeClogged]/edgeNum*100,'LineWidth',2);
            hold on;
            plot(0:n,[0;timeEdgeClogged1]/edgeNum*100,'LineWidth',2,'LineStyle',':');
            plot(0:n,[0;timeEdgeClogged3]/edgeNum*100,'LineWidth',2,'LineStyle',':');
            legend({'All edges','Edges clogged by particle straining','Edges clogged by filled node'}, ...
                'Location','southeast','FontName','Helvetica','FontSize',16);
            xlim([0,n]);
            ylim([0,100]);
            set(gca, 'FontName','Helvetica','FontSize',16);
            xlabel('Time step, t','FontName','Helvetica','FontSize',16);
            ylabel('Percentage of edges clogged','FontName','Helvetica','FontSize',16);
        end
        
        % Plot the diagonal elements in fabric tensor to show anisotropy
        function [] = plotTimeFabricTensor(obj, edgeNormal, varargin)
            % Input parser
            parser = inputParser;
            expectedMechanism = [1,3];
            addParameter(parser, 'mechanism', 0, @(x)(any(expectedMechanism==x)));
            parse(parser, varargin{:});
            mechanism = parser.Results.mechanism;
            
            % Number of edges
            edgeNum = size(obj.edgeFlowing,1);
            
            % Number of time steps
            n = size(obj.timeTotalFlow,1);
            
            % Initialise nxx nyy and nzz
            p = zeros(n,3);
            
            % Determine nxx nyy nzz
            edgeClogged = zeros(edgeNum,1);
            for i = 1:n
                % Determine edges that are clogged
                if mechanism==0
                    edgeClogged = or(edgeClogged, obj.edgeFlowing(:,5)==i);
                else
                    edgeClogged = or(edgeClogged, and(obj.edgeFlowing(:,5)==i, obj.edgeFlowing(:,2)==mechanism));
                end
                
                % Determine the fabric tensor
                P = obj.fabricTensor(edgeNormal(edgeClogged==1,:));
                
                % Store nxx nyy nzz
                p(i,1) = P(1,1);
                p(i,2) = P(2,2);
                p(i,3) = P(3,3);
            end
            
            % Plot
            figure;
            set(gcf,'WindowState','maximized');
            plot(1:n,p(:,1),'LineWidth',2);
            hold on;
            plot(1:n,p(:,2),'LineWidth',2);
            plot(1:n,p(:,3),'LineWidth',2);
            legend({'n_x_x','n_y_y','n_z_z'}, ...
                'Location','northeast','FontName','Helvetica','FontSize',16);
            xlim([1,n]);
            ylim([0,0.6]);
            set(gca, 'FontName','Helvetica','FontSize',16);
            xlabel('Time step, t','FontName','Helvetica','FontSize',16);
            ylabel('Percentage of edges clogged','FontName','Helvetica','FontSize',16);
        end
        
        % To determine how the number of edges connected to each node
        % change with time
        % timeNodeClogged is a 2D array
        % row index     = time
        % column index  = number of edges clogged
        % element       = number of nodes with 'column index' number of edges
        %               clogged
        function [timeNodeClogged] = timeNodeClogged(obj, nodeConnectivity)
            % The maximum number of connectivity
            maxConnectivity = max(nodeConnectivity);
            
            % Number of time steps
            timeNum = size(obj.timeEroded,1);
            
            % Initialise array(s)
            timeNodeClogged = zeros(timeNum,maxConnectivity+1);
            nodeClogged = zeros(size(nodeConnectivity));
            
            for i = 1:timeNum
                % Determine the index of nodes clogged at t = i
                idxNodes = obj.edgeFlowing(obj.edgeFlowing(:,5)==i,3);
                
                % Number of affected nodes
                nodeNum = length(idxNodes);
                
                for j = 1:nodeNum
                    node = idxNodes(j);
                    
                    % Store the number of edges clogged
                    nodeClogged(node) = nodeClogged(node)+1;
                end
                
                % Determine number of nodes with specific number of edges clogged
                for j = 1:maxConnectivity+1
                    timeNodeClogged(i,j) = sum(nodeClogged==j-1);
                end
            end
        end
        function [] = plotTimeNodeClogged(obj, nodeConnectivity)
            % Evaluate timeNodeClogged
            [timeNodeClogged] = obj.timeNodeClogged(nodeConnectivity);
            
            % Node number
            nodeNum = length(nodeConnectivity);
            
            % The maximum number of connectivity
            maxConnectivity = max(nodeConnectivity);
            
            % Number of time steps
            timeNum = size(obj.timeEroded,1);
            
            % Adjust timeNodeClogged
            timeNodeClogged = timeNodeClogged/nodeNum*100;
            timeNodeClogged = [100,zeros(1,maxConnectivity); ...
                timeNodeClogged];
            
            % Legend list
            legendList = cell(maxConnectivity+1,1);
            for i = 1:maxConnectivity+1
                legendList{i} = sprintf('Number of nodes with %i edges clogged', i-1);
            end
            
            % Plot
            figure; 
            set(gcf,'WindowState','maximized');
            for i = 1:maxConnectivity+1
                plot(0:timeNum,timeNodeClogged(:,i),'LineWidth',2);
                hold on;
            end
            legend(legendList,'Location','northeast','FontName','Helvetica','FontSize',16);
            set(gca, 'FontName','Helvetica','FontSize',16);
            xlabel('Time step, t','FontName','Helvetica','FontSize',16);
            ylabel('Percentage of edges clogged','FontName','Helvetica','FontSize',16);
        end
        
        % Percentages of edges clogged in each layer
        function [p] = percentageLayerClogged(obj, filter, layerNum, mechanism)
            % Determine the boundaries of layers
            b = filter.layerBoundary(layerNum, filter.networkParameters.flowDirection);
            
            % Initialise array(s)
            p = zeros(layerNum,1);
            
            % Direction considered
            switch filter.networkParameters.flowDirection(1)
                case 'x'
                    idx = 1;
                case 'y' 
                    idx = 2;
                case 'z'
                    idx = 3;
            end
            switch filter.networkParameters.flowDirection(2:4)
                case 'pos'
                    coef = 1;
                case 'neg'
                    coef = -1;
            end
            
            % For each layer
            for i = 1:layerNum
                % Determine the index of nodes that are within layer i
                idxEdge1 = coef*filter.edgeCoord(:,idx)>coef*b(i);
                idxEdge2 = coef*filter.edgeCoord(:,idx)<=coef*b(i+1);
                idxEdgeWithin = and(idxEdge1, idxEdge2);
                
                % Number of nodes
                edgeTotal = sum(idxEdgeWithin);
                
                % Number of full nodes
                edgeCloggedTemp = and(obj.edgeFlowing(:,1)==0, ismember(obj.edgeFlowing(:,2),mechanism));
                edgeCloggedTemp = and(edgeCloggedTemp, idxEdgeWithin);
                edgeCloggedTemp = sum(edgeCloggedTemp);
                
                % Percentage of full nodes
                p(i) = edgeCloggedTemp/edgeTotal*100;
            end
            assert(all(p<=100));
            assert(all(p>=0));
        end
        
        % ===== edgePath ==================================================
        % Number of particles flowing through each edge (for small
        % networks only)
        function [] = plotPath(obj, edgeNum, nodeCoord, edgeList)
            % Initialise figure with colormap. Note that the colormap is
            % hardwired and must be changed here, if required.
            figure;
            set(gcf, 'WindowState', 'maximized');
            hold on;
            
            % Map geometric properties to plot/line spec properties. Note that
            % the minimum and maximum values are hardwired and must be changed
            % here if required.
            mapGeometryToPlot = @(x, ymin, ymax) (ymin + ((ymax - ymin)/(max(x) - min(x)))*(x - min(x)));
            edgeSizeVal = ones(edgeNum,1);
            if ~sum(obj.edgePath)
                edgeColorVal = 0.01*ones(edgeNum,1);
            else
                edgeColorVal = mapGeometryToPlot(obj.edgePath, 0.99, 0.01);
            end
            
            % The edges are represented as lines. The size of the line is
            % proportional to the throat area, while the grey-scale color of
            % the line is proportional to the conductance (i.e. darker
            % indicates higher conductance).
            fprintf('Plotting edges ');
            for i = 1:edgeNum
                plot3(nodeCoord(edgeList(i,:), 1), ...
                    nodeCoord(edgeList(i,:), 2), ...
                    nodeCoord(edgeList(i,:), 3), ...
                    'lineWidth', edgeSizeVal(i), 'Color', edgeColorVal(i)*ones(3,1));
                % text(mean(nodeCoord(edgeList(i,:), 1))-0.1, ...
                %     mean(nodeCoord(edgeList(i,:), 2))-0.1, ...
                %     mean(nodeCoord(edgeList(i,:), 3)), ...
                %     cellstr(num2str(floor(obj.edgePath(i)))), 'FontSize', 12);
                
                % Display progress
                if i~=1
                    fprintf(repmat('\b', 1, 7));
                end
                fprintf('%6.2f%%', i/edgeNum*100);
            end
            fprintf('\n');
            hold off;
            
            set(gca, 'DataAspectRatio', [1,1,1]);
            title('Path followed by particles', 'FontSize', 24);
            axis off;
        end
        
        % Weighted average of edgeDiameter, weighted with the number of
        % particles passing through each edge
        function [d] = weightedAverageEdgeDiameter(obj, edgeDiameter)
            % Determine the weighted average
            d = sum(obj.edgePath.*edgeDiameter)/sum(obj.edgePath);
        end
        
        % Plot all edges traversed by eroded particles (only when
        % 'storeEdgePath' option is used
        function [] = plotEdgePathEroded(obj, filter, varargin)
            % Input parser
            parser = inputParser;
            addParameter(parser, 'plotView', 1:3)
            parse(parser, varargin{:})
            plotView = parser.Results.plotView;
            
            % Plot
            filter.plotEdges(obj.edgePathEroded>0, ones(filter.edgeNum,1), 'plotView', plotView);
        end
        
        % ===== particlePath ==============================================
        % Number of constrictions met by base particles (options available
        % for eroded, or retained, as well as retained by each mechanism)
        function [n] = constrictionNum(obj, varargin)
            % Input parser
            parser = inputParser;
            expectedMechanism = {'eroded','retained','retained1','retained2','retained3'};
            addParameter(parser, 'mechanism', 'eroded',  @(x) any(validatestring(x, expectedMechanism)));
            parse(parser, varargin{:});
            mechanism = parser.Results.mechanism;
            
            % Maximum number of constrictions met (The index-1 represents 
            % the number of constrictions met, while the entry is the
            % number of base particles)
            maxErodedNum = size(obj.particlePathEroded,1)-1;
            maxRetainedNum = size(obj.particlePathRetained,1)-1;
            
            % Determine the average distance travelled
            switch mechanism
                case 'eroded'
                    n = sum((0:maxErodedNum)'.*obj.particlePathEroded)/sum(obj.particlePathEroded);
                case 'retained'
                    n = sum((0:maxRetainedNum)'.*sum(obj.particlePathRetained,2))/sum(sum(obj.particlePathRetained));
                case 'retained1'
                    n = sum((0:maxRetainedNum)'.*obj.particlePathRetained(:,1))/sum(obj.particlePathRetained(:,1));
                case 'retained2'
                    n = sum((0:maxRetainedNum)'.*obj.particlePathRetained(:,2))/sum(obj.particlePathRetained(:,2));
                case 'retained3'
                    n = sum((0:maxRetainedNum)'.*obj.particlePathRetained(:,3))/sum(obj.particlePathRetained(:,3));
            end
        end
        function [] = plotConstrictionNum(obj)
            % Maximum number of constrictions met (The index-1 represents 
            % the number of constrictions met, while the entry is the
            % number of base particles)
            maxErodedNum = size(obj.particlePathEroded,1)-1;
            maxRetainedNum = size(obj.particlePathRetained,1)-1;
            
            % Determine the mean
            meanEroded = sum((0:maxErodedNum)'.*obj.particlePathEroded)/sum(obj.particlePathEroded);
            % meanRetained = sum((0:maxRetainedNum)'.*sum(obj.particlePathRetained,2))/sum(sum(obj.particlePathRetained));
            meanRetained1 = sum((0:maxRetainedNum)'.*obj.particlePathRetained(:,1))/sum(obj.particlePathRetained(:,1));
            meanRetained2 = sum((0:maxRetainedNum)'.*obj.particlePathRetained(:,2))/sum(obj.particlePathRetained(:,2));
            meanRetained3 = sum((0:maxRetainedNum)'.*obj.particlePathRetained(:,3))/sum(obj.particlePathRetained(:,3));
            
            % Colourlist
            colourList = Constants.colourList;
            
            % Plot in a histogram
            figure;
            set(gcf,'WindowState','maximized');
            plot(0:maxErodedNum,obj.particlePathEroded,'LineWidth',2);
            hold on;
            plot(0:maxRetainedNum,obj.particlePathRetained(:,1),'LineWidth',2);
            plot(0:maxRetainedNum,obj.particlePathRetained(:,2),'LineWidth',2);
            plot(0:maxRetainedNum,obj.particlePathRetained(:,3),'LineWidth',2);
            plot(zeros(2,1),zeros(2,1),'LineWidth',1.5,'Color','k','LineStyle','--');
            plot(meanEroded*ones(2,1),ylim,'LineWidth',1.5,'LineStyle','--','Color',colourList(1,:));
            plot(meanRetained1*ones(2,1),ylim,'LineWidth',1.5,'LineStyle','--','Color',colourList(2,:));
            plot(meanRetained2*ones(2,1),ylim,'LineWidth',1.5,'LineStyle','--','Color',colourList(3,:));
            plot(meanRetained3*ones(2,1),ylim,'LineWidth',1.5,'LineStyle','--','Color',colourList(4,:));
            legend({'Particles eroded','Particles retained by particle straining', ...
                'Particles retained by clogged exits','Particles retained by filled node', ...
                'Average number of constrictions met'}, ...
                'Location','northeast','FontName','Helvetica','FontSize',16);
            xlim([0,max([maxErodedNum,maxRetainedNum])]);
            set(gca, 'FontName','Helvetica','FontSize',16);
            xlabel('Number of constrictions met by base particles','FontName','Helvetica','FontSize',16);
            ylabel('Number of base particles','FontName','Helvetica','FontSize',16);
        end
        
        % ===== particleDistance ==========================================
        % Distance travelled (options available for eroded, or retained, as 
        % well as retained by each mechanism)
        function [d] = averageDistance(obj, varargin)
            % Input parser
            parser = inputParser;
            expectedMechanism = {'eroded','retained','retained1','retained2','retained3'};
            addParameter(parser, 'mechanism', 'eroded',  @(x) any(validatestring(x, expectedMechanism)));
            parse(parser, varargin{:});
            mechanism = parser.Results.mechanism;
            
            % Determine the average distance travelled
            switch mechanism
                case 'eroded'
                    d = mean(obj.particleDistanceEroded);
                case 'retained'
                    d = mean([obj.particleDistanceRetained.mechanism1; ...
                        obj.particleDistanceRetained.mechanism2; ...
                        obj.particleDistanceRetained.mechanism3]);
                case 'retained1'
                    d = mean(obj.particleDistanceRetained.mechanism1);
                case 'retained2'
                    d = mean(obj.particleDistanceRetained.mechanism2);
                case 'retained3'
                    d = mean(obj.particleDistanceRetained.mechanism3);
            end
        end
        function [] = plotDistance(obj)
            % Determine the CDF
            cdfEroded = Statistics.createCDF(obj.particleDistanceEroded,1e4);
            cdfRetained.mechanism1 = Statistics.createCDF(obj.particleDistanceRetained.mechanism1,1e4);
            cdfRetained.mechanism2 = Statistics.createCDF(obj.particleDistanceRetained.mechanism2,1e4);
            cdfRetained.mechanism3 = Statistics.createCDF(obj.particleDistanceRetained.mechanism3,1e4);
            
            % Determine the mean distance travelled
            meanEroded = mean(obj.particleDistanceEroded);
            meanRetained.mechanism1 = mean(obj.particleDistanceRetained.mechanism1);
            meanRetained.mechanism2 = mean(obj.particleDistanceRetained.mechanism2);
            meanRetained.mechanism3 = mean(obj.particleDistanceRetained.mechanism3);
            
            % Colourlist
            colourList = Constants.colourList;
            
            % Plot
            figure;
            set(gcf,'WindowState','maximized');
            plot(cdfEroded(:,2),cdfEroded(:,1),'LineWidth',2);
            hold on;
            plot(cdfRetained.mechanism1(:,2),cdfRetained.mechanism1(:,1),'LineWidth',2);
            plot(cdfRetained.mechanism2(:,2),cdfRetained.mechanism2(:,1),'LineWidth',2);
            plot(cdfRetained.mechanism3(:,2),cdfRetained.mechanism3(:,1),'LineWidth',2);
            plot(zeros(2,1),zeros(2,1),'LineWidth',1.5,'Color','k','LineStyle','--');
            plot(meanEroded*ones(2,1),ylim,'LineWidth',1.5,'LineStyle','--','Color',colourList(1,:));
            plot(meanRetained.mechanism1*ones(2,1),ylim,'LineWidth',1.5,'LineStyle','--','Color',colourList(2,:));
            plot(meanRetained.mechanism2*ones(2,1),ylim,'LineWidth',1.5,'LineStyle','--','Color',colourList(3,:));
            plot(meanRetained.mechanism3*ones(2,1),ylim,'LineWidth',1.5,'LineStyle','--','Color',colourList(4,:));
            legend({sprintf('Particles eroded (Mean = %.2f)',meanEroded),sprintf('Particles retained by particle straining (Mean = %.2f)',meanRetained.mechanism1), ...
                sprintf('Particles retained by clogged exits (Mean = %.2f)',meanRetained.mechanism2),sprintf('Particles retained by filled node (Mean = %.2f)',meanRetained.mechanism3), ...
                'Average number of constrictions met'}, ...
                'Location','northeast','FontName','Helvetica','FontSize',16);
            set(gca, 'FontName','Helvetica','FontSize',16);
            xlabel('Number of constrictions met by base particles','FontName','Helvetica','FontSize',16);
            ylabel('Number of base particles','FontName','Helvetica','FontSize',16);
        end
        
        % ===== particleRoute =============================================
        % Actual paths taken by particles (only when 'storePath' is used)
        function [idxEroded] = isEroded(obj, nodeOutlet)
            % Determine if the particle in particleRoute is eroded
            
            % Number of particles
            particleNum = numel(fieldnames(obj.particleRoute));
            
            % Initialise idxEroded
            idxEroded = zeros(particleNum,1);
            
            % Determine if each particle has eroded
            for i = 1:particleNum
                nodeEnd = obj.particleRoute.(sprintf('particle%i',i))(end);
                if nodeOutlet(nodeEnd)
                    idxEroded(i) = 1;
                end
            end
        end
        function [idxEdge, idxColour] = isEdgeTraversed(obj, nodeNodeTo, nodeEdgeTo, idxParticle, plotFirstFew)
            % Determine the indices of edges traversed by particles in
            % particleRoute
            
            % Determine the number of edges
            edgeNum = size(obj.edgeFlowing,1);
            
            % Initialise idxEdge
            idxEdge = zeros(edgeNum,1);
            idxColour = idxEdge;
            
            % Number of particles
            idxParticle = find(idxParticle==1);
            idxParticle = idxParticle(1:plotFirstFew);
            particleNum = length(idxParticle);
            
            % Determine all edges that are traversed
            fprintf('Particle ');
            for i = 1:particleNum
                if i~=1
                    fprintf(repmat('\b', 1, 7));
                end
                fprintf('%6.2f%%', i/particleNum*100);
                particle = idxParticle(i);
                
                % Determine the number of edges
                nodeList = obj.particleRoute.(sprintf('particle%i',particle));
                nodeList = nodeList(:);
                nodeNum = length(nodeList);
                node = nodeList(1);
                
                for j = 2:nodeNum
                    % Determine nodeTo
                    nodeTo = nodeList(j);
                    
                    % Determine the idx of node to
                    idx = find(nodeNodeTo(node,:)==nodeTo,1);
                    idxEdge(nodeEdgeTo(node,idx)) = 1;
                    idxColour(nodeEdgeTo(node,idx)) = i;
                    
                    % Update node
                    node = nodeTo;
                end
            end
            fprintf('\n');
        end
        function [] = plotParticlePath(obj, filter, varargin)
            % Unpack parameters
            [nodeOutlet, nodeNodeTo, nodeEdgeTo] = Simulator.unpack(filter, 'nodeOutlet', 'nodeNodeTo', 'nodeEdgeTo');
            
            % Input parser
            parser = inputParser;
            expectedMechanism = {'eroded','retained'};
            addParameter(parser, 'mechanism', 'eroded',  @(x) any(validatestring(x, expectedMechanism)))
            addParameter(parser, 'plotFirstFew', 3)
            parse(parser, varargin{:});
            mechanism = parser.Results.mechanism;
            plotFirstFew = parser.Results.plotFirstFew;
            
            % Determine the indices of particles that are eroded or retained
            switch mechanism
                case 'eroded'
                    idxParticle = obj.isEroded(nodeOutlet);
                case 'retained'
                    idxParticle = ~obj.isEroded(nodeOutlet);
            end
            
            % Number of particles
            particleNum = sum(idxParticle);
            if particleNum<plotFirstFew
                plotFirstFew = particleNum;
            end
            
            % Determine indices of edges to plot
            [idxEdge, idxColour] = obj.isEdgeTraversed(nodeNodeTo, nodeEdgeTo, idxParticle, plotFirstFew);
            
            % Plot
            filter.plotEdges(idxEdge, idxColour);
        end
        function [r] = ratioDistance(obj, filter)
            % Unpack parameters
            [nodeOutlet, nodeNodeTo, nodeEdgeTo, edgeList, edgeLength] = Simulator.unpack(filter, 'nodeOutlet', 'nodeNodeTo', 'nodeEdgeTo', 'edgeList', 'edgeLength');
            
            % Determine the particles that are eroded
            idxParticle = obj.isEroded(nodeOutlet);
            idxParticle = find(idxParticle==1);
            particleNum = length(idxParticle);
            
            % Create a Graph object
            G = graph(edgeList(:,1),edgeList(:,2),edgeLength);
            
            % Initialise r
            r = zeros(particleNum,1);
            
            for i = 1:particleNum
                % Particle
                particle = idxParticle(i); 
                
                % Nodelist
                nodeList = obj.particleRoute.(sprintf('particle%i',particle));
                
                % Determine the beginning and end nodes
                nodeStart = nodeList(1);
                nodeEnd = nodeList(end);
                
                % Determine the shortest path
                [~,s] = shortestpath(G, nodeStart, nodeEnd);
                
                % Determine the actual distance
                d = 0;
                nodeNum = length(nodeList);
                for j = 1:nodeNum-1
                    node = nodeList(j);
                    nodeTo = nodeList(j+1);
                    idx = find(nodeNodeTo(node,:)==nodeTo,1);
                    idxEdge = nodeEdgeTo(node,idx);
                    d = d+edgeLength(idxEdge);
                end
                
                % Ratio
                r(i) = d/s;
            end
        end
        
        % ===== layerRetained =============================================
        % PSDs of particles in layers
        function [] = plotLayerRetainedPSD(obj, psdFineBucket, varargin)
            % Input parser
            parser = inputParser;
            addParameter(parser, 'plotEvery', 4);
            parse(parser, varargin{:});
            plotEvery = parser.Results.plotEvery;
            
            % Number of layers
            layerNum = size(obj.layerRetained,1);
            legendList = cell(layerNum,1);
            idxPlot = zeros(layerNum,1);
            
            % Plot
            figure;
            set(gcf,'WindowState','maximized');
            for i = 1:layerNum
                if or(i==1,mod(i,plotEvery)==0)
                    [PSD] = Statistics.createPSDVolBucket(obj.layerRetained(i,:),psdFineBucket);
                    semilogx(psdFineBucket,PSD,'LineWidth',2);
                    hold on;
                    legendList{i} = sprintf('Layer %i',i);
                    idxPlot(i) = 1;
                end
            end
            legend(legendList{idxPlot==1},'Location','southeast','FontName','Helvetica','FontSize',16);
            set(gca, 'FontName','Helvetica','FontSize',16);
            xlabel('Base particle diameter','FontName','Helvetica','FontSize',16);
            ylabel('Percentage smaller by volume','FontName','Helvetica','FontSize',16);
        end
        
        % Plot number of particles within a certain particle diameter range
        function [] = plotLayerRetainedVol(obj, psdFineBucket, varargin)
            % Input parser
            parser = inputParser;
            addParameter(parser, 'idxRange', [])
            addParameter(parser, 'layerNum', 10)
            parse(parser, varargin{:})
            idxRange = parser.Results.idxRange;
            layerNum = parser.Results.layerNum;
            
            if isempty(idxRange)
                % Number of ranges
                dividePSD = 5;
                
                % Dividing psd into specified number of ranges
                idxRange = round(linspace(1,length(psdFineBucket),dividePSD+1));
            else
                dividePSD = length(idxRange)-1;
            end
            
            % Dividing into layerNum
            assert(mod(size(obj.layerRetained,1),layerNum)==0)
            combineNum = size(obj.layerRetained,1)/layerNum;
            
            % Initialise array(s) and cell(s)
            legendList = cell(1,dividePSD);
            particleVol = zeros(layerNum,dividePSD);
            
            % Calculate the mean volume
            meanVolume = (pi/6)*(psdFineBucket.^3);
            meanVolume = meanVolume';
            
            % Assemble the legend list
            for j = 1:dividePSD
                % Upper and lower limits of particle size (Index)
                idxLower = idxRange(j);
                idxUpper = idxRange(j+1);
                
                % Legend list
                if j==1
                    legendList{j} = [sprintf('%.3f ',psdFineBucket(idxLower)), '\leq d^n \leq', sprintf(' %.3f',psdFineBucket(idxUpper))];
                else
                    legendList{j} = [sprintf('%.3f ',psdFineBucket(idxLower)), '< d^n \leq', sprintf(' %.3f',psdFineBucket(idxUpper))];
                end
            end
            
            % Determine the number of particles within each layer for each
            % range of size
            for i = 1:layerNum
                idxLowerLayer = (i-1)*combineNum+1;
                idxUpperLayer = i*combineNum;
                for j = 1:dividePSD
                    % Upper and lower limits of particle size (Index)
                    if j==1
                        idxLower = idxRange(j);
                    else
                        idxLower = idxRange(j)+1;
                    end
                    idxUpper = idxRange(j+1);
                    
                    % Aggregate the number of particles
                    particleVol(i,j) = sum(sum(obj.layerRetained(idxLowerLayer:idxUpperLayer,idxLower:idxUpper),1).*meanVolume(idxLower:idxUpper));
                end
            end
            
            % Plot properties
            colourList = Constants.colourList;
            
            % Plot
            figure
            set(gcf,'WindowState','maximized')
            for i = 1:dividePSD
                plot(particleVol(:,i),(1:layerNum)','LineWidth',3,'Color',colourList(i,:))
                hold on
            end
            ylim([0.5,10.5])
            yticks(1:10)
            legend(legendList,'Location','northeast','FontName','Helvetica','FontSize',36)
            set(gca,'FontName','Helvetica','FontSize',36)
            xlabel('Volume of base particles that are retained','FontName','Helvetica','FontSize',36)
            ylabel('Layer number (Flow direction)','FontName','Helvetica','FontSize',36)
        end
        
        % ===== time* =====================================================
        % Largest particle eroded
        function [d] = largestParticle(obj)
            % Extract the largest particle eroded
            d = obj.timeEroded(end,3);
        end
        
        % Largest particle eroded as d95b (Moraci et. al., 2012)
        % Provided better correlation with Kenney et. al. (1985)
        function [d] = largestParticle95(obj, psdFineBucket)
            % Determine PSD by vol
            psdVolEroded = Statistics.createPSDVolBucket(obj.timeEroded(end,4:end),psdFineBucket);
            
            if all(~isnan(psdVolEroded))
                % Determine the indices at which psdVolEroded is unique
                [~,idx,~] = unique(psdVolEroded);
                
                % Interpolate for d95b
                d = interp1(psdVolEroded(idx),psdFineBucket(idx),95);
            else
                d = 0;
            end
        end
        
        % Number of particles introduced for five parts, each part
        % representing a time duration of T_intro/5
        function [n] = particleIntroducedPartsNum(obj, partNum)
            % Number of time steps introducing particles
            maxFinesTimeStep = size(obj.timeIntroduced,1);
            
            % Boundary of each parts
            idxBoundary = zeros(partNum,1);
            for i = 1:partNum
                idxBoundary(i) = round(maxFinesTimeStep/partNum*i);
            end
            
            % Determine the number of particles
            n = zeros(partNum,1);
            for i = 1:partNum
                if i==1
                    n(i) = obj.timeIntroduced(idxBoundary(i));
                else
                    n(i) = obj.timeIntroduced(idxBoundary(i))-obj.timeIntroduced(idxBoundary(i-1));
                end
            end
        end
        
        % Percentage of particles eroded by number or by volume
        function [p, varargout] = percentageErodedNum(obj, nodeInlet, varargin)
            % Input parser
            parser = inputParser;
            addParameter(parser, 'overall', 0)
            addParameter(parser, 'previousStatistics', [])
            parse(parser, varargin{:});
            overall = parser.Results.overall;
            previousStatistics = parser.Results.previousStatistics;
            
            % Current nodeRetained
            if ~isempty(previousStatistics)
                currentNodeRetained = obj.nodeRetained-previousStatistics.nodeRetained;
            else
                currentNodeRetained = obj.nodeRetained;
            end
            
            % Remove particles retained in inlets to improve stability (it
            % is possible to have sink at inlets, thus particles retained
            % unnecessarily)
            particleInlet = sum(currentNodeRetained(nodeInlet==1,1));
            p = obj.timeEroded(end,1)/(obj.timeIntroduced(end,1)-particleInlet)*100;
            
            if overall
                varargout = cell(1,2);
                varargout{1} = obj.timeEroded(end,1);
                varargout{2} = obj.timeIntroduced(end,1)-particleInlet;
            end
        end
        function [p, varargout] = percentageErodedVol(obj, varargin)
            % Input parser
            parser = inputParser;
            addParameter(parser, 'overall', 0)
            parse(parser, varargin{:});
            overall = parser.Results.overall;
            
            % Remove particles retained in inlets to improve stability (it
            % is possible to have sink at inlets, thus particles retained
            % unnecessarily)
            p = obj.timeEroded(end,2)/(obj.timeEroded(end,2)+obj.timeRetained(end,2))*100;
            
            if overall
                varargout = cell(1,2);
                varargout{1} = obj.timeEroded(end,2);
                varargout{2} = obj.timeEroded(end,2)+obj.timeRetained(end,2);
            end
        end
        
        % Percentage (introduced,retained,eroded) with time
        function [] = plotTimeVol(obj)
            % Total volume of particles introduced
            volTotal = obj.timeIntroduced(end,2);
            
            % Number of time steps introducing fines
            nFine = size(obj.timeIntroduced,1);
            
            % Number of time steps
            n = size(obj.timeTotalFlow,1);
            
            % Plot
            figure;
            set(gcf,'WindowState','maximized');
            plot(0:n,[0;obj.timeIntroduced(:,2);ones(n-nFine,1)*obj.timeIntroduced(end,2)]/volTotal*100,'LineWidth',2);
            hold on;
            plot(0:n,[0;obj.timeRetained(:,2)]/volTotal*100,'LineWidth',2);
            plot(0:n,[0;obj.timeEroded(:,2)]/volTotal*100,'LineWidth',2);
            xlim([0,n]);
            ylim([0,100]);
            legend({'Base particles introduced','Base particles retained','Base particles eroded'}, ...
                'Location','southeast','FontName','Helvetica','FontSize',16);
            set(gca, 'FontName','Helvetica','FontSize',16);
            xlabel('Time step, t','FontName','Helvetica','FontSize',16);
            ylabel('Percentage of base particles','FontName','Helvetica','FontSize',16);
        end
        function [] = plotTimeNum(obj, nodeInlet, varargin)
            % Input parser
            parser = inputParser;
            addParameter(parser, 'previousStatistics', [])
            parse(parser, varargin{:});
            previousStatistics = parser.Results.previousStatistics;
            
            % Current nodeRetained
            if ~isempty(previousStatistics)
                currentNodeRetained = obj.nodeRetained-previousStatistics.nodeRetained;
            else
                currentNodeRetained = obj.nodeRetained;
            end
            
            % Total number of particles stucked in inlets
            numDiff = sum(currentNodeRetained(nodeInlet==1,1));
            numIntroduced = obj.timeIntroduced(:,1);
            numIntroduced = numIntroduced-numDiff;
            
            % Number of time steps introducing fines
            nFine = size(obj.timeIntroduced,1);
            
            % Number of time steps
            n = size(obj.timeTotalFlow,1);
            
            % Plot
            figure;
            set(gcf,'WindowState','maximized');
            plot(0:n,[0;numIntroduced;ones(n-nFine,1)*numIntroduced(end,1)],'LineWidth',2);
            hold on;
            plot(0:n,[0;obj.timeRetained(:,1)],'LineWidth',2);
            plot(0:n,[0;obj.timeEroded(:,1)],'LineWidth',2);
            xlim([0,n]);
            ylim([0,obj.timeIntroduced(end,1)-numDiff]);
            legend({'Base particles introduced','Base particles retained','Base particles eroded'}, ...
                'Location','southeast','FontName','Helvetica','FontSize',16);
            set(gca, 'FontName','Helvetica','FontSize',16);
            xlabel('Time step, t','FontName','Helvetica','FontSize',16);
            ylabel('Number of base particles','FontName','Helvetica','FontSize',16);
        end
        
        % Time variation of PSDs of particles (retained,eroded)
        function [] = plotPSD(obj, psdFineBucket)
            [PSDEroded] = Statistics.createPSDVolBucket(obj.timeEroded(end,4:end),psdFineBucket);
            [PSDIntroduced] = Statistics.createPSDVolBucket(obj.timeIntroduced(end,4:end),psdFineBucket);
            [PSDRetained] = Statistics.createPSDVolBucket(obj.timeRetained(end,4:end),psdFineBucket);
            
            % Plot
            figure;
            set(gcf,'WindowState','maximized');
            semilogx(psdFineBucket,PSDIntroduced,'LineWidth',2);
            hold on;
            semilogx(psdFineBucket,PSDEroded,'LineWidth',2);
            semilogx(psdFineBucket,PSDRetained,'LineWidth',2);
            legend({'Base particles introduced','Base particles eroded','Base particles retained'}, ...
                'Location','southeast','FontName','Helvetica','FontSize',16);
            set(gca, 'FontName','Helvetica','FontSize',16);
            xlabel('Base particle diameter','FontName','Helvetica','FontSize',16);
            ylabel('Percentage smaller by volume','FontName','Helvetica','FontSize',16);
        end
        function [] = plotTimeRetainedPSD(obj, psdFineBucket, varargin)
            % Input parser
            parser = inputParser;
            addParameter(parser, 'idxPlot', [])
            addParameter(parser, 'csdNum', [])
            addParameter(parser, 'psdVol', [])
            addParameter(parser, 'xl', [])
            parse(parser, varargin{:})
            idxPlot = parser.Results.idxPlot;
            csdNum = parser.Results.csdNum;
            psdVol = parser.Results.psdVol;
            xl = parser.Results.xl;
            
            % Number of time steps
            n = size(obj.timeTotalFlow,1);
            
            % Default idxPlot
            if isempty(idxPlot)
                idxPlot = [1,100:100:n,n];
            end
            
            % Initialise array(s)
            legendList = cell(n,1);
            
            % Plot
            figure;
            set(gcf,'WindowState','maximized');
            for i = idxPlot
                % [PSD] = Statistics.createPSDNumBucket(obj.timeRetained(i,4:end));
                [PSD] = Statistics.createPSDVolBucket(obj.timeRetained(i,4:end),psdFineBucket);
                semilogx(psdFineBucket,PSD,'LineWidth',3)
                hold on;
                legendList{i} = sprintf('Time step t = %i',i);
            end
            if ~isempty(csdNum)
                semilogx(csdNum(:,2),csdNum(:,1),'LineWidth',3,'Color','k','LineStyle','-.')
                legendList{n+1} = 'CSD by number';
                idxPlot(end+1) = n+1;
            end
            if ~isempty(psdVol)
                semilogx(psdVol(:,2),psdVol(:,1),'LineWidth',3,'Color','k','LineStyle',':')
                legendList{n+2} = 'PSD by volume';
                idxPlot(end+1) = n+2;
            end
            % if ~isempty(csdNum)
            %     D35c = interp1(csdNum(:,1),csdNum(:,2),35);
            %     semilogx(D35c*ones(1,2),ylim,'LineWidth',2,'Color','k','LineStyle','--')
            %     Dminc = min(csdNum(:,2));
            %     semilogx(Dminc*ones(1,2),ylim,'LineWidth',2,'Color','k','LineStyle','--')
            % end
            if ~isempty(xl)
                xlim(xl)
            end
            xticks([0.05,0.1,0.5,1,5,10])
            legend(legendList{idxPlot},'Location','northwest','FontName','Helvetica','FontSize',36);
            set(gca, 'FontName','Helvetica','FontSize',36);
            xlabel('Normalised diameter of base particles [-]','FontName','Helvetica','FontSize',36);
            ylabel('Percentage smaller','FontName','Helvetica','FontSize',36);
        end
        function [] = plotTimeErodedPSD(obj, psdFineBucket, varargin)
            % Input parser
            parser = inputParser;
            addParameter(parser, 'idxPlot', [])
            addParameter(parser, 'csdNum', [])
            addParameter(parser, 'psdVol', [])
            addParameter(parser, 'xl', [])
            parse(parser, varargin{:})
            idxPlot = parser.Results.idxPlot;
            csdNum = parser.Results.csdNum;
            psdVol = parser.Results.psdVol;
            xl = parser.Results.xl;
            
            % Number of time steps
            n = size(obj.timeTotalFlow,1);
            
            % Default idxPlot
            if isempty(idxPlot)
                idxPlot = [1,100:100:n,n];
            end
            
            % Initialise array(s)
            legendList = cell(n,1);
            
            % Plot
            figure;
            set(gcf,'WindowState','maximized');
            for i = idxPlot
                [PSD] = Statistics.createPSDVolBucket(obj.timeEroded(i,4:end),psdFineBucket);
                semilogx(psdFineBucket,PSD,'LineWidth',3);
                hold on;
                legendList{i} = sprintf('Time step t = %i',i);
            end
            if ~isempty(csdNum)
                semilogx(csdNum(:,2),csdNum(:,1),'LineWidth',3,'Color','k','LineStyle','-.')
                legendList{n+1} = 'CSD by number';
                idxPlot(end+1) = n+1;
            end
            if ~isempty(psdVol)
                semilogx(psdVol(:,2),psdVol(:,1),'LineWidth',3,'Color','k','LineStyle',':')
                legendList{n+2} = 'PSD by volume';
                idxPlot(end+1) = n+2;
            end
            if ~isempty(xl)
                xlim(xl)
            end
            xticks([0.05,0.1,0.5,1,5,10])
            legend(legendList{idxPlot},'Location','southeast','FontName','Helvetica','FontSize',36);
            set(gca, 'FontName','Helvetica','FontSize',36);
            xlabel('Normalised diameter of base particles [-]','FontName','Helvetica','FontSize',36);
            ylabel('Percentage smaller','FontName','Helvetica','FontSize',36);
        end
        
        % Variation of total flow with time
        function [] = plotTimeTotalFlow(obj, simulator)
            % Number of time steps
            n = size(obj.timeTotalFlow,1);
            
            % Determine the total flow without any particles
            % Unpack parameters
            [nodeNum, nodeConnectivity, nodeCoord,nodeEdgeTo, nodeInlet, nodeNodeTo, nodeOutlet, nodeSurfaceArea, nodeVolume] = Simulator.unpack(simulator.filter, 'nodeNum', 'nodeConnectivity', 'nodeCoord', 'nodeEdgeTo', 'nodeInlet', 'nodeNodeTo', 'nodeOutlet', 'nodeSurfaceArea', 'nodeVolume');
            [edgeNum, edgeArea, edgeLength, edgeList, edgeReduced] = Simulator.unpack(simulator.filter, 'edgeNum', 'edgeArea', 'edgeLength', 'edgeList', 'edgeReduced');
            [flowDirection, networkType] = Simulator.unpack(simulator.filter.networkParameters, 'flowDirection', 'networkType');
            [pressureEntry, pressureExit, flowMethod] = Simulator.unpack(simulator.flowSolverParameters, 'pressureEntry', 'pressureExit', 'flowMethod');
            nodeActive = zeros(nodeNum,3);
            edgeFlowingWithoutParticles = ones(edgeNum,1);
            nodeRetainedTemp = 0;
            
            % Determine local node pressures with Matlab and Python
            [~, nodeSource, ~] = Simulator.determinePressure(nodeNum, edgeNum, nodeActive, nodeConnectivity, nodeCoord,nodeEdgeTo, nodeInlet, nodeNodeTo, nodeOutlet, nodeRetainedTemp, nodeSurfaceArea, nodeVolume, edgeArea, edgeFlowingWithoutParticles, edgeLength, edgeList, edgeReduced, pressureEntry, pressureExit, flowMethod, flowDirection, networkType);
            initialTotalFlow = sum(nodeSource(nodeInlet==1));
            
            % Determine L and A
            switch flowDirection(1)
                case 'x'
                    idxFlowDirection = 1;
                case 'y'
                    idxFlowDirection = 2;
                case 'z'
                    idxFlowDirection = 3;
            end
            L = abs(max(nodeCoord(:,idxFlowDirection))-min(nodeCoord(:,idxFlowDirection)));
            idxOtherDirections = find(~ismember(1:3,idxFlowDirection)==1);
            A = 1;
            for l = 1:2
                A = A*(abs(max(nodeCoord(:,idxOtherDirections(l)))-min(nodeCoord(:,idxOtherDirections(l)))));
            end
            
            % Determine permeability
            dynamicViscosity = Constants.dynamicViscosity;
            initialPermeability = (initialTotalFlow*dynamicViscosity*L)/(A*(pressureEntry-pressureExit));
            permeability = (obj.timeTotalFlow*dynamicViscosity*L)/(A*(pressureEntry-pressureExit));
            
            % Plot
            figure;
            set(gcf,'WindowState','maximized');
            plot(0:n,[initialPermeability;permeability]/initialPermeability,'LineWidth',2);
            xlim([0,n]);
            ylim([0,1]);
            set(gca, 'FontName','Helvetica','FontSize',16);
            xlabel('Time step, t','FontName','Helvetica','FontSize',16);
            ylabel('Permeability, k/k_0 [-]','FontName','Helvetica','FontSize',16);
        end
    end
    
    methods (Static)
        % Create CDF given an array of data
        function [CDF] = createCDF(y,interval,varargin)
            % The function takes x as input, which is a vector of scalar values, and
            % output its CDF at points 0:1e-4:1
            
            % Input parser
            parser = inputParser;
            addParameter(parser, 'frequency', [])
            parse(parser, varargin{:})
            frequency = parser.Results.frequency;
            
            if isempty(frequency)
                [f,x] = ecdf(y);
            else
                [f,x] = ecdf(y,'Frequency',frequency);
            end
            CDF = zeros(interval+1,2);
            CDF(:,1) = 0:1/interval:1;
            CDF(:,2) = interp1(f,x,CDF(:,1));
            CDF(:,1) = CDF(:,1)*100;
        end
        
        % Convertion or creation of PSD
        function [psdNum] = convertPSDVolToNum(psdVol)
            % Determine psdNum
            
            % Convert PSD of fines by vol to number
            PSD_num_temp = zeros(length(psdVol)-1,2);
            for n = 1:(length(psdVol)-1)
                PSD_num_temp(n,1) = mean(psdVol(n:n+1,2));
                PSD_num_temp(n,2) = psdVol(n+1,1)-psdVol(n,1);
            end
            Vol_temp = PSD_num_temp(:,1).^3 * pi * (4/3);
            Num_temp = (PSD_num_temp(:,2)./Vol_temp) / (PSD_num_temp(1,2)./Vol_temp(1));
            
            % Num_temp = Vol_temp ./ min(Vol_temp);
            PP_num = 100 * cumsum(Num_temp) ./ sum(Num_temp);
            PP_num = [0; PP_num];
            PP_num(end) = 100;
            psdNum = zeros(10001,2);
            psdNum(:,1) = 0:0.01:100;
            for k = 1:10001
                ind = find(PP_num >= (k-1)/100, 1, 'first');
                psdNum(k,2) = psdVol(ind,2);
            end
        end
        function [psdArea] = convertPSDVolToArea(psdVol)
            % Determine psdArea
            
            % Convert PSD of fines by vol to area
            PSD_num_temp = zeros(length(psdVol)-1,2);
            for n = 1:(length(psdVol)-1)
                PSD_num_temp(n,1) = mean(psdVol(n:n+1,2));
                PSD_num_temp(n,2) = psdVol(n+1,1)-psdVol(n,1);
            end
            Vol_temp = (PSD_num_temp(:,1).^3 * pi * (1/6))./(PSD_num_temp(:,1).^2 * pi);
            Num_temp = (PSD_num_temp(:,2)./Vol_temp) / (PSD_num_temp(1,2)./Vol_temp(1));
            
            % Num_temp = Vol_temp ./ min(Vol_temp);
            PP_num = 100 * cumsum(Num_temp) ./ sum(Num_temp);
            PP_num = [0; PP_num];
            PP_num(end) = 100;
            psdArea = zeros(10001,2);
            psdArea(:,1) = 0:0.01:100;
            for k = 1:10001
                ind = find(PP_num >= (k-1)/100, 1, 'first');
                psdArea(k,2) = psdVol(ind,2);
            end
        end
        function [psdVol] = convertDataToPSDVol(data)
            % Must be a column vector
            if length(data)~=size(data,1)
                data = data';
            end
            
            % Determine psdVol
            psdVol = sortrows(data);
            vol_coarse = psdVol.^3 * pi * (1/6);
            psdVol(:,2) = psdVol;
            psdVol(:,1) = 100 * cumsum(vol_coarse)./sum(vol_coarse);
        end
        function [PSD] = createPSDVolBucket(psdNumber, psdFineBucket)
            PSD = zeros(length(psdNumber),1);
            PSD(1) = (1/6*pi*psdFineBucket(1)^3)*psdNumber(1);
            for i = 2:length(psdNumber)
                PSD(i) = PSD(i-1)+(1/6*pi*psdFineBucket(i)^3)*psdNumber(i);
            end
            PSD = PSD/PSD(end)*100;
        end
        function [PSD] = createPSDNumBucket(psdNumber)
            PSD = zeros(length(psdNumber),1);
            PSD(1) = psdNumber(1);
            for i = 2:length(psdNumber)
                PSD(i) = PSD(i-1)+psdNumber(i);
            end
            PSD = PSD/PSD(end)*100;
        end
        
        % Obtain the Terzaghi's ratio (Another method is implmented in
        % Simulator class)
        function [ratio] = TerzaghiRatio(simulator, run)
            % psdFineVol and psdCourseVol
            psdFineVol = simulator.base.(sprintf('base%i',run)).psdFineVol;
            psdCoarseVol = Statistics.convertDataToPSDVol(2*simulator.filter.particleRadius);
            
            % Determine Terzaghi's ratio
            d85 = interp1(psdFineVol(:,1),psdFineVol(:,2),85);
            D15 = interp1(psdCoarseVol(:,1),psdCoarseVol(:,2),15);
            ratio = D15/d85;
        end
        
        % Determine the fabric tensor given a list of vectors
        function [P] = fabricTensor(n, varargin)
            if ~isempty(varargin)
                scale = varargin{1};
                assert(size(scale,1)>size(scale,2))
                n = n.*scale;
            end
            
            % Number of vectors
            vectorNum = size(n,1);
            
            % Initialise array(s)
            P = zeros(3);
            
            % Determine the sum
            for i = 1:vectorNum
                % nx*nx
                P(1,1) = P(1,1)+n(i,1)^2;
                
                % ny*ny
                P(2,2) = P(2,2)+n(i,2)^2;
                
                % nz*nz
                P(3,3) = P(3,3)+n(i,3)^2;
                
                % nx*ny
                P(1,2) = P(1,2)+n(i,1)*n(i,2);
                P(2,1) = P(2,1)+n(i,2)*n(i,1);
                
                % ny*nz
                P(2,3) = P(2,3)+n(i,2)*n(i,3);
                P(3,2) = P(3,2)+n(i,3)*n(i,2);
                
                % nx*ny
                P(1,3) = P(1,3)+n(i,1)*n(i,3);
                P(3,1) = P(3,1)+n(i,3)*n(i,1);
            end
            
            % Divide by number of vector
            P = P/vectorNum;
        end
        
        % Comparing similarity of arrays; 'similarityIndex' is based on
        % Sufian et. al. (2018)
        function [r] = similarityBool(v1, v2)
            r = sum(v1==v2)/length(v1);
        end
        function [r] = similarityIndex(v1, v2)
            if and(sum(v1)==0,sum(v2)==0)
                r = 2;
            else
                % Average number of input
                num = mean([sum(v1),sum(v2)]);
                
                % Number of similar ones
                idxBothOne = and(v1,v2);
                
                % Measure of similarity
                r = sum(idxBothOne)/num;
            end
        end
        
        % Pythagoras theorem
        function [d] = pythagoras(v1,v2)
            d = sqrt(sum((v1-v2).^2));
        end
    end
end