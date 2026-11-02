classdef Filter
    % The main function of the class is to generate networks based on the
    % output from the DEM and the DTM
    %
    % =====================================================================
    % AVAILABLE DATA (see Constants class)
    %
    % cuCoarse (Default to 1.2)
    % density (Default to dense)
    % size (Default to [1,1,1], cut from [1,1,4]) (Size of [0,0,0] is for
    %     [1,1,1] without cutting)
    % triangulationMethod (Default to 'Del', i.e. wDTM)
    % =====================================================================
    % NETWORK PARAMETERS
    %
    % networkType (Default to actual)
    % assignmentMethod (Default to actual)
    % directedOrUndirected (Default to undirected)
    % [nx,ny,nz] (Default to such that the number of nodes is comparable)
    % flowDirection (Default to zpos, i.e. in positive z-direction)
    % =====================================================================
    
    properties
        % Network parameters
        filterParameters = struct('minCoarse', [], ...
            'cuCoarse', [], ...
            'density', [], ...
            'size', [], ...
            'triangulationMethod', []);
        networkParameters = struct('networkType', [], ...
            'assignmentMethod', [], ...
            'directedOrUndirected', [], ...
            'nx', [], 'ny', [], 'nz', [], ...
            'flowDirection', []);
        
        % Node properties
        nodeNum
        nodeBoundary
        nodeConnectivity
        nodeCoord
        nodeEdgeDirection
        nodeEdgeTo
        nodeEdgeToPos
        nodeInlet
        nodeLattice
        nodeNodeTo
        nodeOutlet
        nodeSurfaceArea
        nodeVolume
        
        % Edge properties
        edgeNum
        edgeArea
        edgeCoord
        edgeDiameter
        edgeLength
        edgeList
        edgeNormal
        edgeParticles
        edgeReduced
        
        % Particle properties
        particleNum
        particleCoord
        particleRadius
        
        % CSD(s)
        csdNum
        analyticalCSD
        
        % Boundary
        boundary
        
        % Indices of nodes and edges within the boundary (For smaller
        % networks only)
        idxNodesInside
        idxEdgesOutside
        
        % Shortest path
        shortestPathZpos = struct;
        
        % List of node pairs at a certain distance away
        nodeList
    end
    
    methods
        % ===== Creating Filter Objects ===================================
        % Creation of an object of this class, define parameters for
        % filters
        function [obj] = Filter(varargin)
            if ~isempty(varargin)
                obj.filterParameters.minCoarse = varargin{1};
                obj.filterParameters.cuCoarse = varargin{2};
                obj.filterParameters.density = varargin{3};
                obj.filterParameters.size = varargin{4};
                obj.filterParameters.triangulationMethod = varargin{5};
                obj.networkParameters.networkType = varargin{6};
                if strcmp(obj.networkParameters.networkType, 'actual')
                    obj.networkParameters.assignmentMethod = varargin{7};
                else
                    obj.networkParameters.assignmentMethod = 'random';
                end
                nxyz = varargin{8};
                if length(nxyz)==3
                    obj.networkParameters.nx = nxyz(1);
                    obj.networkParameters.ny = nxyz(2);
                    obj.networkParameters.nz = nxyz(3);
                end
                obj.networkParameters.flowDirection = varargin{9};
                obj.networkParameters.directedOrUndirected = varargin{10};
            end
        end
        
        % Load data from Delaunay triangulation method (DTM)
        function [obj] = loadDTM(obj, loadMethod, varargin)
            % Main steps:
            %  (1) If processed data is available, skip to (7)
            %  (2) Load raw data output from DTM
            %  (3) Remove unconnected nodes [see method 
            %      'removeUnconnectedNodes']
            %  (4) Normalise all geometric properties by such that D0F = 1 
            %      [see method 'normaliseGeometricProperties']
            %  (5) Calculate constriction sizes [see method
            %      'calculateEdgeDiameter']
            %  (6) Calculate the reduced areas for all edges [see method 
            %      'calculateEdgeReduced']
            %  (7) Segment filter from larger sizes, this is to ensure the
            %      inlets and outlets are defined correctly at zlo and zhi
            %      [see method 'segmentFilter']
            %  (8) Pre-generate analytical CSD as proposed by Silveira et. 
            %      al. (1975) Locke et. al. (2001) [see method 
            %      'analyticalCSDDensest' and 'analyticalCSDLoosest']
            %  (9) Normalise all geometric properties by minCoarse (default  
            %      to 1) [see method 'normaliseGeometricProperties']
            % (10) Determine the CSD by number [see method
            %      'determineCSDNum']
            
            % Options for calculating edgeDiameter:
            % originalRadius: Use the original radius (or projected radius)
            %                 instead of the radius of the circle
            %                 intersecting the plane
            % moreParticles: Consider particles within a larger boundary
            %                when checking the whether the inscribed circle
            %                is within filter void space
            % smallestConstriction: Use the smallest valid (valid meaning
            %                       within filter void space) instead of
            %                       the largest one
            % visualiser: Use for debugging, visualise to ensure correct
            %             implementation (Must couple with idxEdges to plot
            %             selected edges only)
            % idxEdges: Use with visualiser, define the edges to be plotted
            parser = inputParser;
            addParameter(parser, 'originalRadius', 0)
            addParameter(parser, 'moreParticles', 1)
            addParameter(parser, 'smallestConstriction', 0)
            addParameter(parser, 'visualiser', 0)
            addParameter(parser, 'idxEdges', []);
            parse(parser, varargin{:})
            originalRadius = parser.Results.originalRadius;
            moreParticles = parser.Results.moreParticles;
            smallestConstriction = parser.Results.smallestConstriction;
            visualiser = parser.Results.visualiser;
            idxEdges = parser.Results.idxEdges;
            
            % Filter to be segment from
            currentSize = obj.filterParameters.size;
            obj.filterParameters.size = Constants.cutFrom(currentSize);
            
            if and(strcmp(loadMethod, 'read'), isfile(FileLocation.filterDTMCleaned(obj)))
                % Load cleaned data if available
                obj = obj.loadState(FileLocation.filterDTMCleaned(obj), ...
                    'variables', {'nodeNum', 'nodeCoord', 'nodeBoundary', 'nodeVolume', 'nodeSurfaceArea', ...
                    'edgeNum', 'edgeArea', 'edgeCoord', 'edgeDiameter', 'edgeLength', 'edgeList', 'edgeNormal', 'edgeParticles', 'edgeReduced', ...
                    'particleNum', 'particleCoord', 'particleRadius', ...
                    'boundary'});
            else
                % Load raw dtm data from dataFolder
                nodeData = load(FileLocation.filterDTMRaw(obj,'node'));
                edgeData = load(FileLocation.filterDTMRaw(obj,'edge'));
                particleData = load(FileLocation.filterDTMRaw(obj,'particle'));
                nodeData = nodeData.(sprintf('%s',FileLocation.filterDTMRawVariable(obj,'node')));
                edgeData = edgeData.(sprintf('%s',FileLocation.filterDTMRawVariable(obj,'edge')));
                particleData = particleData.(sprintf('%s',FileLocation.filterDTMRawVariable(obj,'particle')));
                
                % Store data
                obj.nodeCoord = nodeData.centroid;
                obj.nodeVolume = nodeData.pore_volume;
                obj.nodeBoundary = nodeData.boundary;
                obj.nodeSurfaceArea = nodeData.surface_area;
                obj.edgeList = edgeData.index;
                obj.edgeParticles = edgeData.face_index;
                obj.edgeArea = edgeData.pore_area;
                obj.edgeLength = edgeData.length;
                obj.particleCoord = particleData(:,size(particleData,2)-3:size(particleData,2)-1);
                obj.particleRadius = particleData(:,size(particleData,2));
                
                % Read the boundary coordinates
                s = fileread(FileLocation.filterBoundarySize(obj));
                idxStart = regexp(s,']');
                idxEnd = regexp(s,'N');
                idxEnd = idxEnd(end);
                s = s(idxStart+1:idxEnd-1);
                s = regexp(s,'\s','split');
                obj.boundary = zeros(6,1);
                j = 0;
                for i = 1:length(s)
                    if str2double(s{i})>0
                        j = j+1;
                        obj.boundary(j) = str2double(s{i});
                    end
                end
                
                % Determine nodeNum, edgeNum and particleNum
                obj.nodeNum = size(obj.nodeCoord,1);
                obj.edgeNum = size(obj.edgeList,1);
                obj.particleNum = size(obj.particleCoord,1);
                
                % Remove unconnected nodes
                obj = obj.removeUnconnectedNodes();
                
                % Normalise geometric properties
                obj = obj.normaliseGeometricProperties(1/(2*min(obj.particleRadius)));
                
                % Calculate edgeDiameter and edgeReduced
                obj = obj.calculateEdgeDiameter(originalRadius, moreParticles, smallestConstriction, visualiser, idxEdges);
                obj = obj.calculateEdgeReduced(visualiser, idxEdges);
                
                % Save cleaned data
                if ~visualiser
                    obj.saveState(FileLocation.filterDTMCleaned(obj), ...
                        'variables', {'nodeNum', 'nodeCoord', 'nodeBoundary', 'nodeVolume', 'nodeSurfaceArea', ...
                        'edgeNum', 'edgeArea', 'edgeCoord', 'edgeDiameter', 'edgeLength', 'edgeList', 'edgeNormal', 'edgeParticles', 'edgeReduced', ...
                        'particleNum', 'particleCoord', 'particleRadius', ...
                        'boundary'});
                end
            end
            
            % Revert to required size
            obj.filterParameters.size = currentSize;
            
            % Cut for smaller filters
            availableData = Constants.availableSize;
            if ~any(all(availableData==currentSize, 2))
                obj = obj.segmentFilter(loadMethod);
            end
            
            % Determine analytical CSD
            if and(strcmp(loadMethod, 'read'), isfile(FileLocation.filterAnalyticalCSD(obj)))
                % Load cleaned data if available
                obj = obj.loadState(FileLocation.filterAnalyticalCSD(obj), ...
                    'variables', {'analyticalCSD'});
            else
                % Generate the analytical CSDs
                obj = obj.analyticalCSDDensest;
                obj = obj.analyticalCSDLoosest;
                
                % Save data
                obj.saveState(FileLocation.filterAnalyticalCSD(obj), ...
                    'variables', {'analyticalCSD'});
            end
            
            % Normalise geometric properties (default to D0F)
            obj = obj.normaliseGeometricProperties(obj.filterParameters.minCoarse);
            
            % Determine CSD by number
            obj = obj.determineCSDNum();
        end
        
        % Generate network
        function [obj] = generateNetwork(obj)
            % Generate network, either based on actual void space structure
            % [see method 'generateNetworkTopologyActual'] or regular 
            % lattices with the geometric properties randomly distributed
            % [see method 'generateNetworkTopologyRegular']
            
            % Generate network topology
            if strcmp(obj.networkParameters.networkType, 'actual')
                obj = obj.generateNetworkTopologyActual();
            else
                obj = obj.generateNetworkTopologyRegular();
            end
            
            % Assign network geometry
            if strcmp(obj.networkParameters.assignmentMethod, 'random')
                [obj] = assignNetworkGeometry(obj);
            end
            
            % Shortest path
            if strcmp(obj.networkParameters.networkType, 'actual')
                if strcmp(obj.networkParameters.flowDirection, 'zpos')
                    obj = obj.shortestPathActual();
                end
            end
        end
        
        % Generate network topology (regular)
        function [obj] = generateNetworkTopologyRegular(obj)
            if isfile(FileLocation.filterNetwork(obj))
                obj = obj.loadState(FileLocation.filterNetwork(obj), ...
                    'variables', {'nodeNum', 'edgeNum', 'nodeInlet', 'nodeLattice', 'nodeOutlet', 'edgeList', ...
                    'nodeConnectivity', 'nodeEdgeTo', 'nodeNodeTo', 'nodeEdgeDirection'});
                
                % Update nodeEdgeDirection
                if and(isempty(obj.nodeEdgeDirection), strcmp(obj.networkParameters.directedOrUndirected,'undirected'))
                    % Generate nodeEdgeDirection
                    obj = obj.generateNodeEdgeDirection();
                    
                    % Save
                    obj.saveState(FileLocation.filterNetwork(obj), ...
                        'variables', {'nodeNum', 'edgeNum', 'nodeInlet', 'nodeLattice', 'nodeOutlet', 'edgeList', ...
                        'nodeConnectivity', 'nodeEdgeTo', 'nodeNodeTo', 'nodeEdgeDirection'});
                end
            else
                switch obj.networkParameters.networkType
                    case 'square'
                        [obj.nodeNum, obj.edgeNum, obj.edgeList, obj.nodeEdgeTo, obj.nodeNodeTo, ...
                            obj.nodeConnectivity, obj.nodeLattice, obj.nodeInlet, obj.nodeOutlet] ...
                            = generateSquareNetwork(obj.networkParameters.nx, obj.networkParameters.nz, obj.networkParameters.directedOrUndirected, 1);
                    case 'slantedsquare'
                        [obj.nodeNum, obj.edgeNum, obj.edgeList, obj.nodeEdgeTo, obj.nodeNodeTo, ...
                            obj.nodeConnectivity, obj.nodeLattice, obj.nodeInlet, obj.nodeOutlet] ...
                            = generateSlantedSquareNetwork(obj.networkParameters.nx, obj.networkParameters.nz, obj.networkParameters.directedOrUndirected, 1);
                    case 'cubic'
                        [obj.nodeNum, obj.edgeNum, obj.edgeList, obj.nodeEdgeTo, obj.nodeNodeTo, ...
                            obj.nodeConnectivity, obj.nodeLattice, obj.nodeInlet, obj.nodeOutlet] ...
                            = generateCubicNetwork(obj.networkParameters.nx, obj.networkParameters.ny, obj.networkParameters.nz, obj.networkParameters.directedOrUndirected, 1);
                    case 'slantedcubic'
                        [obj.nodeNum, obj.edgeNum, obj.edgeList, obj.nodeEdgeTo, obj.nodeNodeTo, ...
                            obj.nodeConnectivity, obj.nodeLattice, obj.nodeInlet, obj.nodeOutlet] ...
                            = generateSlantedCubicNetwork(obj.networkParameters.nx, obj.networkParameters.ny, obj.networkParameters.nz, obj.networkParameters.directedOrUndirected, 1);
                    case 'bcc'
                        [obj.nodeNum, obj.edgeNum, obj.edgeList, obj.nodeEdgeTo, obj.nodeNodeTo, ...
                            obj.nodeConnectivity, obj.nodeLattice, obj.nodeInlet, obj.nodeOutlet] ...
                            = generateBCCNetwork(obj.networkParameters.nx, obj.networkParameters.ny, obj.networkParameters.nz, obj.networkParameters.directedOrUndirected, 1);
                    case 'fcc'
                        [obj.nodeNum, obj.edgeNum, obj.edgeList, obj.nodeEdgeTo, obj.nodeNodeTo, ...
                            obj.nodeConnectivity, obj.nodeLattice, obj.nodeInlet, obj.nodeOutlet] ...
                            = generateFCCNetwork(obj.networkParameters.nx, obj.networkParameters.ny, obj.networkParameters.nz, obj.networkParameters.directedOrUndirected, 1);
                    case 'squarenp'
                        [obj.nodeNum, obj.edgeNum, obj.edgeList, obj.nodeEdgeTo, obj.nodeNodeTo, ...
                            obj.nodeConnectivity, obj.nodeLattice, obj.nodeInlet, obj.nodeOutlet] ...
                            = generateSquareNetwork(obj.networkParameters.nx, obj.networkParameters.nz, obj.networkParameters.directedOrUndirected, 0);
                    case 'slantedsquarenp'
                        [obj.nodeNum, obj.edgeNum, obj.edgeList, obj.nodeEdgeTo, obj.nodeNodeTo, ...
                            obj.nodeConnectivity, obj.nodeLattice, obj.nodeInlet, obj.nodeOutlet] ...
                            = generateSlantedSquareNetwork(obj.networkParameters.nx, obj.networkParameters.nz, obj.networkParameters.directedOrUndirected, 0);
                    case 'cubicnp'
                        [obj.nodeNum, obj.edgeNum, obj.edgeList, obj.nodeEdgeTo, obj.nodeNodeTo, ...
                            obj.nodeConnectivity, obj.nodeLattice, obj.nodeInlet, obj.nodeOutlet] ...
                            = generateCubicNetwork(obj.networkParameters.nx, obj.networkParameters.ny, obj.networkParameters.nz, obj.networkParameters.directedOrUndirected, 0);
                    case 'slantedcubicnp'
                        [obj.nodeNum, obj.edgeNum, obj.edgeList, obj.nodeEdgeTo, obj.nodeNodeTo, ...
                            obj.nodeConnectivity, obj.nodeLattice, obj.nodeInlet, obj.nodeOutlet] ...
                            = generateSlantedCubicNetwork(obj.networkParameters.nx, obj.networkParameters.ny, obj.networkParameters.nz, obj.networkParameters.directedOrUndirected, 0);
                    case 'bccnp'
                        [obj.nodeNum, obj.edgeNum, obj.edgeList, obj.nodeEdgeTo, obj.nodeNodeTo, ...
                            obj.nodeConnectivity, obj.nodeLattice, obj.nodeInlet, obj.nodeOutlet] ...
                            = generateBCCNetwork(obj.networkParameters.nx, obj.networkParameters.ny, obj.networkParameters.nz, obj.networkParameters.directedOrUndirected, 0);
                    case 'fccnp'
                        [obj.nodeNum, obj.edgeNum, obj.edgeList, obj.nodeEdgeTo, obj.nodeNodeTo, ...
                            obj.nodeConnectivity, obj.nodeLattice, obj.nodeInlet, obj.nodeOutlet] ...
                            = generateFCCNetwork(obj.networkParameters.nx, obj.networkParameters.ny, obj.networkParameters.nz, obj.networkParameters.directedOrUndirected, 0);
                end
                
                if strcmp(obj.networkParameters.directedOrUndirected,'undirected')
                    % Generate nodeEdgeDirection
                    obj = obj.generateNodeEdgeDirection();
                end
                
                % Save to improve efficiency
                obj.saveState(FileLocation.filterNetwork(obj), ...
                    'variables', {'nodeNum', 'edgeNum', 'nodeInlet', 'nodeLattice', 'nodeOutlet', 'edgeList', ...
                    'nodeConnectivity', 'nodeEdgeTo', 'nodeNodeTo', 'nodeEdgeDirection'});
            end
        end
        
        % Assign constant geometric properties (Only used for validation)
        function [obj] = assignConstantConductance(obj, edgeDiameter)
            obj.nodeSurfaceArea = ones(obj.nodeNum,1);
            obj.nodeVolume = ones(obj.nodeNum,1);
            obj.edgeArea = ones(obj.edgeNum,1);
            obj.edgeDiameter = ones(obj.edgeNum,1)*edgeDiameter;
            obj.edgeLength = ones(obj.edgeNum,1);
            obj.edgeReduced = ones(obj.edgeNum,1);
            obj.nodeCoord = obj.nodeLattice;
        end
        
        % ===== Statistics ================================================
        % Determine the percolation threshold, either: 
        % (1) Assuming geometric properties are randomly distributed or
        %     with spatial correlation
        % (2) With directed or undirected network
        function [p] = percolationThreshold(obj, varargin)
            % Input parser
            parser = inputParser;
            addParameter(parser, 'minP', 0);
            addParameter(parser, 'every', max([1,floor(obj.edgeNum/1e4)]));
            addParameter(parser, 'display', 'off');
            parse(parser, varargin{:});
            minP = parser.Results.minP;
            every = parser.Results.every;
            display = parser.Results.display;
            
            % Initialise array(s)
            isPercolating = 0;
            isFirstLoop = 1;
            edgeOpen = rand(obj.edgeNum,1)<minP;
            idxClosed = find(edgeOpen==0);
            edgeOpenNum = sum(edgeOpen);
            
            % Indices and numbers of inlets and outlets
            idxInlet = find(obj.nodeInlet==1);
            idxOutlet = find(obj.nodeOutlet==1);
            inletNum = sum(obj.nodeInlet);
            outletNum = sum(obj.nodeOutlet);
            
            if strcmp(display, 'on')
                fprintf('Determining percolation threshold ');
            end
            while ~isPercolating
                % Randomly open an edge
                randNum = ceil(rand(1)*length(idxClosed));
                idxToBeOpened = idxClosed(randNum);
                idxClosed(randNum) = [];
                edgeOpen(idxToBeOpened) = 1;
                
                % Determine the percolation threshold
                edgeOpenNum = edgeOpenNum+1;
                p = edgeOpenNum/obj.edgeNum;
                
                if and(p>=minP, mod(edgeOpenNum,every)==0)
                    % Create a Graph object
                    G = graph(obj.edgeList(edgeOpen==1,1),obj.edgeList(edgeOpen==1,2),[],obj.nodeNum);
                    
                    % Initialise
                    isLastInlet = 0;
                    idxInletCurrent = idxInlet;
                    inletNumCurrent = inletNum;
                    
                    % Determine if the network is percolating
                    while ~isLastInlet
                        inlet = idxInletCurrent(1);
                        
                        % Perform a breadth-first search
                        [v] = bfsearch(G, inlet);
                        
                        % If the network percolates
                        for j = 1:outletNum
                            outlet = idxOutlet(j);
                            if any(v==outlet)
                                isPercolating = 1;
                                break;
                            end
                        end
                        
                        % Remove inlet if it is in v
                        idxRemove = zeros(inletNumCurrent,1);
                        for j = 1:inletNumCurrent
                            inlet = idxInletCurrent(j);
                            if any(v==inlet)
                                idxRemove(j) = 1;
                            end
                        end
                        idxInletCurrent(idxRemove==1) = [];
                        inletNumCurrent = length(idxInletCurrent);
                        
                        % Break if network is percolating
                        if isPercolating
                            break;
                        end
                        
                        % Is this the last inlet
                        if inletNumCurrent==0
                            isLastInlet = 1;
                        end
                    end
                end
                
                % Display progress
                if strcmp(display, 'on')
                    if ~isFirstLoop
                        fprintf(repmat('\b', 1, 6));
                    else
                        isFirstLoop = 0;
                    end
                    fprintf('%6.4f', p);
                end
            end
            if strcmp(display, 'on')
                fprintf('\n');
            end
        end
        function [p] = percolationThresholdAssigned(obj, varargin)
            % Input parser
            parser = inputParser;
            addParameter(parser, 'minP', 0);
            addParameter(parser, 'display', 'off');
            parse(parser, varargin{:});
            minP = parser.Results.minP;
            display = parser.Results.display;
            
            % Initialisation
            p = minP;
            increment = 1e-4;
            isPercolating = 0;
            isFirstLoop = 1;
            CSD = obj.csdNum;
            
            % Indices and numbers of inlets and outlets
            idxInlet = find(obj.nodeInlet==1);
            idxOutlet = find(obj.nodeOutlet==1);
            inletNum = sum(obj.nodeInlet);
            outletNum = sum(obj.nodeOutlet);
            
            if strcmp(display, 'on')
                fprintf('Determining percolation threshold ');
            end
            while ~isPercolating
                % Determine the edge diameter corresponding to p
                d = interp1(CSD(:,1),CSD(:,2),(1-p)*100);
                
                % Determine edges that are open for a certain p
                edgeOpen = obj.edgeDiameter>d;
                
                % Create a Graph object
                G = graph(obj.edgeList(edgeOpen==1,1),obj.edgeList(edgeOpen==1,2),[],obj.nodeNum);
                
                % Initialise
                isLastInlet = 0;
                idxInletCurrent = idxInlet;
                inletNumCurrent = inletNum;
                
                % Determine if the network is percolating
                while ~isLastInlet
                    inlet = idxInletCurrent(1);
                    
                    % Perform a breadth-first search
                    [v] = bfsearch(G, inlet);
                    
                    % If the network percolates
                    for j = 1:outletNum
                        outlet = idxOutlet(j);
                        if any(v==outlet)
                            isPercolating = 1;
                            break;
                        end
                    end
                    
                    % Remove inlet if it is in v
                    idxRemove = zeros(inletNumCurrent,1);
                    for j = 1:inletNumCurrent
                        inlet = idxInletCurrent(j);
                        if any(v==inlet)
                            idxRemove(j) = 1;
                        end
                    end
                    idxInletCurrent(idxRemove==1) = [];
                    inletNumCurrent = length(idxInletCurrent);
                    
                    % Break if network is percolating
                    if isPercolating
                        break;
                    end
                    
                    % Is this the last inlet
                    if inletNumCurrent==0
                        isLastInlet = 1;
                    end
                end
                
                % Increase p
                if ~isPercolating
                    p = p+increment;
                end
                
                % Display progress
                if strcmp(display, 'on')
                    if ~isFirstLoop
                        fprintf(repmat('\b', 1, 6));
                    else
                        isFirstLoop = 0;
                    end
                    fprintf('%6.4f', p);
                end
            end
            if strcmp(display, 'on')
                fprintf('\n');
            end
        end
        function [p, distance] = percolationThresholdAssignedDirected(obj, idxPositive, idxNegative, varargin)
            % Input parser
            parser = inputParser;
            addParameter(parser, 'minP', 0);
            addParameter(parser, 'display', 'off');
            parse(parser, varargin{:});
            minP = parser.Results.minP;
            display = parser.Results.display;
            
            % Initialisation
            p = minP;
            increment = 1e-4;
            isPercolating = 0;
            isFirstLoop = 1;
            CSD = obj.csdNum;
            
            % Indices and numbers of inlets and outlets
            idxInlet = find(obj.nodeInlet==1);
            idxOutlet = find(obj.nodeOutlet==1);
            inletNum = sum(obj.nodeInlet);
            outletNum = sum(obj.nodeOutlet);
            
            if strcmp(display, 'on')
                fprintf('Determining percolation threshold ');
            end
            while ~isPercolating
                % Determine the edge diameter corresponding to p
                d = interp1(CSD(:,1),CSD(:,2),(1-p)*100);
                
                % Determine edges that are open for a certain p
                edgeOpen = obj.edgeDiameter>d;
                idxPositiveCurrent = and(edgeOpen, idxPositive);
                idxNegativeCurrent = and(edgeOpen, idxNegative);
                
                % Create a Graph object
                G = digraph([obj.edgeList(idxPositiveCurrent==1,1);obj.edgeList(idxNegativeCurrent==1,2)], ...
                    [obj.edgeList(idxPositiveCurrent==1,2);obj.edgeList(idxNegativeCurrent==1,1)], ...
                    [obj.edgeLength(idxPositiveCurrent==1);obj.edgeLength(idxNegativeCurrent==1)],obj.nodeNum);
                
                % Initialise
                isLastInlet = 0;
                idxInletCurrent = idxInlet;
                inletNumCurrent = inletNum;
                
                % Determine if the network is percolating
                while ~isLastInlet
                    inlet = idxInletCurrent(1);
                    
                    % Perform a breadth-first search
                    [v] = bfsearch(G, inlet);
                    
                    % If the network percolates
                    for j = 1:outletNum
                        outlet = idxOutlet(j);
                        if any(v==outlet)
                            isPercolating = 1;
                            [~,distance] = shortestpath(G,inlet,outlet);
                            break;
                        end
                    end
                    
                    % Remove inlet if it is in v
                    idxRemove = zeros(inletNumCurrent,1);
                    for j = 1:inletNumCurrent
                        inlet = idxInletCurrent(j);
                        if any(v==inlet)
                            idxRemove(j) = 1;
                        end
                    end
                    idxInletCurrent(idxRemove==1) = [];
                    inletNumCurrent = length(idxInletCurrent);
                    
                    % Break if network is percolating
                    if isPercolating
                        break;
                    end
                    
                    % Is this the last inlet
                    if inletNumCurrent==0
                        isLastInlet = 1;
                    end
                end
                
                % Increase p
                if ~isPercolating
                    p = p+increment;
                end
                
                % Display progress
                if strcmp(display, 'on')
                    if ~isFirstLoop
                        fprintf(repmat('\b', 1, 6));
                    else
                        isFirstLoop = 0;
                    end
                    fprintf('%6.4f', p);
                end
            end
            if strcmp(display, 'on')
                fprintf('\n');
            end
        end
        
        % Determine Dxc* and DxF given x
        function [D] = percentileEdgeDiameter(obj, p)
            D = interp1(obj.csdNum(:,1),obj.csdNum(:,2),p);
        end
        function [D] = percentileParticleDiameter(obj, p)
            % Create CDF for particle diameter
            psdCoarseVol = Statistics.convertDataToPSDVol(2*obj.particleRadius);
            
            % Determine the correponding D
            D = interp1(psdCoarseVol(:,1),psdCoarseVol(:,2),p);
        end
        function [D] = percentileParticleDiameterVol(obj, p)
            % PSD by vol
            psdCoarseVol = Statistics.convertDataToPSDVol(2*obj.particleRadius);
            
            % Interpolation
            D = interp1(psdCoarseVol(:,1),psdCoarseVol(:,2),p);
        end
        
        % Determine the shortest path from inlets to outlets
        function [P, d] = shortestPath(obj, varargin)
            % Input parser
            parser = inputParser;
            expectedWeight = {'length','number'};
            addParameter(parser, 'weight', 'length',  @(x) any(validatestring(x, expectedWeight)));
            parse(parser, varargin{:});
            weight = parser.Results.weight;
            
            if numel(fieldnames(obj.shortestPathZpos))==0
                [P,d] = obj.determineShortestPath(weight);
            else
                P = obj.shortestPathZpos.P;
                d = obj.shortestPathZpos.d;
            end
        end
        function [] = plotShortestPath(obj, varargin)
            % Only for zpos
            assert(strcmp(obj.networkParameters.flowDirection, 'zpos'));
            
            % Input parser
            parser = inputParser;
            expectedWeight = {'length','number'};
            addParameter(parser, 'weight', 'length',  @(x) any(validatestring(x, expectedWeight)));
            parse(parser, varargin{:});
            weight = parser.Results.weight;
            
            [P,~] = obj.shortestPath('weight',weight);
            
            % Number of edges
            plotNum = length(P)-1;
            
            % Plot in three directions
            xlabelList = {'x','y','x'};
            ylabelList = {'z (Flow direction)','z (Flow direction)','y'};
            plotList = [1,3; 2,3; 1,2];
            for i = 1:3
                figure;
                set(gcf,'WindowState','maximized');
                for j = 1:plotNum
                    plot([obj.nodeCoord(P(j),plotList(i,1)),obj.nodeCoord(P(j+1),plotList(i,1))], ...
                        [obj.nodeCoord(P(j),plotList(i,2)),obj.nodeCoord(P(j+1),plotList(i,2))], ...
                        'LineWidth',Constants.lineWidth,'Color','k');
                    hold on;
                end
                xlim([floor(min(obj.nodeCoord(:,plotList(i,1)))),ceil(max(obj.nodeCoord(:,plotList(i,1))))]);
                ylim([floor(min(obj.nodeCoord(:,plotList(i,2)))),ceil(max(obj.nodeCoord(:,plotList(i,2))))]);
                set(gca, 'DataAspectRatio', [1,1,1])
                set(gca, 'FontName','Helvetica','FontSize',16);
                xlabel(xlabelList{i},'FontName','Helvetica','FontSize',16);
                ylabel(ylabelList{i},'FontName','Helvetica','FontSize',16);
            end
        end
        
        % Analytical CSD proposed by Silveira et. al. 1975 Locke et. al.
        % 2001 Indraratna et. al. (2007)
        function [csd] = analyticalCSDCombined(obj, particleNum, Rd)
            % Rd is a fraction
            assert(and(Rd>=0, Rd<=1))
            
            if particleNum==4
                % Loosest state
                csd = obj.analyticalCSD.loosest;
            elseif particleNum==3
                % Three particle configuration, combine with loosest state
                % if Rd<1
                if Rd==1
                    csd = obj.analyticalCSD.densest;
                else
                    % Equation 6 in Indraratna et. al. 2007
                    csd = obj.analyticalCSD.densest;
                    csd(:,2) = csd(:,2)+(csd(:,1)/100).*(1-Rd).*(obj.analyticalCSD.loosest(:,2)-csd(:,2));
                end
            end
        end
        
        % Calculate total length given the indices of nodes using node
        % coordinates
        function [d] = pathLength(obj, idxNodes)
            % Number of nodes
            pathNodeNum = length(idxNodes);
            
            % Initialise d
            d = 0;
            
            for i = 1:pathNodeNum-1
                % Determine the length
                d = d+Statistics.pythagoras(obj.nodeCoord(idxNodes(i),:),obj.nodeCoord(idxNodes(i+1),:));
            end
        end
        
        % Given list of monodisperse base particle diameters, determine how 
        % many particles can fit in each node
        function [m] = meanNormalisedVolume(obj, d)
            % Number of bases
            baseNum = length(d);
            
            % Initialise array
            m = zeros(baseNum,1);
            
            % Determine the proportion of nodes
            for i = 1:baseNum
                % Volume of base particle
                baseVol = (pi/6)*(d(i)^3);
                
                % Normalise nodeVolume by the volume of one base particle
                nodeVolumeNormalised = obj.nodeVolume*0.5/baseVol;
                
                % Determine the min and max
                m(i) = mean(nodeVolumeNormalised);
            end
        end
        function [p] = proportionVolume(obj, d)
            % Number of bases
            baseNum = length(d);
            
            % Initialise struct
            p = struct;
            
            % Determine the proportion of nodes
            for i = 1:baseNum
                % Volume of base particle
                baseVol = (pi/6)*(d(i)^3);
                
                % Normalise nodeVolume by the volume of one base particle
                nodeVolumeNormalised = obj.nodeVolume*0.5/baseVol;
                
                % Determine the min and max
                minNorm = floor(min(nodeVolumeNormalised));
                maxNorm = ceil(max(nodeVolumeNormalised));
                
                % Initialise array
                if minNorm==0
                    minNorm = 1;
                end
                p.(sprintf('d%i',round(d(i)*1000))) = zeros(maxNorm-minNorm,2);
                
                for j = minNorm:maxNorm-1
                    if j==1
                        p.(sprintf('d%i',round(d(i)*1000)))(j,1) = 2;
                        p.(sprintf('d%i',round(d(i)*1000)))(j,2) = sum(nodeVolumeNormalised<2);
                    else
                        p.(sprintf('d%i',round(d(i)*1000)))(j,1) = j+1;
                        p.(sprintf('d%i',round(d(i)*1000)))(j,2) = sum(and(nodeVolumeNormalised>=j, nodeVolumeNormalised<j+1));
                    end
                end
                
                % Obtain percentage
                assert(sum(p.(sprintf('d%i',round(d(i)*1000)))(:,2))==obj.nodeNum);
                p.(sprintf('d%i',round(d(i)*1000)))(:,2) = p.(sprintf('d%i',round(d(i)*1000)))(:,2)/obj.nodeNum*100;
            end
        end
        function [] = plotProportionVolume(obj, d)
            % Number of bases
            baseNum = length(d);
            
            % Determine the percentages
            [p] = obj.proportionVolume(d);
            
            % Plot
            figure; 
            set(gcf,'WindowState','maximized');
            for i = 1:baseNum
                plot(p.(sprintf('d%i',round(d(i)*1000)))(:,1),p.(sprintf('d%i',round(d(i)*1000)))(:,2));
                hold on;
            end
        end
        function [] = plotHistProportionVolume(obj, d)
            % Number of bases
            baseNum = length(d);
            
            % Determine the percentages
            [p] = obj.proportionVolume(d);
            
            for i = 1:baseNum
                % Plot
                figure;
                set(gcf,'WindowState','maximized');
                histogram('BinEdges',[p.(sprintf('d%i',round(d(i)*1000)))(1,1)-1,p.(sprintf('d%i',round(d(i)*1000)))(:,1)'], ...
                    'BinCounts',p.(sprintf('d%i',round(d(i)*1000)))(:,2))
            end
        end
        function [] = plotMeanNormalisedVolume(obj, d)
            % Determine the mean normalised node volume
            [m] = obj.meanNormalisedVolume(d);
            
            % Plot
            figure;
            set(gcf,'WindowState','maximized');
            plot(d,m);
        end
        
        % Average edgeDiameter in each layer (Not useful as the average is
        % the same)
        function [d] = layerEdgeDiameter(obj)
            assert(strcmp(obj.networkParameters.networkType, 'actual'));
            
            % Number of layers
            layerNum = max(obj.nodeLattice(:,end));
            
            % Initialise d
            d = zeros(layerNum,1);
            
            % Consider the flow direction
            switch obj.networkParameters.flowDirection(1)
                case 'x'
                    idx = 1;
                case 'y'
                    idx = 2;
                case 'z'
                    idx = 3;
            end
            
            % Minimum and maximum nodeCoord
            minCoord = min(obj.nodeCoord(:,idx));
            maxCoord = max(obj.nodeCoord(:,idx));
            
            % Create bucket
            coordBucket = linspace(minCoord,maxCoord,layerNum+1);
            
            for i = 1:layerNum
                % Edges in each layer
                idxEdge = and(obj.edgeCoord(:,idx)>coordBucket(i), obj.edgeCoord(:,idx)<coordBucket(i+1));
                
                % Average diameter
                d(i) = mean(obj.edgeDiameter(idxEdge==1));
            end
        end
        
        % Min and max nodeCoord for inlets and outlets (Used for debugging,
        % i.e. whether inlets and outlets are defined correctly)
        function [m] = minMaxCoord(obj)
            % Initialise m
            m = zeros(4,1);
            
            % Determine the min and max
            m(1) = min(obj.nodeCoord(obj.nodeInlet==1,3));
            m(2) = max(obj.nodeCoord(obj.nodeInlet==1,3));
            m(3) = min(obj.nodeCoord(obj.nodeOutlet==1,3));
            m(4) = max(obj.nodeCoord(obj.nodeOutlet==1,3));
        end
        
        % Size of clusters with nodes having largest edge diameter smaller
        % than d
        function [s] = nodeSmaller(obj, d)
            % Determine the indices of nodes that are not inlet or outlet
            idxNodes = ~or(obj.nodeInlet, obj.nodeOutlet);
            idxNodes = find(idxNodes==1);
            
            % Initialise idxNodeSmaller
            idxNodeSmaller = zeros(length(idxNodes),1);
            
            % Determine the largest edgeDiameter
            for i = 1:length(idxNodes)
                node = idxNodes(i);
                
                % Determine the indices of edges
                n = obj.nodeConnectivity(node);
                idxEdges = obj.nodeEdgeTo(node,1:n);
                
                % The largest diameter
                largestEdgeDiameter = max(obj.edgeDiameter(idxEdges));
                
                if largestEdgeDiameter<d
                    idxNodeSmaller(i) = 1;
                end
            end
            
            % Nodes that are smaller than d
            idxNodes = idxNodes(idxNodeSmaller==1);
            
            % Determine the indices of edges that connect two smaller nodes
            idxEdges1 = ismember(obj.edgeList(:,1), idxNodes);
            idxEdges2 = ismember(obj.edgeList(:,2), idxNodes);
            idxEdges = and(idxEdges1, idxEdges2);
            
            % Create a Graph object
            G = graph(obj.edgeList(idxEdges==1,1), obj.edgeList(idxEdges==1,2));
                        
            % Determine the clusters and the average size
            [~, binsize] = conncomp(G);
            binsize = binsize(binsize~=1);
            binsize = [binsize, ones(1,length(idxNodes)-sum(binsize))];
            assert(length(idxNodes)==sum(binsize));
            
            % Output
            s = binsize;
        end
        
        % Determine the boundaries of layers (Used for plotting, e.g.
        % determining how a specific property varies with depth)
        function [b] = layerBoundary(obj, layerNum, direction)
            % Note that (1,2,3) represents (x,y,z) directions recpectively
            % Determine the index corresponding to each direction
            switch direction(1)
                case 'x'
                    idx = ones(1,2)*1;
                case 'y'
                    idx = ones(1,2)*2;
                case 'z'
                    idx = ones(1,2)*3;
            end
            switch direction(2:4)
                case 'pos'
                    idx(1) = idx(1)*2-1;
                    idx(2) = idx(2)*2;
                case 'neg'
                    idx(1) = idx(1)*2;
                    idx(2) = idx(2)*2-1;
            end
            
            % The boundaries
            lo = obj.boundary(idx(1));
            hi = obj.boundary(idx(2));
            
            % Output the boundaries of each layer
            b = linspace(lo,hi,layerNum+1);
        end
        function [nodeLayer] = layerNode(obj, layerNum, direction)
            % Initialise nodeLayer
            nodeLayer = zeros(obj.nodeNum,1);
            
            % Determine the boundaries
            b = obj.layerBoundary(layerNum, direction);
            
            % Index to consider in nodeCoord
            switch direction(1)
                case 'x'
                    idxCoord = 1;
                case 'y'
                    idxCoord = 2;
                case 'z'
                    idxCoord = 3;
            end
            switch direction(2:4)
                case 'pos'
                    coef = 1;
                case 'neg'
                    coef = -1;
            end
            
            % For each layer
            for i = 1:layerNum
                % For each node
                for j = 1:obj.nodeNum
                    % Store the index of layer
                    if and(coef*obj.nodeCoord(j,idxCoord)>coef*b(i), coef*obj.nodeCoord(j,idxCoord)<coef*b(i+1))
                        nodeLayer(j) = i;
                    end
                end
            end
        end
        
        % Two-point correlation given a node property
        function [r] = spatialCorrelation(obj, nodeProperty, distance)
            % Node list
            nodeComb = obj.nodeList.(sprintf('distance%i',distance));
            
            % Determine the correlation
            r = corrcoef(nodeProperty(nodeComb(:,1)),nodeProperty(nodeComb(:,2)));
            r = r(1,2);
        end
        
        % Correlation of node volume
        function [c] = spatialCorrelationNodeVolume(obj, varargin)
            % Input parser
            parser = inputParser;
            addParameter(parser, 'distanceUpTo', 3)
            parse(parser, varargin{:});
            distanceUpTo = parser.Results.distanceUpTo;
            
            % Initialise c
            c = zeros(distanceUpTo,1);
            
            % Determine the spatial correlation
            for i = 1:distanceUpTo
                c(i) = obj.spatialCorrelation(obj.nodeVolume,i);
            end
        end
        
        % ===== Plotting ==================================================
        % Methods to plot the average or sum of a node- or an 
        % edge-associated variable (Only for actual networks)
        % Inlets and outlets ignored
        function [] = plotAverageNodeContour(obj, nodeVariable, varargin)
            % Input parser
            parser = inputParser;
            addParameter(parser, 'layerNum', 10)
            addParameter(parser, 'gridNum', 50)
            addParameter(parser, 'layerDirection', 1)
            addParameter(parser, 'idxPlotLayer', [])
            parse(parser, varargin{:})
            layerNum = parser.Results.layerNum;
            gridNum = parser.Results.gridNum;
            layerDirection = parser.Results.layerDirection;
            idxPlotLayer = parser.Results.idxPlotLayer;
            if isempty(idxPlotLayer)
                idxPlotLayer = 1:layerNum;
            end
            
            % Initialise array
            averageVariable = zeros(layerNum,gridNum,gridNum);
            boundaryLayer = zeros(layerNum+1,1);
            boundaryGrid1 = zeros(gridNum+1,1);
            boundaryGrid2 = zeros(gridNum+1,1);
            midLayer = zeros(layerNum,1);
            midGrid1 = zeros(gridNum,1);
            midGrid2 = zeros(gridNum,1);
            % direction = {'x','y','z'};
            
            % Determine the indices of nodes that are not inlet or outlet
            idxInternal = ~or(obj.nodeInlet, obj.nodeOutlet);
            
            % The boundary in layerDirection
            lowerLayerBoundary = obj.boundary(layerDirection*2-1);
            upperLayerBoundary = obj.boundary(layerDirection*2);
            
            % The boundary in perpendicular directions
            perpendicularDirection = 1:3;
            perpendicularDirection = perpendicularDirection(perpendicularDirection~=layerDirection);
            lowerGrid1 = obj.boundary(perpendicularDirection(1)*2-1);
            upperGrid1 = obj.boundary(perpendicularDirection(1)*2);
            lowerGrid2 = obj.boundary(perpendicularDirection(2)*2-1);
            upperGrid2 = obj.boundary(perpendicularDirection(2)*2);
            
            % Determine the boundaries and mid
            boundaryLayer(1) = lowerLayerBoundary;
            boundaryGrid1(1) = lowerGrid1;
            boundaryGrid2(1) = lowerGrid2;
            for i = 1:layerNum
                boundaryLayer(i+1) = lowerLayerBoundary+(upperLayerBoundary-lowerLayerBoundary)*(i)/layerNum;
                midLayer(i) = mean(boundaryLayer(i:i+1));
            end
            for i = 1:gridNum
                boundaryGrid1(i+1) = lowerGrid1+(upperGrid1-lowerGrid1)*(i)/gridNum;
                boundaryGrid2(i+1) = lowerGrid2+(upperGrid2-lowerGrid2)*(i)/gridNum;
                midGrid1(i) = mean(boundaryGrid1(i:i+1));
                midGrid2(i) = mean(boundaryGrid2(i:i+1));
            end
            
            for i = idxPlotLayer
                % Determine the boundary of the layer
                boundary1 = boundaryLayer(i);
                boundary2 = boundaryLayer(i+1);
                
                % Determine the nodes in layer i
                idxLayer = and(obj.nodeCoord(idxInternal==1,layerDirection)>boundary1, obj.nodeCoord(idxInternal==1,layerDirection)<=boundary2);
                
                for j = 1:gridNum
                    % Determine the nodes in grid1 j
                    boundary1 = boundaryGrid1(j);
                    boundary2 = boundaryGrid1(j+1);
                    
                    % Determine the nodes in grid1 j
                    idxGrid1 = and(obj.nodeCoord(idxInternal==1,perpendicularDirection(1))>boundary1, obj.nodeCoord(idxInternal==1,perpendicularDirection(1))<=boundary2);
                    
                    for k = 1:gridNum
                        % Determine the nodes in grid2 k
                        boundary1 = boundaryGrid2(k);
                        boundary2 = boundaryGrid2(k+1);
                        
                        % Determine the nodes in grid1 j
                        idxGrid2 = and(obj.nodeCoord(idxInternal==1,perpendicularDirection(2))>boundary1, obj.nodeCoord(idxInternal==1,perpendicularDirection(2))<=boundary2);
                        
                        % Nodes in layer i, grid1 j, grid2 k
                        idxNodes = and(idxLayer, idxGrid1);
                        idxNodes = and(idxNodes, idxGrid2);
                        
                        % Determine the average
                        averageVariable(i,j,k) = mean(nodeVariable(idxNodes==1));
                    end
                end
            end
            
            % Fill missing values
            % averageVariable = fillmissing(averageVariable,'linear');
            % averageVariable(averageVariable<0) = 0;
            averageVariable(isnan(averageVariable)) = 0;
            
            % Plot contour plots
            for i = idxPlotLayer
                figure;
                set(gcf,'WindowState','maximized');
                contourf(midGrid1, midGrid2, reshape(averageVariable(i,:,:),gridNum,gridNum));
                colorbar;
                set(gca, 'DataAspectRatio', [1,1,1])
            end
        end
        function [] = plotSumNodeContour(obj, nodeVariable, varargin)
            % Input parser
            parser = inputParser;
            addParameter(parser, 'layerNum', 10)
            addParameter(parser, 'gridNum', 50)
            addParameter(parser, 'layerDirection', 1)
            addParameter(parser, 'idxPlotLayer', [])
            parse(parser, varargin{:})
            layerNum = parser.Results.layerNum;
            gridNum = parser.Results.gridNum;
            layerDirection = parser.Results.layerDirection;
            idxPlotLayer = parser.Results.idxPlotLayer;
            if isempty(idxPlotLayer)
                idxPlotLayer = 1:layerNum;
            end
            
            % Initialise array
            sumVariable = zeros(layerNum,gridNum,gridNum);
            boundaryLayer = zeros(layerNum+1,1);
            boundaryGrid1 = zeros(gridNum+1,1);
            boundaryGrid2 = zeros(gridNum+1,1);
            midLayer = zeros(layerNum,1);
            midGrid1 = zeros(gridNum,1);
            midGrid2 = zeros(gridNum,1);
            % direction = {'x','y','z'};
            
            % Determine the indices of nodes that are not inlet or outlet
            idxInternal = ~or(obj.nodeInlet, obj.nodeOutlet);
            
            % The boundary in layerDirection
            lowerLayerBoundary = obj.boundary(layerDirection*2-1);
            upperLayerBoundary = obj.boundary(layerDirection*2);
            
            % The boundary in perpendicular directions
            perpendicularDirection = 1:3;
            perpendicularDirection = perpendicularDirection(perpendicularDirection~=layerDirection);
            lowerGrid1 = obj.boundary(perpendicularDirection(1)*2-1);
            upperGrid1 = obj.boundary(perpendicularDirection(1)*2);
            lowerGrid2 = obj.boundary(perpendicularDirection(2)*2-1);
            upperGrid2 = obj.boundary(perpendicularDirection(2)*2);
            
            % Determine the boundaries and mid
            boundaryLayer(1) = lowerLayerBoundary;
            boundaryGrid1(1) = lowerGrid1;
            boundaryGrid2(1) = lowerGrid2;
            for i = 1:layerNum
                boundaryLayer(i+1) = lowerLayerBoundary+(upperLayerBoundary-lowerLayerBoundary)*(i)/layerNum;
                midLayer(i) = mean(boundaryLayer(i:i+1));
            end
            for i = 1:gridNum
                boundaryGrid1(i+1) = lowerGrid1+(upperGrid1-lowerGrid1)*(i)/gridNum;
                boundaryGrid2(i+1) = lowerGrid2+(upperGrid2-lowerGrid2)*(i)/gridNum;
                midGrid1(i) = mean(boundaryGrid1(i:i+1));
                midGrid2(i) = mean(boundaryGrid2(i:i+1));
            end
            
            for i = idxPlotLayer
                % Determine the boundary of the layer
                boundary1 = boundaryLayer(i);
                boundary2 = boundaryLayer(i+1);
                
                % Determine the nodes in layer i
                idxLayer = and(obj.nodeCoord(idxInternal==1,layerDirection)>boundary1, obj.nodeCoord(idxInternal==1,layerDirection)<=boundary2);
                
                for j = 1:gridNum
                    % Determine the nodes in grid1 j
                    boundary1 = boundaryGrid1(j);
                    boundary2 = boundaryGrid1(j+1);
                    
                    % Determine the nodes in grid1 j
                    idxGrid1 = and(obj.nodeCoord(idxInternal==1,perpendicularDirection(1))>boundary1, obj.nodeCoord(idxInternal==1,perpendicularDirection(1))<=boundary2);
                    
                    for k = 1:gridNum
                        % Determine the nodes in grid2 k
                        boundary1 = boundaryGrid2(k);
                        boundary2 = boundaryGrid2(k+1);
                        
                        % Determine the nodes in grid1 j
                        idxGrid2 = and(obj.nodeCoord(idxInternal==1,perpendicularDirection(2))>boundary1, obj.nodeCoord(idxInternal==1,perpendicularDirection(2))<=boundary2);
                        
                        % Nodes in layer i, grid1 j, grid2 k
                        idxNodes = and(idxLayer, idxGrid1);
                        idxNodes = and(idxNodes, idxGrid2);
                        
                        % Determine the average
                        sumVariable(i,j,k) = sum(nodeVariable(idxNodes==1));
                    end
                end
            end
            
            % Fill missing values
            % averageVariable = fillmissing(averageVariable,'linear');
            % averageVariable(averageVariable<0) = 0;
            sumVariable(isnan(sumVariable)) = 0;
            
            % Plot contour plots
            for i = idxPlotLayer
                figure;
                set(gcf,'WindowState','maximized');
                contourf(midGrid1, midGrid2, reshape(sumVariable(i,:,:),gridNum,gridNum));
                colorbar;
                set(gca, 'DataAspectRatio', [1,1,1])
            end
        end
        function [] = plotSumEdgeContour(obj, edgeVariable, varargin)
            % Input parser
            parser = inputParser;
            addParameter(parser, 'layerNum', 10)
            addParameter(parser, 'gridNum', 50)
            addParameter(parser, 'layerDirection', 1)
            addParameter(parser, 'idxPlotLayer', [])
            parse(parser, varargin{:})
            layerNum = parser.Results.layerNum;
            gridNum = parser.Results.gridNum;
            layerDirection = parser.Results.layerDirection;
            idxPlotLayer = parser.Results.idxPlotLayer;
            if isempty(idxPlotLayer)
                idxPlotLayer = 1:layerNum;
            end
            
            % Initialise array
            sumVariable = zeros(layerNum,gridNum,gridNum);
            boundaryLayer = zeros(layerNum+1,1);
            boundaryGrid1 = zeros(gridNum+1,1);
            boundaryGrid2 = zeros(gridNum+1,1);
            midLayer = zeros(layerNum,1);
            midGrid1 = zeros(gridNum,1);
            midGrid2 = zeros(gridNum,1);
            % direction = {'x','y','z'};
            
            % The boundary in layerDirection
            lowerLayerBoundary = obj.boundary(layerDirection*2-1);
            upperLayerBoundary = obj.boundary(layerDirection*2);
            
            % The boundary in perpendicular directions
            perpendicularDirection = 1:3;
            perpendicularDirection = perpendicularDirection(perpendicularDirection~=layerDirection);
            lowerGrid1 = obj.boundary(perpendicularDirection(1)*2-1);
            upperGrid1 = obj.boundary(perpendicularDirection(1)*2);
            lowerGrid2 = obj.boundary(perpendicularDirection(2)*2-1);
            upperGrid2 = obj.boundary(perpendicularDirection(2)*2);
            
            % Determine the boundaries and mid
            boundaryLayer(1) = lowerLayerBoundary;
            boundaryGrid1(1) = lowerGrid1;
            boundaryGrid2(1) = lowerGrid2;
            for i = 1:layerNum
                boundaryLayer(i+1) = lowerLayerBoundary+(upperLayerBoundary-lowerLayerBoundary)*(i)/layerNum;
                midLayer(i) = mean(boundaryLayer(i:i+1));
            end
            for i = 1:gridNum
                boundaryGrid1(i+1) = lowerGrid1+(upperGrid1-lowerGrid1)*(i)/gridNum;
                boundaryGrid2(i+1) = lowerGrid2+(upperGrid2-lowerGrid2)*(i)/gridNum;
                midGrid1(i) = mean(boundaryGrid1(i:i+1));
                midGrid2(i) = mean(boundaryGrid2(i:i+1));
            end
            
            for i = idxPlotLayer
                % Determine the boundary of the layer
                boundary1 = boundaryLayer(i);
                boundary2 = boundaryLayer(i+1);
                
                % Determine the nodes in layer i
                idxLayer = and(obj.edgeCoord(:,layerDirection)>boundary1, obj.edgeCoord(:,layerDirection)<=boundary2);
                
                for j = 1:gridNum
                    % Determine the nodes in grid1 j
                    boundary1 = boundaryGrid1(j);
                    boundary2 = boundaryGrid1(j+1);
                    
                    % Determine the nodes in grid1 j
                    idxGrid1 = and(obj.edgeCoord(:,perpendicularDirection(1))>boundary1, obj.edgeCoord(:,perpendicularDirection(1))<=boundary2);
                    
                    for k = 1:gridNum
                        % Determine the nodes in grid2 k
                        boundary1 = boundaryGrid2(k);
                        boundary2 = boundaryGrid2(k+1);
                        
                        % Determine the nodes in grid1 j
                        idxGrid2 = and(obj.edgeCoord(:,perpendicularDirection(2))>boundary1, obj.edgeCoord(:,perpendicularDirection(2))<=boundary2);
                        
                        % Nodes in layer i, grid1 j, grid2 k
                        idxEdges = and(idxLayer, idxGrid1);
                        idxEdges = and(idxEdges, idxGrid2);
                        
                        % Determine the average
                        sumVariable(i,j,k) = sum(edgeVariable(idxEdges==1));
                    end
                end
            end
            
            % Fill missing values
            % averageVariable = fillmissing(averageVariable,'linear');
            % averageVariable(averageVariable<0) = 0;
            sumVariable(isnan(sumVariable)) = 0;
            
            % Plot contour plots
            for i = idxPlotLayer
                figure;
                set(gcf,'WindowState','maximized');
                contourf(midGrid1, midGrid2, reshape(sumVariable(i,:,:),gridNum,gridNum));
                colorbar;
                set(gca, 'DataAspectRatio', [1,1,1])
            end
        end
        
        % Method to plot from three directions given arrays of 
        % (1) idxEdge: Boolean for which edges to plot
        % (2) idxColour: Which colour to plot in
        function [] = plotEdges(obj, idxEdge, idxColour, varargin)
            % Input parser
            parser = inputParser;
            addParameter(parser, 'plotView', 1:3)
            parse(parser, varargin{:})
            plotView = parser.Results.plotView;
            
            % Number of edges to plot
            edgesToPlotNum = sum(idxEdge);
            idxEdge = find(idxEdge==1);
            
            % Plot properties
            colourList = Constants.colourList;
            
            % Plot in three directions
            xlabelList = {'x','y','x'};
            ylabelList = {'z (Flow direction)','z (Flow direction)','y'}; 
            plotList = [1,3; 2,3; 1,2];
            for i = plotView
                figure;
                set(gcf,'WindowState','maximized');
                for j = 1:edgesToPlotNum
                    edge = idxEdge(j);
                    plot([obj.nodeCoord(obj.edgeList(edge,1),plotList(i,1)),obj.nodeCoord(obj.edgeList(edge,2),plotList(i,1))], ...
                        [obj.nodeCoord(obj.edgeList(edge,1),plotList(i,2)),obj.nodeCoord(obj.edgeList(edge,2),plotList(i,2))], ...
                        'LineWidth',Constants.lineWidth,'Color',colourList(idxColour(edge),:));
                    hold on;
                end
                xlim([floor(min(obj.nodeCoord(:,plotList(i,1)))),ceil(max(obj.nodeCoord(:,plotList(i,1))))]);
                ylim([floor(min(obj.nodeCoord(:,plotList(i,2)))),ceil(max(obj.nodeCoord(:,plotList(i,2))))]);
                set(gca, 'DataAspectRatio', [1,1,1])
                set(gca, 'FontName','Helvetica','FontSize',36);
                xlabel(xlabelList{i},'FontName','Helvetica','FontSize',36);
                ylabel(ylabelList{i},'FontName','Helvetica','FontSize',36);
            end
        end
        
        % ===== Saving & Loading ==========================================
        % Remove unnecessary properties to reduce size
        function [obj] = removeProperties(obj)
            if strcmp(obj.networkParameters.assignmentMethod, 'actual')
                variables = {'filterParameters', 'networkParameters'};
            else
                variables = {'filterParameters', 'networkParameters', ...
                    'nodeSurfaceArea', 'nodeVolume', 'edgeArea', 'edgeDiameter'};
            end
            
            filter = Filter();
            for i = 1:length(variables)
                filter.(sprintf('%s', variables{i})) = obj.(sprintf('%s', variables{i}));
            end
            obj = filter;
        end
        
        % Re-generate required properties
        function [obj] = regenerateProperties(obj, loadMethod)
            % Load data
            filter = obj;
            filter = filter.loadDTM(loadMethod);
            filter = filter.generateNetwork();
            
            % Use defined properties
            if strcmp(obj.networkParameters.assignmentMethod, 'random')
                variables = {'filterParameters', 'networkParameters', ...
                    'nodeSurfaceArea', 'nodeVolume', 'edgeArea', 'edgeDiameter'};
                for i = 1:length(variables)
                    filter.(sprintf('%s', variables{i})) = obj.(sprintf('%s', variables{i}));
                end
            end
            
            % Output
            obj = filter;
        end
        
        % Load a state of object
        function [obj] = loadState(obj, saveLocation, varargin)
            % Input parser
            parser = inputParser;
            addParameter(parser, 'variables', 'all')
            parse(parser, varargin{:})
            variables = parser.Results.variables;
            
            % Load a state of object
            if all(strcmp(variables, 'all'))
                filter = load(saveLocation);
                obj = filter.filter;
            else
                filter = load(saveLocation);
                filter = filter.filter;
                for i = 1:length(variables)
                    obj.(sprintf('%s',variables{i})) = filter.(sprintf('%s',variables{i}));
                end
            end
        end
        
        % Store the current state of object
        function [] = saveState(obj, saveLocation, varargin)
            % Input parser
            parser = inputParser;
            addParameter(parser, 'variables', 'all')
            parse(parser, varargin{:})
            variables = parser.Results.variables;
            
            % Save current state of object
            if all(strcmp(variables, 'all'))
                filter = obj;
                save(saveLocation, 'filter');
            else
                filter = Filter();
                for i = 1:length(variables)
                    filter.(sprintf('%s',variables{i})) = obj.(sprintf('%s',variables{i}));
                end
                save(saveLocation, 'filter');
            end
        end
    end
    
    methods (Access = private)
        % ===== Pre-processing ============================================
        % Remove nodes that are not connected to the overall network
        function [obj] = removeUnconnectedNodes(obj)
            % Determine the nodes that are not connected to the system
            idxNodes = sort(unique(obj.edgeList(:)));
            idxNotConnected = find(~ismember(1:obj.nodeNum,idxNodes));
            
            % Remove these nodes
            obj.nodeCoord = obj.nodeCoord(idxNodes,:);
            obj.nodeVolume = obj.nodeVolume(idxNodes,:);
            obj.nodeBoundary = obj.nodeBoundary(idxNodes,:);
            obj.nodeSurfaceArea = obj.nodeSurfaceArea(idxNodes,:);
            obj.nodeNum = size(obj.nodeCoord,1);
            
            % Update the indices in edgeList
            for i = 1:length(idxNotConnected)
                idxNode = idxNotConnected(i);
                obj.edgeList(obj.edgeList(:,1)>idxNode,1) = obj.edgeList(obj.edgeList(:,1)>idxNode,1)-1;
                obj.edgeList(obj.edgeList(:,2)>idxNode,2) = obj.edgeList(obj.edgeList(:,2)>idxNode,2)-1;
                idxNotConnected(idxNotConnected>idxNode) = idxNotConnected(idxNotConnected>idxNode)-1;
            end
        end
        
        % Normalise geometric properties
        function [obj] = normaliseGeometricProperties(obj, factor)
            % Normalise geometric properties
            obj.edgeArea = obj.edgeArea*(factor^2);
            obj.edgeCoord = obj.edgeCoord*factor;
            obj.edgeDiameter = obj.edgeDiameter*factor;
            obj.edgeLength = obj.edgeLength*factor;
            obj.nodeCoord = obj.nodeCoord*factor;
            obj.nodeSurfaceArea = obj.nodeSurfaceArea*(factor^2);
            obj.nodeVolume = obj.nodeVolume*(factor^3);
            obj.particleCoord = obj.particleCoord*factor;
            obj.particleRadius = obj.particleRadius*factor;
            obj.boundary = obj.boundary*factor;
        end
        
        % Recalculate edgeArea (remove negative values)
        function [obj] = recalculateEdgeArea(obj)
            for i = 1:obj.edgeNum
                % Determine the coordinates of the three particles forming the
                % constriction
                a = obj.particleCoord(obj.edgeParticles(i,1),:);
                b = obj.particleCoord(obj.edgeParticles(i,2),:);
                c = obj.particleCoord(obj.edgeParticles(i,3),:);
                
                % Determine the radii
                r = zeros(3,1);
                r(1) = obj.particleRadius(obj.edgeParticles(i,1));
                r(2) = obj.particleRadius(obj.edgeParticles(i,2));
                r(3) = obj.particleRadius(obj.edgeParticles(i,3));
                
                % Determine the area of the triangle
                A = [a(1), b(1), c(1); a(2), b(2), c(2); 1, 1, 1];
                B = [a(2), b(2), c(2); a(3), b(3), c(3); 1, 1, 1];
                C = [a(3), b(3), c(3); a(1), b(1), c(1); 1, 1, 1];
                areaTri = 0.5*sqrt(det(A)^2+det(B)^2+det(C)^2);
                
                % Determine the angle for each particle
                t1 = acos(dot(b-a,c-a)/norm(b-a)/norm(c-a));
                t2 = acos(dot(a-b,c-b)/norm(a-b)/norm(c-b));
                t3 = acos(dot(a-c,b-c)/norm(a-c)/norm(b-c));
                
                % Determine the area of particles
                areaPar = 0.5*t1*r(1)^2 + 0.5*t2*r(2)^2 + 0.5*t3*r(3)^2;
                
                % Determine the normal distance of each vertex to the opposite side
                coords = [a; b; c];
                d = zeros(3,1);
                areaExclude = 0;
                for j = 1:3
                    % Determine the index of other particles
                    idxOpposite = find(~ismember(1:3,j));
                    
                    % Determine the unit normal to the constriction face
                    v1 = coords(1,:)-coords(2,:);
                    v2 = coords(1,:)-coords(3,:);
                    n = cross(v1, v2);
                    n = n/norm(n);
                    
                    % Determine the line perpendicular to opposite side passing through
                    % particle j
                    v1 = coords(idxOpposite(1),:)-coords(idxOpposite(2),:);
                    v3 = cross(v1, n);
                    v3 = v3/norm(v3);
                    
                    % Determine the coordinates of intersecting point
                    v4 = coords(j,:);
                    coef = (dot(coords(idxOpposite(1),:),v3)-dot(v4,v3))/dot(v3,v3);
                    intersect = v4+coef*v3;
                    
                    % The distance
                    d(j) = abs(coef);
                    
                    % Calculate area if radius is larger than normal distance
                    if and(r(j)>d(j), norm(v1)>norm(intersect-coords(idxOpposite(1),:)))
                        theta = acos(d(j)/r(j))*2;
                        areaExclude = areaExclude+0.5*theta*r(j)^2-0.5*r(j)^2*sin(theta);
                    end
                end
                
                % edgeArea
                obj.edgeArea(i) = areaTri-areaPar+areaExclude;
            end
        end
        
        % Calculate edgeDiameter and edgeReduced
        function [obj] = calculateEdgeDiameter(obj, originalRadius, moreParticles, smallestConstriction, visualiser, idxEdges)
            % The script determines the constriction size distribution of random
            % sphere packings, in which the particles forming a constriction are
            % determined through Delauney triangulation method (both classical and
            % modified)
            
            % References:
            % (1) Shire, T., 2014. Micro-scale modelling of granular filters (Doctoral
            % dissertation, Imperial College London).
            % (2) Reboul, N., Vincens, E. and Cambou, B., 2010. A computational
            % procedure to assess the distribution of constriction sizes for an
            % assembly of spheres. Computers and Geotechnics, 37(1-2), pp.195-206.
            
            % idxEdges must be a row vector
            if size(idxEdges,2)~=length(idxEdges)
                idxEdges = idxEdges';
            end
            
            if ~visualiser
                idxEdges = 1:obj.edgeNum;
            end
            
            % Determine the maximum particle size
            maxParticleRadius = max(obj.particleRadius);
            
            % Determine if each edge is connected to a node at boundary
            if moreParticles
                edgeBoundary = ones(obj.edgeNum,1);
            else
                idxIsBoundary = sum(obj.nodeBoundary,2);
                idxIsBoundary = find(idxIsBoundary>0);
                edgeBoundary = zeros(obj.edgeNum,1);
                for i = 1:length(idxIsBoundary)
                    idxNode = idxIsBoundary(i);
                    edgeBoundary = or(edgeBoundary, obj.edgeList(:,1)==idxNode);
                    edgeBoundary = or(edgeBoundary, obj.edgeList(:,2)==idxNode);
                end
            end
            
            % Initialise matrix to store centroids, radii and normals of constrictions
            edgeCentroid = zeros(obj.edgeNum,3);
            edgeRadius = zeros(obj.edgeNum,1);
            edgeDirection = zeros(obj.edgeNum,3);
            edgeOverlapNum = zeros(obj.edgeNum,1);
            % obj.edgeOverlapNum = zeros(obj.edgeNum,1);
            edgeSamePlaneNum = zeros(obj.edgeNum,1);
            % obj.edgeOverlap = zeros(obj.edgeNum,1);
            
            % For each constriction
            fprintf('Calculating edgeDiameter ');
            for i = idxEdges
                % Display progress
                if i~=1
                    fprintf(repmat('\b', 1, 7));
                end
                fprintf('%6.2f%%', i/obj.edgeNum*100);
                
                % Determine the centroid and radius by solving simultaneous equation
                [x] = inscribedCircle([obj.particleCoord(obj.edgeParticles(i,1),:),obj.particleRadius(obj.edgeParticles(i,1))], ...
                    [obj.particleCoord(obj.edgeParticles(i,2),:),obj.particleRadius(obj.edgeParticles(i,2))], ...
                    [obj.particleCoord(obj.edgeParticles(i,3),:),obj.particleRadius(obj.edgeParticles(i,3))]);
                edgeCentroid(i,:) = x(1:3);
                edgeRadius(i) = x(4);
                
                % Function for the constriction plane: [x] = [v1]+a[v2]+b[v3]
                v1 = obj.particleCoord(obj.edgeParticles(i,1),:);
                v2 = obj.particleCoord(obj.edgeParticles(i,1),:)-obj.particleCoord(obj.edgeParticles(i,2),:);
                v3 = obj.particleCoord(obj.edgeParticles(i,1),:)-obj.particleCoord(obj.edgeParticles(i,3),:);
                
                % Determine the normal to the plane
                n = cross(v2,v3);
                edgeDirection(i,:) = n/norm(n);
                
                if visualiser
                    [xs,ys,zs] = sphere;
                    figure;
                    set(gcf, 'WindowState', 'maximized')
                    for j = 1:3
                        p = obj.particleCoord(obj.edgeParticles(i,j),:);
                        q = obj.particleRadius(obj.edgeParticles(i,j));
                        surf(xs*q+p(1), ys*q+p(2), zs*q+p(3));
                        hold on;
                    end
                    p = edgeCentroid(i,:);
                    q = edgeRadius(i);
                    surf(xs*q+p(1), ys*q+p(2), zs*q+p(3));
                    set(gca, 'DataAspectRatio', [1,1,1])
                end
                
                % Determine a list of neighbouring particles
                if edgeBoundary(i)==0
                    neighbourList = zeros(obj.edgeNum,1);
                    for j = 1:3
                        for k = 1:3
                            neighbourList = or(neighbourList==1,obj.edgeParticles(:,k)==obj.edgeParticles(i,j));
                        end
                    end
                    neighbourList = obj.edgeParticles(neighbourList==1,:);
                    neighbourList = neighbourList(:);
                    neighbourList = unique(neighbourList);
                    neighbourList = neighbourList(neighbourList~=obj.edgeParticles(i,1));
                    neighbourList = neighbourList(neighbourList~=obj.edgeParticles(i,2));
                    neighbourList = neighbourList(neighbourList~=obj.edgeParticles(i,3));
                else
                    a = 2;
                    minX = edgeCentroid(i,1)-a*maxParticleRadius-edgeRadius(i);
                    maxX = edgeCentroid(i,1)+a*maxParticleRadius+edgeRadius(i);
                    minY = edgeCentroid(i,2)-a*maxParticleRadius-edgeRadius(i);
                    maxY = edgeCentroid(i,2)+a*maxParticleRadius+edgeRadius(i);
                    minZ = edgeCentroid(i,3)-a*maxParticleRadius-edgeRadius(i);
                    maxZ = edgeCentroid(i,3)+a*maxParticleRadius+edgeRadius(i);
                    neighbourList = ones(obj.particleNum,1);
                    neighbourList = and(neighbourList, obj.particleCoord(:,1)>minX);
                    neighbourList = and(neighbourList, obj.particleCoord(:,1)<maxX);
                    neighbourList = and(neighbourList, obj.particleCoord(:,2)>minY);
                    neighbourList = and(neighbourList, obj.particleCoord(:,2)<maxY);
                    neighbourList = and(neighbourList, obj.particleCoord(:,3)>minZ);
                    neighbourList = and(neighbourList, obj.particleCoord(:,3)<maxZ);
                    neighbourList = find(neighbourList==1);
                    neighbourList = neighbourList(neighbourList~=obj.edgeParticles(i,1));
                    neighbourList = neighbourList(neighbourList~=obj.edgeParticles(i,2));
                    neighbourList = neighbourList(neighbourList~=obj.edgeParticles(i,3));
                end
                neighbourNum = length(neighbourList);
                
                % Initialise empty arrays
                neighbourOverlap = zeros(neighbourNum,1);
                neighbourData = zeros(neighbourNum,4);
                neighbourSamePlane = zeros(neighbourNum,1);
                
                % For each neighbouring particle
                for j = 1:neighbourNum
                    neighbour = neighbourList(j);
                    
                    % Check if the particle overlaps with the inscribed circle, if yes,
                    % store the centre and radius of the circle
                    % intersecting the constriction plane
                    
                    % Determine the centroid of the projected image on the plane
                    % Simultaneous equations
                    % (1) Equation of coordinates on plane [x].[n] = [v1].[n]
                    % (2) Normal line passing through particle [x] = [v4]+c[n]
                    v4 = obj.particleCoord(neighbour,1:3);
                    n = edgeDirection(i,:);
                    coef = (dot(v1,n)-dot(v4,n))/dot(n,n);
                    neighbourData(j,1:3) = v4+coef*n;
                    
                    % Conditions for overlapping are:
                    % (1) Particle intersects with the plane parallel to constriction
                    %     face
                    % (2) Intersecting circle overlaps with the inscribed circle
                    distanceToPlane = norm(dot(obj.particleCoord(neighbour,:)-obj.particleCoord(obj.edgeParticles(i,1),:), edgeDirection(i,:)));
                    neighbourRadius = obj.particleRadius(neighbour);
                    if neighbourRadius>distanceToPlane
                        neighbourSamePlane(j) = 1;
                        if originalRadius
                            neighbourData(j,4) = neighbourRadius;
                        else
                            neighbourData(j,4) = sqrt(obj.particleRadius(neighbour)^2-distanceToPlane^2);
                        end
                        distanceToCircle = sqrt(sum((neighbourData(j,1:3)-edgeCentroid(i,:)).^2));
                        if distanceToCircle<(neighbourData(j,4)+edgeRadius(i))
                            neighbourOverlap(j) = 1;
                        end
                    end
                end
                edgeOverlapNum(i) = sum(neighbourOverlap);
                edgeSamePlaneNum(i) = sum(neighbourSamePlane);
                
                % If there is overlap
                if edgeOverlapNum(i)~=0
                    % % Determine the corresponding edge
                    % idxNeighbour = neighbourList(neighbourOverlap==1);
                    % for j = 1:obj.edgeOverlapNum(i)
                    %     neighbour = idxNeighbour(j);
                    %     combList = combnk(1:3,2);
                    %     for k = 1:size(combList,1)
                    %         % Edges with the overlapping particle as one of
                    %         % the three particles
                    %         idxCorresponding = obj.edgeParticles(:,1)==neighbour;
                    %         idxCorresponding = or(idxCorresponding, obj.edgeParticles(:,2)==neighbour);
                    %         idxCorresponding = or(idxCorresponding, obj.edgeParticles(:,3)==neighbour);
                    %
                    %         % Edges formed by the overlapping particle plus
                    %         % any two of the three particles forming edge i
                    %         for l = 1:2
                    %             idxCorresponding2 = obj.edgeParticles(:,1)==obj.edgeParticles(i,combList(k,l));
                    %             idxCorresponding2 = or(idxCorresponding2, obj.edgeParticles(:,2)==obj.edgeParticles(i,combList(k,l)));
                    %             idxCorresponding2 = or(idxCorresponding2, obj.edgeParticles(:,3)==obj.edgeParticles(i,combList(k,l)));
                    %             idxCorresponding = and(idxCorresponding, idxCorresponding2);
                    %         end
                    %
                    %         % Determine the index of edge
                    %         if any(idxCorresponding==1)
                    %             % Corresponding edge must be from another
                    %             % node
                    %             if obj.edgeOverlap(i)~=0
                    %                 i
                    %                 neighbour
                    %             end
                    %             obj.edgeOverlap(i) = find(idxCorresponding==1);
                    %         end
                    %     end
                    % end
                    
                    % Determine the list of combinations of three particles
                    if moreParticles
                        neighbourData = [neighbourData(neighbourSamePlane==1,:);[obj.particleCoord(obj.edgeParticles(i,:)',:),obj.particleRadius(obj.edgeParticles(i,:)')]];
                        neighbourOverlap = neighbourOverlap(neighbourSamePlane==1);
                    else
                        neighbourData = [neighbourData(neighbourOverlap==1,:);[obj.particleCoord(obj.edgeParticles(i,:)',:),obj.particleRadius(obj.edgeParticles(i,:)')]];
                        neighbourOverlap = neighbourOverlap(neighbourOverlap==1);
                    end
                    n = length(neighbourOverlap);
                    idxOverlap = [find(neighbourOverlap==1);n+1;n+2;n+3];
                    combList = combnk(idxOverlap,3);
                    % combList = zeros(n*3,3);
                    % for j = 1:n
                    %     combList(j*3-2,:) = [j,n+1,n+2];
                    %     combList(j*3-1,:) = [j,n+1,n+3];
                    %     combList(j*3,:) = [j,n+2,n+3];
                    % end
                    
                    if visualiser
                        figure;
                        set(gcf, 'WindowState', 'maximized')
                        for j = 1:size(neighbourData,1)
                            p = neighbourData(j,:);
                            surf(xs*p(4)+p(1), ys*p(4)+p(2), zs*p(4)+p(3));
                            hold on;
                        end
                        p = [edgeCentroid(i,:),edgeRadius(i)];
                        surf(xs*p(4)+p(1), ys*p(4)+p(2), zs*p(4)+p(3));
                        set(gca, 'DataAspectRatio', [1,1,1])
                    end
                    
                    % For each combination of three particles
                    if ~smallestConstriction
                        edgeRadius(i) = 0;
                    end
                    for j = 1:size(combList,1)
                        idx = combList(j,:);
                        idxNeighbour = find(~ismember(1:n+3,combList(j,:)));
                        
                        % Determine the centroid and radius
                        [x] = inscribedCircle(neighbourData(idx(1),:), neighbourData(idx(2),:), neighbourData(idx(3),:));
                        edgeCentroidCurrent = x(1:3);
                        edgeRadiusCurrent = x(4);
                        
                        % Initialise array
                        neighbourOverlap = zeros(size(idxNeighbour));
                        
                        % For each particle
                        for k = 1:length(idxNeighbour)
                            % Check if any particle overlaps with the inscribed circle
                            
                            % Condition for overlapping is:
                            % (1) Any particle intersect with the inscribed circle
                            neighbourRadius = neighbourData(idxNeighbour(k),end);
                            distanceToCircle = sqrt(sum((neighbourData(idxNeighbour(k),1:3)-edgeCentroidCurrent).^2));
                            if distanceToCircle<(neighbourRadius+edgeRadiusCurrent)
                                neighbourOverlap(k) = 1;
                            end
                        end
                        
                        % If there is no overlap, store the data
                        if ~sum(neighbourOverlap)
                            if smallestConstriction
                                if edgeRadiusCurrent<edgeRadius(i)
                                    edgeCentroid(i,:) = edgeCentroidCurrent;
                                    edgeRadius(i) = edgeRadiusCurrent;
                                end
                            else
                                if edgeRadiusCurrent>edgeRadius(i)
                                    edgeCentroid(i,:) = edgeCentroidCurrent;
                                    edgeRadius(i) = edgeRadiusCurrent;
                                end
                            end
                        end
                    end
                    
                    if visualiser
                        figure;
                        set(gcf, 'WindowState', 'maximized')
                        for j = 1:size(neighbourData,1)
                            p = neighbourData(j,:);
                            surf(xs*p(4)+p(1), ys*p(4)+p(2), zs*p(4)+p(3));
                            hold on;
                        end
                        p = [edgeCentroid(i,:),edgeRadius(i)];
                        surf(xs*p(4)+p(1), ys*p(4)+p(2), zs*p(4)+p(3));
                        set(gca, 'DataAspectRatio', [1,1,1])
                    end
                end
            end
            fprintf('\n');
            
            % Store data
            obj.edgeCoord = edgeCentroid;
            obj.edgeDiameter = edgeRadius*2;
            obj.edgeNormal = edgeDirection;
        end
        function [obj] = calculateEdgeReduced(obj, visualiser, idxEdges)
            % idxEdges must be a row vector
            if size(idxEdges,2)~=length(idxEdges)
                idxEdges = idxEdges';
            end
            if ~visualiser
                idxEdges = 1:obj.edgeNum;
            end
            
            % Colourlist
            colourList = Constants.colourList;
            
            % edgeRadius
            edgeRadius = obj.edgeDiameter/2;
            
            % Initialise arrays
            obj.edgeReduced = zeros(obj.edgeNum,1);
            edgeIsInside = zeros(obj.edgeNum,1);
            edgeCircleArea = zeros(obj.edgeNum,1);
            
            % For each constriction
            for i = idxEdges
                % The coordinates of the triangle
                P = zeros(3);
                P(1,:) = obj.particleCoord(obj.edgeParticles(i,1),:);
                P(2,:) = obj.particleCoord(obj.edgeParticles(i,2),:);
                P(3,:) = obj.particleCoord(obj.edgeParticles(i,3),:);
                
                % The coordinates of the centroid
                circleCentroid = obj.edgeCoord(i,:);
                
                if visualiser
                    % Sphere plot
                    [xs,ys,zs] = sphere;
                    figure;
                    set(gcf, 'WindowState', 'maximized')
                    for j = 1:3
                        p = obj.particleCoord(obj.edgeParticles(i,j),:);
                        q = obj.particleRadius(obj.edgeParticles(i,j));
                        surf(xs*q+p(1), ys*q+p(2), zs*q+p(3));
                        hold on;
                    end
                    p = obj.edgeCoord(i,:);
                    q = edgeRadius(i);
                    surf(xs*q+p(1), ys*q+p(2), zs*q+p(3));
                    set(gca, 'DataAspectRatio', [1,1,1])
                    
                    % Show triangle
                    figure;
                    set(gcf, 'WindowState', 'maximized')
                    plot3(P(:,1),P(:,2),P(:,3),'Color',colourList(1,:),'LineWidth',3);
                    hold on;
                    plot3([P(3,1),P(1,1)],[P(3,2),P(1,2)],[P(3,3),P(1,3)],'Color',colourList(1,:),'LineWidth',3);
                    scatter3(circleCentroid(1),circleCentroid(2),circleCentroid(3),'LineWidth',3);
                    set(gca, 'DataAspectRatio', [1,1,1])
                end
                
                % Determine if the edgeCentroid is within the triangle
                isInside = zeros(3,1);
                combList = [1, 2, 3; 2, 3, 1; 3, 1, 2];
                for j = 1:3
                    % Determine the normal with third point of triangle
                    v1 = P(combList(j,2),:)-P(combList(j,1),:);
                    v2 = P(combList(j,3),:)-P(combList(j,1),:);
                    n1 = cross(v1, v2);
                    
                    % Determine the normal with the circle centroid
                    v3 = circleCentroid-P(combList(j,1),:);
                    n2 = cross(v1, v3);
                    
                    % Determine dot product
                    direction = dot(n1, n2);
                    
                    if direction>0
                        isInside(j) = 1;
                    end
                end
                edgeIsInside(i) = all(isInside);
                
                % Determine the area of insribed circle that coincides with the
                % triangle
                if edgeIsInside(i)
                    % Determine the area of circle
                    edgeCircleArea(i) = pi*(edgeRadius(i)^2);
                    
                    % Determine the normal distance of insribed circle to the opposite side
                    d = zeros(3,1);
                    areaExclude = 0;
                    for j = 1:3
                        % Determine the index of other particles
                        idxOpposite = find(~ismember(1:3,j));
                        
                        % Determine the unit normal to the constriction face
                        v1 = P(1,:)-P(2,:);
                        v2 = P(1,:)-P(3,:);
                        n = cross(v1, v2);
                        n = n/norm(n);
                        
                        % Determine the line perpendicular to opposite side passing through
                        % particle j
                        v1 = P(idxOpposite(1),:)-P(idxOpposite(2),:);
                        v3 = cross(v1, n);
                        v3 = v3/norm(v3);
                        
                        % Determine the coordinates of intersecting point
                        v4 = circleCentroid;
                        coef = (dot(P(idxOpposite(1),:),v3)-dot(v4,v3))/dot(v3,v3);
                        intersect = v4+coef*v3;
                        if visualiser
                            scatter3(intersect(1),intersect(2),intersect(3),'LineWidth',3);
                        end
                        
                        % The distance
                        d(j) = abs(coef);
                        
                        % Calculate area if radius is larger than normal distance
                        if and(edgeRadius(i)>d(j), norm(v1)>norm(intersect-P(idxOpposite(1),:)))
                            theta = acos(d(j)/edgeRadius(i))*2;
                            areaExclude = areaExclude+0.5*theta*edgeRadius(i)^2-0.5*edgeRadius(i)^2*sin(theta);
                        end
                    end
                    
                    % The area of circle in triangle
                    edgeCircleArea(i) = edgeCircleArea(i)-areaExclude;
                else
                    % Determine which side the centroid of circle is on
                    if length(find(isInside==0))~=1
                        fprintf('Edge %i\n',i);
                    end
                    assert(length(find(isInside==0))==1);
                    idxOpposite = combList(isInside==0,1:2);
                    
                    % Determine the unit normal to the constriction face
                    v1 = P(1,:)-P(2,:);
                    v2 = P(1,:)-P(3,:);
                    n = cross(v1, v2);
                    n = n/norm(n);
                    
                    % Determine the line perpendicular to opposite side passing through
                    % particle j
                    v1 = P(idxOpposite(1),:)-P(idxOpposite(2),:);
                    v3 = cross(v1, n);
                    v3 = v3/norm(v3);
                    
                    % Determine the coordinates of intersecting point
                    v4 = circleCentroid;
                    coef = (dot(P(idxOpposite(1),:),v3)-dot(v4,v3))/dot(v3,v3);
                    intersect = v4+coef*v3;
                    if visualiser
                        scatter3(intersect(1),intersect(2),intersect(3),'LineWidth',3);
                        v5 = v4+sign(coef)*edgeRadius(i)*v3;
                        scatter3(v5(1),v5(2),v5(3),'LineWidth',3);
                    end
                    
                    % The distance
                    d = abs(coef);
                    
                    % Calculate area
                    if and(edgeRadius(i)>d, norm(v1)>norm(intersect-P(idxOpposite(1),:)))
                        theta = acos(d/edgeRadius(i))*2;
                        edgeCircleArea(i) = 0.5*theta*edgeRadius(i)^2-0.5*edgeRadius(i)^2*sin(theta);
                    end
                end
                
                % Determine the ratio of reduced area to edgeArea
                obj.edgeReduced(i) = (obj.edgeArea(i)-edgeCircleArea(i))/obj.edgeArea(i);
            end
            assert(all(obj.edgeReduced(idxEdges)>0));
        end
        
        % Segment smaller filter samples from the larger ones
        function [obj] = cutFilter(obj, loadMethod)
            % Boundary in z direction
            zlo = obj.boundary(5);
            zhi = obj.boundary(6);
            
            % Determine new zhi
            zhi = zlo+(zhi-zlo)/4*obj.filterParameters.size(3);
            obj.boundary(6) = zhi;
            
            if and(strcmp(loadMethod, 'read'), isfile(FileLocation.filterCut(obj)))
                % Load data if available
                obj = obj.loadState(FileLocation.filterCut(obj), ...
                    'variables', {'idxNodesInside','idxEdgesOutside','nodeBoundary','edgeList'});
            else
                % Determine the nodes within of the new boundary
                obj.idxNodesInside = obj.nodeCoord(:,3)<zhi;
                
                % Determine the edges connected to nodes outside the boundary
                idxNodesOutside = find(obj.idxNodesInside==0)';
                firstNodeOutside = zeros(obj.edgeNum,1);
                secondNodeOutside = zeros(obj.edgeNum,1);
                fprintf('Determining edges outside ');
                for i = 1:obj.edgeNum
                    firstNodeOutside(i) = any(obj.edgeList(i,1)==idxNodesOutside);
                    secondNodeOutside(i) = any(obj.edgeList(i,2)==idxNodesOutside);
                    
                    if i~=1
                        fprintf(repmat('\b', 1, 7));
                    end
                    fprintf('%6.2f%%', i/obj.edgeNum*100);
                end
                fprintf('\n');
                obj.idxEdgesOutside = or(firstNodeOutside, secondNodeOutside);
                
                % Among the remaining nodes, determine the number of edges
                % removed; add to nodeBoundary
                idxEdgesJustOutside = and(firstNodeOutside, ~secondNodeOutside);
                idxNodesBoundary = obj.edgeList(idxEdgesJustOutside==1,2);
                for i = 1:length(idxNodesBoundary)
                    % Current node
                    node = idxNodesBoundary(i);
                    
                    % Add to nodeBoundary
                    obj.nodeBoundary(node,obj.nodeBoundary(node,:)~=0) = obj.nodeBoundary(node,obj.nodeBoundary(node,:)~=0)+1;
                    obj.nodeBoundary(node,6) = max([max(obj.nodeBoundary(node,:)),1]);
                end
                idxEdgesJustOutside = and(~firstNodeOutside, secondNodeOutside);
                idxNodesBoundary = obj.edgeList(idxEdgesJustOutside==1,1);
                for i = 1:length(idxNodesBoundary)
                    % Current node
                    node = idxNodesBoundary(i);
                    
                    % Add to nodeBoundary
                    obj.nodeBoundary(node,obj.nodeBoundary(node,:)~=0) = obj.nodeBoundary(node,obj.nodeBoundary(node,:)~=0)+1;
                    obj.nodeBoundary(node,6) = max([max(obj.nodeBoundary(node,:)),1]);
                end
                
                % Update the indices in edgeList
                obj.edgeList = obj.edgeList(obj.idxEdgesOutside==0,:);
                fprintf('Updating edgeList ');
                isFirst = 1;
                nodeOutsideNum = length(idxNodesOutside);
                while ~isempty(idxNodesOutside)
                    % Current node
                    node = idxNodesOutside(1);
                    i = 1;
                    isConsecutive = 1;
                    while isConsecutive
                        if length(idxNodesOutside)>i
                            i = i+1;
                            nextNode = idxNodesOutside(i);
                            if nextNode-node~=1
                                isConsecutive = 0;
                                deductBy = i-1;
                            else
                                node = nextNode;
                            end
                        else
                            isConsecutive = 0;
                            deductBy = i;
                        end
                    end
                    
                    obj.edgeList(obj.edgeList(:,1)>node,1) = obj.edgeList(obj.edgeList(:,1)>node,1)-deductBy;
                    obj.edgeList(obj.edgeList(:,2)>node,2) = obj.edgeList(obj.edgeList(:,2)>node,2)-deductBy;
                    
                    % Remove entry
                    idxNodesOutside = idxNodesOutside(idxNodesOutside>node);
                    if ~isempty(idxNodesOutside)
                        idxNodesOutside = idxNodesOutside-deductBy;
                    end
                    
                    % Print progress
                    if ~isFirst
                        fprintf(repmat('\b', 1, 7));
                    else
                        isFirst = 0;
                    end
                    fprintf('%6.2f%%', (1-length(idxNodesOutside)/nodeOutsideNum)*100);
                end
                fprintf('\n');
                
                % Save data
                obj.saveState(FileLocation.filterCut(obj), ...
                    'variables', {'idxNodesInside','idxEdgesOutside','nodeBoundary','edgeList'});
            end
            
            % Remove the nodes outside the new boundary
            obj.nodeBoundary = obj.nodeBoundary(obj.idxNodesInside==1,:);
            obj.nodeCoord = obj.nodeCoord(obj.idxNodesInside==1,:);
            obj.nodeSurfaceArea = obj.nodeSurfaceArea(obj.idxNodesInside==1);
            obj.nodeVolume = obj.nodeVolume(obj.idxNodesInside==1);
            
            % Remove edges outside the boundary
            obj.edgeArea = obj.edgeArea(obj.idxEdgesOutside==0);
            obj.edgeCoord = obj.edgeCoord(obj.idxEdgesOutside==0,:);
            obj.edgeDiameter = obj.edgeDiameter(obj.idxEdgesOutside==0);
            obj.edgeLength = obj.edgeLength(obj.idxEdgesOutside==0);
            obj.edgeNormal = obj.edgeNormal(obj.idxEdgesOutside==0,:);
            obj.edgeParticles = obj.edgeParticles(obj.idxEdgesOutside==0,:);
            obj.edgeReduced = obj.edgeReduced(obj.idxEdgesOutside==0);
            
            % Update nodeNum and edgeNum
            obj.nodeNum = size(obj.nodeBoundary,1);
            obj.edgeNum = size(obj.edgeArea,1);
        end
        function [obj] = segmentFilter(obj, loadMethod)
            % Boundary in z direction
            zlo = obj.boundary(5);
            zhi = obj.boundary(6);
            
            % Determine new zlo and zhi
            zloNew = zlo+(zhi-zlo)/4*0.5;
            switch obj.filterParameters.size(3)
                case 1
                    zhiNew = zlo+(zhi-zlo)/4*1.5;
                case 2
                    zhiNew = zlo+(zhi-zlo)/4*2.5;
                case 3
                    zhiNew = zlo+(zhi-zlo)/4*3.5;
            end
            obj.boundary(5) = zloNew;
            obj.boundary(6) = zhiNew;
            
            if and(strcmp(loadMethod, 'read'), isfile(FileLocation.filterCut(obj)))
                % Load data if available
                obj = obj.loadState(FileLocation.filterCut(obj), ...
                    'variables', {'idxNodesInside','idxEdgesOutside','nodeBoundary','edgeList'});
            else
                % Determine the nodes within of the new boundary
                idxNodesInsideLo = obj.nodeCoord(:,3)>zloNew;
                idxNodesInsideHi = obj.nodeCoord(:,3)<zhiNew;
                obj.idxNodesInside = and(idxNodesInsideLo, idxNodesInsideHi);
                
                % Determine the edges connected to nodes outside the boundary
                idxNodesOutside = find(obj.idxNodesInside==0)';
                firstNodeOutside = zeros(obj.edgeNum,1);
                secondNodeOutside = zeros(obj.edgeNum,1);
                fprintf('Determining edges outside ');
                for i = 1:obj.edgeNum
                    firstNodeOutside(i) = any(obj.edgeList(i,1)==idxNodesOutside);
                    secondNodeOutside(i) = any(obj.edgeList(i,2)==idxNodesOutside);
                    
                    if i~=1
                        fprintf(repmat('\b', 1, 7));
                    end
                    fprintf('%6.2f%%', i/obj.edgeNum*100);
                end
                fprintf('\n');
                obj.idxEdgesOutside = or(firstNodeOutside, secondNodeOutside);
                
                % Clear all data for zlo and zhi in nodeBoundary
                obj.nodeBoundary(:,5:6) = 0;
                
                % Among the remaining nodes, determine the number of edges
                % removed; add to nodeBoundary
                zmid = mean([zloNew,zhiNew]);
                idxEdgesJustOutside = and(firstNodeOutside, ~secondNodeOutside);
                idxNodesBoundary = obj.edgeList(idxEdgesJustOutside==1,2);
                for i = 1:length(idxNodesBoundary)
                    % Current node
                    node = idxNodesBoundary(i);
                    
                    % Determine if it's an inlet or outlet
                    if obj.nodeCoord(node,3)<zmid
                        % zlo
                        idxBoundary = 5;
                    else
                        % zhi
                        idxBoundary = 6;
                    end
                    
                    % Add to nodeBoundary
                    obj.nodeBoundary(node,obj.nodeBoundary(node,:)~=0) = obj.nodeBoundary(node,obj.nodeBoundary(node,:)~=0)+1;
                    obj.nodeBoundary(node,idxBoundary) = max([max(obj.nodeBoundary(node,:)),1]);
                end
                idxEdgesJustOutside = and(~firstNodeOutside, secondNodeOutside);
                idxNodesBoundary = obj.edgeList(idxEdgesJustOutside==1,1);
                for i = 1:length(idxNodesBoundary)
                    % Current node
                    node = idxNodesBoundary(i);
                    
                    % Determine if it's an inlet or outlet
                    if obj.nodeCoord(node,3)<zmid
                        % zlo
                        idxBoundary = 5;
                    else
                        % zhi
                        idxBoundary = 6;
                    end
                    
                    % Add to nodeBoundary
                    obj.nodeBoundary(node,obj.nodeBoundary(node,:)~=0) = obj.nodeBoundary(node,obj.nodeBoundary(node,:)~=0)+1;
                    obj.nodeBoundary(node,idxBoundary) = max([max(obj.nodeBoundary(node,:)),1]);
                end
                
                % Update the indices in edgeList
                obj.edgeList = obj.edgeList(obj.idxEdgesOutside==0,:);
                fprintf('Updating edgeList ');
                isFirst = 1;
                nodeOutsideNum = length(idxNodesOutside);
                while ~isempty(idxNodesOutside)
                    % Current node
                    node = idxNodesOutside(1);
                    i = 1;
                    isConsecutive = 1;
                    while isConsecutive
                        if length(idxNodesOutside)>i
                            i = i+1;
                            nextNode = idxNodesOutside(i);
                            if nextNode-node~=1
                                isConsecutive = 0;
                                deductBy = i-1;
                            else
                                node = nextNode;
                            end
                        else
                            isConsecutive = 0;
                            deductBy = i;
                        end
                    end
                    
                    obj.edgeList(obj.edgeList(:,1)>node,1) = obj.edgeList(obj.edgeList(:,1)>node,1)-deductBy;
                    obj.edgeList(obj.edgeList(:,2)>node,2) = obj.edgeList(obj.edgeList(:,2)>node,2)-deductBy;
                    
                    % Remove entry
                    idxNodesOutside = idxNodesOutside(idxNodesOutside>node);
                    if ~isempty(idxNodesOutside)
                        idxNodesOutside = idxNodesOutside-deductBy;
                    end
                    
                    % Print progress
                    if ~isFirst
                        fprintf(repmat('\b', 1, 7));
                    else
                        isFirst = 0;
                    end
                    fprintf('%6.2f%%', (1-length(idxNodesOutside)/nodeOutsideNum)*100);
                end
                fprintf('\n');
                
                % Save data
                obj.saveState(FileLocation.filterCut(obj), ...
                    'variables', {'idxNodesInside','idxEdgesOutside','nodeBoundary','edgeList'});
            end
            
            % Remove the nodes outside the new boundary
            obj.nodeBoundary = obj.nodeBoundary(obj.idxNodesInside==1,:);
            obj.nodeCoord = obj.nodeCoord(obj.idxNodesInside==1,:);
            obj.nodeSurfaceArea = obj.nodeSurfaceArea(obj.idxNodesInside==1);
            obj.nodeVolume = obj.nodeVolume(obj.idxNodesInside==1);
            
            % Remove edges outside the boundary
            obj.edgeArea = obj.edgeArea(obj.idxEdgesOutside==0);
            obj.edgeCoord = obj.edgeCoord(obj.idxEdgesOutside==0,:);
            obj.edgeDiameter = obj.edgeDiameter(obj.idxEdgesOutside==0);
            obj.edgeLength = obj.edgeLength(obj.idxEdgesOutside==0);
            obj.edgeNormal = obj.edgeNormal(obj.idxEdgesOutside==0,:);
            obj.edgeParticles = obj.edgeParticles(obj.idxEdgesOutside==0,:);
            obj.edgeReduced = obj.edgeReduced(obj.idxEdgesOutside==0);
            
            % Update nodeNum and edgeNum
            obj.nodeNum = size(obj.nodeBoundary,1);
            obj.edgeNum = size(obj.edgeArea,1);
        end
        
        % Shortest path
        function [P,d] = determineShortestPath(obj, weight)
            % Only for zpos
            assert(strcmp(obj.networkParameters.flowDirection, 'zpos'));
            
            % Number of inlets and outlets
            inletNum = sum(obj.nodeInlet);
            
            % Indices of inlets and outlets
            idxInlet = find(obj.nodeInlet==1);
            idxOutlet = find(obj.nodeOutlet==1);
            
            % Create a graph object
            switch weight
                case 'length'
                    G = graph(obj.edgeList(:,1),obj.edgeList(:,2),obj.edgeLength);
                case 'number'
                    G = graph(obj.edgeList(:,1),obj.edgeList(:,2));
            end
            
            % Initialise d
            [P,d] = shortestpath(G, idxInlet(1), idxOutlet(1));
            
            if strcmp(obj.networkParameters.networkType,'actual')
                % Min and max of x- and y-coordinates
                minX = min(obj.nodeCoord(:,1));
                maxX = max(obj.nodeCoord(:,1));
                minY = min(obj.nodeCoord(:,2));
                maxY = max(obj.nodeCoord(:,2));
                if or(obj.filterParameters.size(1)==1, obj.filterParameters.size(1)==0)
                    fraction = 0.05;
                elseif obj.filterParameters.size(1)==2
                    fraction = 0.025;
                end
                intervalX = fraction*(maxX-minX);
                intervalY = fraction*(maxY-minY);
                
                % Determine the shortest path from inlet to outlet
                fprintf('Finding shortest path ');
                for i = 1:inletNum
                    % Current inlet
                    inlet = idxInlet(i);
                    
                    % Determine the outlets close to x and y coordinates
                    inletX = obj.nodeCoord(inlet,1);
                    inletY = obj.nodeCoord(inlet,2);
                    idxOutlet = obj.nodeOutlet;
                    idxOutlet = and(idxOutlet, obj.nodeCoord(:,1)>inletX-intervalX);
                    idxOutlet = and(idxOutlet, obj.nodeCoord(:,1)<inletX+intervalX);
                    idxOutlet = and(idxOutlet, obj.nodeCoord(:,2)>inletY-intervalY);
                    idxOutlet = and(idxOutlet, obj.nodeCoord(:,2)<inletY+intervalY);
                    outletNum = sum(idxOutlet);
                    idxOutlet = find(idxOutlet==1);
                    
                    for j = 1:outletNum
                        % Current outlet
                        outlet = idxOutlet(j);
                        
                        % Determine the shortest path for this combination
                        [PTemp,dTemp] = shortestpath(G, inlet, outlet);
                        
                        % Store if dTemp is smaller than d
                        if dTemp<d
                            d = dTemp;
                            P = PTemp;
                        end
                    end
                    
                    if i~=1
                        fprintf(repmat('\b', 1, 7));
                    end
                    fprintf('%6.2f%%', i/inletNum*100);
                end
                fprintf('\n');
            else
                P = [];
                d = obj.edgeLength(1)*(obj.networkParameters.nz-1);
            end
        end
        function [obj] = shortestPathActual(obj)
            if isfile(FileLocation.filterShortestPath(obj))
                % Load data if available
                obj = obj.loadState(FileLocation.filterShortestPath(obj), ...
                    'variables', {'shortestPathZpos'});
                
                % Scale
                obj.shortestPathZpos.d = obj.shortestPathZpos.d*obj.filterParameters.minCoarse;
            else
                % Determine the shortest path
                [P,d] = obj.determineShortestPath('length');
                
                % Save
                obj.shortestPathZpos.P = P;
                obj.shortestPathZpos.d = d/obj.filterParameters.minCoarse;
                obj.saveState(FileLocation.filterShortestPath(obj), ...
                    'variables', {'shortestPathZpos'});
                obj.shortestPathZpos.d = d;
            end
        end
        
        % Generate network topology (actual)
        function [obj] = generateNetworkTopologyActual(obj)
            if isfile(FileLocation.filterNetwork(obj))
                obj = obj.loadState(FileLocation.filterNetwork(obj), ...
                    'variables', {'nodeInlet', 'nodeLattice', 'nodeOutlet', ...
                    'nodeConnectivity', 'nodeEdgeTo', 'nodeNodeTo', 'nodeEdgeToPos', 'nodeEdgeDirection','nodeList'});
                
                % Update nodeEdgeDirection
                if and(isempty(obj.nodeEdgeDirection), strcmp(obj.networkParameters.directedOrUndirected,'undirected'))
                    % Generate nodeEdgeDirection
                    obj = obj.generateNodeEdgeDirection();
                    
                    % Save
                    obj.saveState(FileLocation.filterNetwork(obj), 'variables', {'nodeInlet', 'nodeLattice', 'nodeOutlet', ...
                        'nodeConnectivity', 'nodeEdgeTo', 'nodeNodeTo', 'nodeEdgeToPos', 'nodeEdgeDirection', 'nodeList'});
                end
                
                % Update nodeList
                if 0 %and(isempty(obj.nodeList), all(obj.filterParameters.size==[1,1,1]))
                    % Generate nodeList
                    obj = obj.determineNodeList(); %#ok<UNRCH>
                    
                    % Save
                    obj.saveState(FileLocation.filterNetwork(obj), 'variables', {'nodeInlet', 'nodeLattice', 'nodeOutlet', ...
                        'nodeConnectivity', 'nodeEdgeTo', 'nodeNodeTo', 'nodeEdgeToPos', 'nodeEdgeDirection', 'nodeList'});
                end
            else
                if strcmp(obj.networkParameters.directedOrUndirected, 'undirected')
                    % Determine the connectivity of each node
                    obj.nodeConnectivity = zeros(obj.nodeNum,1);
                    if strcmp(obj.filterParameters.triangulationMethod, 'Mod')
                        fprintf('Determining node connectivity ');
                        for i = 1:obj.nodeNum
                            idxEdges = [find(obj.edgeList(:,1)==i);find(obj.edgeList(:,2)==i)];
                            obj.nodeConnectivity(i) = length(idxEdges);
                            if i~=1
                                fprintf(repmat('\b', 1, 7));
                            end
                            fprintf('%6.2f%%', i/obj.nodeNum*100);
                        end
                        fprintf('\n');
                        obj.nodeEdgeTo = zeros(obj.nodeNum,max(obj.nodeConnectivity));
                        obj.nodeNodeTo = obj.nodeEdgeTo;
                    else
                        obj.nodeEdgeTo = zeros(obj.nodeNum,4);
                        obj.nodeNodeTo = obj.nodeEdgeTo;
                    end
                    
                    % Determine the neighbours of each node
                    fprintf('Generating neighbour lists ');
                    for i = 1:obj.nodeNum
                        idxEdges = [find(obj.edgeList(:,1)==i);find(obj.edgeList(:,2)==i)];
                        idxNodes = [obj.edgeList(obj.edgeList(:,1)==i,2);obj.edgeList(obj.edgeList(:,2)==i,1)];
                        if strcmp(obj.filterParameters.triangulationMethod, 'Del')
                            obj.nodeConnectivity(i) = length(idxEdges);
                        end
                        obj.nodeEdgeTo(i,1:obj.nodeConnectivity(i)) = idxEdges;
                        obj.nodeNodeTo(i,1:obj.nodeConnectivity(i)) = idxNodes;
                        if i~=1
                            fprintf(repmat('\b', 1, 7));
                        end
                        fprintf('%6.2f%%', i/obj.nodeNum*100);
                    end
                    fprintf('\n');
                    
                    % Generate nodeEdgeDirection
                    obj = obj.generateNodeEdgeDirection();
                    
                    % Determine node list
                    obj = obj.determineNodeList;
                    
                    % Current flow direction
                    currentFlowDirection = obj.networkParameters.flowDirection;
                    
                    % For each flowDirection
                    flowDirectionList = {'xpos','xneg','ypos','yneg','zpos','zneg'};
                    for i = 1:length(flowDirectionList)
                        % Determine the indices
                        switch flowDirectionList{i}(1)
                            case 'x'
                                idxCoord = 1;
                            case 'y'
                                idxCoord = 2;
                            case 'z'
                                idxCoord = 3;
                        end
                        idxOdd = idxCoord*2-1;
                        switch flowDirectionList{i}(2:4)
                            case 'pos'
                                idxBoundary = [idxOdd, idxOdd+1];
                                coef = 1;
                            case 'neg'
                                idxBoundary = [idxOdd+1, idxOdd];
                                coef = -1;
                        end
                        
                        % Determine the inlets and the outlets
                        obj.nodeInlet = obj.nodeBoundary(:,idxBoundary(1))>0;
                        obj.nodeOutlet = obj.nodeBoundary(:,idxBoundary(2))>0;
                        
                        % Determine nodeLattice
                        obj.nodeLattice = zeros(obj.nodeNum,1);
                        n = 20;
                        minCoord = min(obj.nodeCoord(:,idxCoord));
                        maxCoord = max(obj.nodeCoord(:,idxCoord));
                        intervalCoord = linspace(minCoord, maxCoord, n+1);
                        for idxBoundary = 1:obj.nodeNum
                            obj.nodeLattice(idxBoundary) = find(obj.nodeCoord(idxBoundary,idxCoord)<=intervalCoord(2:end),1);
                        end
                        if strcmp(flowDirectionList{i}(2:4),'neg')
                            obj.nodeLattice = 21-obj.nodeLattice;
                        end
                        
                        % Determine nodeEdgeToPos
                        obj.nodeEdgeToPos = zeros(size(obj.nodeEdgeTo));
                        for j = 1:obj.nodeNum
                            idxNodes = obj.nodeNodeTo(j,1:obj.nodeConnectivity(j));
                            for k = 1:length(idxNodes)
                                if (coef*(obj.nodeCoord(idxNodes(k),idxCoord)-obj.nodeCoord(j,idxCoord)))>0
                                    obj.nodeEdgeToPos(j,k) = 1;
                                end
                            end
                        end
                        
                        % Save data
                        obj.networkParameters.flowDirection = flowDirectionList{i};
                        obj.saveState(FileLocation.filterNetwork(obj), 'variables', {'nodeInlet', 'nodeLattice', 'nodeOutlet', ...
                            'nodeConnectivity', 'nodeEdgeTo', 'nodeNodeTo', 'nodeEdgeToPos', 'nodeEdgeDirection', 'nodeList'});
                    end
                    
                    % Load the required flowDirection
                    obj.networkParameters.flowDirection = currentFlowDirection;
                    obj = obj.loadState(FileLocation.filterNetwork(obj), ...
                        'variables', {'nodeInlet', 'nodeLattice', 'nodeOutlet', ...
                        'nodeConnectivity', 'nodeEdgeTo', 'nodeNodeTo', 'nodeEdgeToPos', 'nodeEdgeDirection', 'nodeList'});
                else
                    % Determine the indices
                    switch obj.networkParameters.flowDirection(1)
                        case 'x'
                            idxCoord = 1;
                        case 'y'
                            idxCoord = 2;
                        case 'z'
                            idxCoord = 3;
                    end
                    idxBoundary = idxCoord*2-1;
                    switch obj.networkParameters.flowDirection(2:4)
                        case 'pos'
                            idxBoundary = [idxBoundary, idxBoundary+1];
                            coef = 1;
                        case 'neg'
                            idxBoundary = [idxBoundary+1, idxBoundary];
                            coef = -1;
                    end
                    
                    % Determine the connectivity of each node
                    obj.nodeConnectivity = zeros(obj.nodeNum,1);
                    fprintf('Determining node connectivity ');
                    for i = 1:obj.nodeNum
                        idxNodes = [obj.edgeList(obj.edgeList(:,1)==i,2);obj.edgeList(obj.edgeList(:,2)==i,1)];
                        isParallel = zeros(length(idxNodes),1);
                        for j = 1:length(idxNodes)
                            if (coef*(obj.nodeCoord(idxNodes(j),idxCoord)-obj.nodeCoord(i,idxCoord)))>0
                                isParallel(j) = 1;
                            end
                        end
                        obj.nodeConnectivity(i) = sum(isParallel);
                        if i~=1
                            fprintf(repmat('\b', 1, 7));
                        end
                        fprintf('%6.2f%%', i/obj.nodeNum*100);
                    end
                    fprintf('\n');
                    
                    % Determine the neighbours of each node
                    obj.nodeEdgeTo = zeros(obj.nodeNum,max(obj.nodeConnectivity));
                    obj.nodeNodeTo = obj.nodeEdgeTo;
                    fprintf('Generating neighbour lists ');
                    for i = 1:obj.nodeNum
                        idxEdges = [find(obj.edgeList(:,1)==i);find(obj.edgeList(:,2)==i)];
                        idxNodes = [obj.edgeList(obj.edgeList(:,1)==i,2);obj.edgeList(obj.edgeList(:,2)==i,1)];
                        isParallel = zeros(length(idxNodes),1);
                        for j = 1:length(idxNodes)
                            if (coef*(obj.nodeCoord(idxNodes(j),idxCoord)-obj.nodeCoord(i,idxCoord)))>0
                                isParallel(j) = 1;
                            end
                        end
                        obj.nodeEdgeTo(i,1:obj.nodeConnectivity(i)) = idxEdges(isParallel==1);
                        obj.nodeNodeTo(i,1:obj.nodeConnectivity(i)) = idxNodes(isParallel==1);
                        if i~=1
                            fprintf(repmat('\b', 1, 7));
                        end
                        fprintf('%6.2f%%', i/obj.nodeNum*100);
                    end
                    fprintf('\n');
                    
                    % Determine node list
                    obj = obj.determineNodeList;
                    
                    % Determine the inlets and outlets
                    obj.nodeInlet = obj.nodeBoundary(:,idxBoundary(1))>0;
                    obj.nodeOutlet = obj.nodeBoundary(:,idxBoundary(2))>0;
                    
                    % Determine nodeLattice
                    obj.nodeLattice = zeros(obj.nodeNum,1);
                    n = 20;
                    minCoord = min(obj.nodeCoord(:,idxCoord));
                    maxCoord = max(obj.nodeCoord(:,idxCoord));
                    intervalCoord = linspace(minCoord, maxCoord, n+1);
                    for idxBoundary = 1:obj.nodeNum
                        obj.nodeLattice(idxBoundary) = find(obj.nodeCoord(idxBoundary,idxCoord)<=intervalCoord(2:end),1);
                    end
                    if strcmp(obj.networkParameters.flowDirection(2:4),'neg')
                        obj.nodeLattice = 21-obj.nodeLattice;
                    end
                    
                    % Save to improve efficiency
                    obj.saveState(FileLocation.filterNetwork(obj), 'variables', {'nodeInlet', 'nodeLattice', 'nodeOutlet', ...
                        'nodeConnectivity', 'nodeEdgeTo', 'nodeNodeTo', 'nodeList'});
                end
            end
        end
        
        % Generate nodeEdgeDirection to be used for flow calculations
        function [obj] = generateNodeEdgeDirection(obj)
            % Initialise nodeEdgeDirection
            obj.nodeEdgeDirection = zeros(size(obj.nodeEdgeTo));
            
            fprintf('Generating nodeEdgeDirection ');
            for i = 1:obj.nodeNum
                for j = 1:obj.nodeConnectivity(i)
                    if obj.edgeList(obj.nodeEdgeTo(i,j),1)==i
                        obj.nodeEdgeDirection(i,j) = 1;
                    else
                        obj.nodeEdgeDirection(i,j) = -1;
                    end
                end
                if i~=1
                    fprintf(repmat('\b', 1, 7));
                end
                fprintf('%6.2f%%', i/obj.nodeNum*100);
            end
            fprintf('\n');
        end
        
        % Randomly assign network geometry
        function [obj] = assignNetworkGeometry(obj)
            % Parameters
            cdfInterval = 1e4;
            
            % Create CDFs
            cdfEdgeDiameter = Statistics.createCDF(obj.edgeDiameter,cdfInterval);
            cdfEdgeArea = Statistics.createCDF(obj.edgeArea,cdfInterval);
            cdfNodeVolume = Statistics.createCDF(obj.nodeVolume,cdfInterval);
            cdfNodeSurfaceArea = Statistics.createCDF(obj.nodeSurfaceArea,cdfInterval);
            cdfEdgeDiameter = cdfEdgeDiameter(2:end,:);
            cdfEdgeArea = cdfEdgeArea(2:end,:);
            cdfNodeVolume = cdfNodeVolume(2:end,:);
            cdfNodeSurfaceArea = cdfNodeSurfaceArea(2:end,:);
            
            % Random generation of network properties
            r1 = rand(obj.nodeNum,1);
            r2 = rand(obj.edgeNum,1);
            obj.nodeSurfaceArea = cdfNodeSurfaceArea(ceil(length(cdfNodeSurfaceArea)*r1),2);
            obj.nodeVolume = cdfNodeVolume(ceil(length(cdfNodeVolume)*r1),2);
            obj.edgeArea = cdfEdgeArea(ceil(length(cdfEdgeArea)*r2),2);
            obj.edgeDiameter = cdfEdgeDiameter(ceil(length(cdfEdgeDiameter)*r2),2);
            
            % Determine the median of edgeLength
            if ~strcmp(obj.networkParameters.networkType, 'actual')
                edgeLength50 = median(obj.edgeLength);
                obj.edgeLength = ones(obj.edgeNum,1)*edgeLength50;
            end
            
            % Determine the median of edgeReduced
            edgeReduced50 = median(obj.edgeReduced);
            obj.edgeReduced = ones(obj.edgeNum,1)*edgeReduced50;
            
            % Determine nodeCoord
            if ~strcmp(obj.networkParameters.networkType, 'actual')
                % Determine the median of edgeLength
                l = edgeLength50;
                switch obj.networkParameters.networkType
                    case 'cubic'
                        distance = [l, l, l];
                    case 'square'
                        distance = [l, l, l];
                    case 'bcc'
                        distance = [l/2, l/2, l/sqrt(2)];
                    case 'fcc'
                        distance = [l/sqrt(2), l/sqrt(2), l/sqrt(2)];
                    case 'slantedcubic'
                        distance = [l*cosd(atand(1/sqrt(2)))*cosd(30), l*cosd(atand(1/sqrt(2)))*sind(30), l*sind(atand(1/sqrt(2)))];
                    case 'slantedsquare'
                        distance = [l/sqrt(2), l, l/sqrt(2)];
                    case 'cubicnp'
                        distance = [l, l, l];
                    case 'squarenp'
                        distance = [l, l, l];
                    case 'bccnp'
                        distance = [l/2, l/2, l/sqrt(2)];
                    case 'fccnp'
                        distance = [l/sqrt(2), l/sqrt(2), l/sqrt(2)];
                    case 'slantedcubicnp'
                        distance = [l*cosd(atand(1/sqrt(2)))*cosd(30), l*cosd(atand(1/sqrt(2)))*sind(30), l*sind(atand(1/sqrt(2)))];
                    case 'slantedsquarenp'
                        distance = [l/sqrt(2), l, l/sqrt(2)];
                end
                obj.nodeCoord = (obj.nodeLattice).*repmat(distance, obj.nodeNum, 1);
            end
        end
        
        % Determine CSD by num
        function [obj] = determineCSDNum(obj)
            obj.csdNum = Statistics.createCDF(obj.edgeDiameter, 1e4);
        end
        
        % Obtain pairs of nodes that are a certain distance away, used for
        % two-point spatial correlation analysis
        function [obj] = determineNodeList(obj, varargin)
            % Input parser
            parser = inputParser;
            addParameter(parser, 'distanceUpTo', 5);
            parse(parser, varargin{:});
            distanceUpTo = parser.Results.distanceUpTo;
            
            % Create a graph object
            G = graph(obj.edgeList(:,1),obj.edgeList(:,2));
            
            for i = 1:distanceUpTo
                % Initialise array
                nodeComb = [];
                
                % For each node
                fprintf('Determining node list for distance %i ', i);
                for j = 1:obj.nodeNum
                    % Nodes that are a specific distance away
                    idxNodesCurrent = nearest(G,j,i);
                    if i~=1
                        idxNodesLess = nearest(G,j,i-1);
                        idxNodesCurrent = idxNodesCurrent(~ismember(idxNodesCurrent,idxNodesLess));
                    end
                    
                    % Add the combinations of nodes
                    nodeComb = [nodeComb; [j*ones(length(idxNodesCurrent),1),idxNodesCurrent]]; %#ok<AGROW>
                    
                    if j~=1
                        fprintf(repmat('\b', 1, 7));
                    end
                    fprintf('%6.2f%%', j/obj.nodeNum*100);
                end
                fprintf('\n');
                
                % Sort and search for unique combination
                nodeComb = sort(nodeComb,2);
                nodeComb = unique(nodeComb,'rows');
                
                obj.nodeList.(sprintf('distance%i',i)) = nodeComb;
            end
        end
        
        % Determine the CSD from PSD based on PSD by surface area and
        % relative density (Silveira et. al. 1975, Locke et. al. 2001,
        % Indraratna et. al. 2007)
        function [obj] = analyticalCSDDensest(obj)
            % Determine psd by surface area
            psdVol = Statistics.createCDF(obj.particleRadius*2,1e4);
            psdSATemp = Statistics.convertPSDVolToArea(psdVol);
            
            % Lower resolution
            psdSA = zeros(101,2);
            idxTemp = (1:100:length(psdSATemp))';
            psdSA(:,1) = psdSATemp(idxTemp,1);
            psdSA(:,2) = psdSATemp(idxTemp,2);
            
            % Determine the probability and mean diameter
            pdfSA = zeros(100,2);
            pdfSA(:,1) = diff(psdSA(:,1));
            pdfSA(:,2) = diff(psdSA(:,2))/2;
            pdfSA(:,2) = pdfSA(:,2)+psdSA(1:100,2);
            
            % Determine a list of index (combination, can repeat)
            l = size(pdfSA,1);
            idxList = (combvec(1:l,1:l,1:l))';
            idxList = unique(sort(idxList,2),'rows'); %#ok<UDIM>
            idxNum = size(idxList,1);
            
            % Initialise array(s)
            csList = zeros(idxNum,2);
            
            % For each combination
            fprintf('Generating densest analytical CSD ');
            for i = 1:idxNum
                % Current indices
                idx = idxList(i,:);
                
                % Current combination of diameter
                d1 = pdfSA(idx(1),2);
                d2 = pdfSA(idx(2),2);
                d3 = pdfSA(idx(3),2);
                
                % Solve for the constriction size
                csList(i,1) = inscribedCircleDensest(d1, d2, d3);
                
                % Determine the number of times each element is repeated
                uniqueIdx = unique(idx);
                uniqueNum = length(uniqueIdx);
                coef = factorial(3);
                for j = 1:uniqueNum
                    coef = coef/factorial(sum(idx==uniqueIdx(j)));
                end
                
                % Determine the probability
                csList(i,2) = coef*((pdfSA(1,1))^3);
                
                if i~=1
                    fprintf(repmat('\b', 1, 7));
                end
                fprintf('%6.2f%%', i/idxNum*100);
            end
            fprintf('\n');
            
            % Determine csd
            obj.analyticalCSD.densest = Statistics.createCDF(csList(:,1),1e4,'frequency',csList(:,2));
        end
        function [obj] = analyticalCSDLoosest(obj)
            % Determine psd by surface area
            psdVol = Statistics.createCDF(obj.particleRadius*2,1e4);
            psdSATemp = Statistics.convertPSDVolToArea(psdVol);
            
            % Lower the resolution/increment
            psdSA = zeros(51,2);
            idxTemp = (1:200:length(psdSATemp))';
            psdSA(:,1) = psdSATemp(idxTemp,1);
            psdSA(:,2) = psdSATemp(idxTemp,2);
            
            % Determine the probability and mean diameter
            pdfSA = zeros(50,2);
            pdfSA(:,1) = diff(psdSA(:,1));
            pdfSA(:,2) = diff(psdSA(:,2))/2;
            pdfSA(:,2) = pdfSA(:,2)+psdSA(1:50,2);
            
            % Determine a list of index (combination, can repeat)
            l = size(pdfSA,1);
            idxList = (combvec(1:l,1:l,1:l,1:l))';
            idxList = unique(sort(idxList,2),'rows'); %#ok<UDIM>
            idxNum = size(idxList,1);
            
            % Initialise array(s)
            csList = zeros(idxNum,2);
            
            % For each combination
            fprintf('Generating loosest analytical CSD ');
            for i = 1:idxNum
                % Current indices
                idx = idxList(i,:);
                
                % Current combination of diameter
                d1 = pdfSA(idx(1),2);
                d2 = pdfSA(idx(2),2);
                d3 = pdfSA(idx(3),2);
                d4 = pdfSA(idx(4),2);
                
                % Solve for the constriction size
                csList(i,1) = inscribedCircleLoosest(d1, d2, d3, d4);
                
                % Determine the number of times each element is repeated
                uniqueIdx = unique(idx);
                uniqueNum = length(uniqueIdx);
                coef = factorial(4);
                for j = 1:uniqueNum
                    coef = coef/factorial(sum(idx==uniqueIdx(j)));
                end
                
                % Determine the probability
                csList(i,2) = coef*((pdfSA(1,1))^4);
                
                if i~=1
                    fprintf(repmat('\b', 1, 7));
                end
                fprintf('%6.2f%%', i/idxNum*100);
            end
            fprintf('\n');
            
            % Determine csd
            obj.analyticalCSD.loosest = Statistics.createCDF(csList(:,1),1e4,'frequency',csList(:,2));
        end
    end
end

% PRIVATE FUNCTIONS
% ===== Constriction Size =================================================
% Largest insribed sphere in 3D
function [x] = inscribedCircle(x1, x2, x3)
% The function solves for the centroid and radius of inscribed circle given
% the centroids of surrounding three particles

f = @(x)f1(x, x1, x2, x3);
x0 = (x1+x2+x3)/3;
x0(4) = x0(4)/2;
x = fsolve(f, x0, optimoptions('fsolve','Display','off'));
end
function [y] = f1(x, x1, x2, x3)
% Simultaneous equations to solve for centroid and radius of inscribed
% circle

% Refer to equations 1&2 in Reboul, N., Vincens, E. and Cambou, B., 2010.
% A computational procedure to assess the distribution of constriction
% sizes for an assembly of spheres. Computers and Geotechnics, 37(1-2),
% pp.195-206.

y = zeros(4,1);
y(1) = sum((x1(1:3)-x(1:3)).^2)-((x1(4)+x(4))^2);
y(2) = sum((x2(1:3)-x(1:3)).^2)-((x2(4)+x(4))^2);
y(3) = sum((x3(1:3)-x(1:3)).^2)-((x3(4)+x(4))^2);
A = [x(1)-x1(1), x2(1)-x1(1), x3(1)-x1(1); x(2)-x1(2), x2(2)-x1(2), x3(2)-x1(2); x(3)-x1(3), x2(3)-x1(3), x3(3)-x1(3)];
y(4) = det(A);
end

% Largest inscribed circle in 2D (Three particles arrangement) (Silveira
% et. al. 1975)
function [d] = inscribedCircleDensest(d1, d2, d3)

f = @(d)f2(d, d1, d2, d3);
x0 = 0.15;
d = fzero(f, x0);
end
function [y] = f2(d, d1, d2, d3)
% Refer to Equation 1 in Indraratna et. al. (2007) Constriction-Based
% Retention Criterion for Granular Filter Design

y = (((2/d1)^2)+((2/d2)^2)+((2/d3)^2)+((2/d)^2)) - (0.5*(((2/d1)+(2/d2)+(2/d3)+(2/d))^2));
end

% Largest inscribed circle in 2D (Four particles arrangement) (Silveira
% et. al. 1975)
function [d] = inscribedCircleLoosest(d1, d2, d3, d4)
% Determine the min alpha (using Law of Cosine)
a = d1/2+d2/2;
b = d1/2+d4/2;
c = d2/2+d4/2;
alphaMin = acos(((a^2)+(b^2)-(c^2))/(2*a*b));

% Determine the max alpha (using Law of Cosine)
a = d1/2+d3/2;
b = d1/2+d2/2;
c = d2/2+d3/2;
alphaMax1 = acos(((a^2)+(b^2)-(c^2))/(2*a*b));
a = d1/2+d3/2;
b = d1/2+d4/2;
c = d4/2+d3/2;
alphaMax2 = acos(((a^2)+(b^2)-(c^2))/(2*a*b));
alphaMax = alphaMax1+alphaMax2;

% Number of alpha values
trialNum = 1e3;
alphaList = linspace(alphaMin, alphaMax, trialNum);

% Initialise array(s)
Sc = zeros(trialNum,1);
a = d1/2+d2/2;
b = d1/2+d4/2;
d = d2/2+d3/2;
e = d4/2+d3/2;

% For each value of alpha
for i = 1:trialNum
    % Determine beta, gamma and delta
    alpha = alphaList(i);
    c = sqrt((a^2)+(b^2)-2*a*b*cos(alpha));
    beta = acos(((a^2)+(c^2)-(b^2))/(2*a*c))+acos(((c^2)+(d^2)-(e^2))/(2*c*d));
    delta = acos(((b^2)+(c^2)-(a^2))/(2*b*c))+acos(((c^2)+(e^2)-(d^2))/(2*c*e));
    gamma = acos(((d^2)+(e^2)-(c^2))/(2*d*e));
    
    % Determine Sc 
    Sc(i) = f3(d1, d2, d3, d4, alpha, beta, gamma, delta);
end

% Maximum Sc
Sc = max(Sc);
assert(Sc>0)

% Determine the equivalent d (Equation 3 in Indraratna et. al. 2007)
d = sqrt(4*Sc/pi);
assert(isreal(d))
end
function [Sc] = f3(d1, d2, d3, d4, alpha, beta, gamma, delta)
% Equation 2 in Indraratna et. al. 2007

Sc = (1/8)*((d1+d2)*(d1+d4)*sin(alpha)+(d2+d3)*(d3+d4)*sin(gamma)-(alpha*(d1^2)+beta*(d2^2)+gamma*(d3^2)+delta*(d4^2)));
end

% ===== Generate Regular Networks =========================================
function [nodeNum, edgeNum, edgeList, nodeEdgeTo, nodeNodeTo, nodeConnectivity, nodeLattice, nodeInlet, nodeOutlet] = generateCubicNetwork(nx, ny, nz, directedOrUndirected, periodic)
% This function generates a cubic lattice

% Define the number of nodes, and the list of nodes (with one indexing,
% i.e. node index starts from 1). Loop through and assign a lattice
% coordinate to each node. Determine the lattice coordiantes of all edges
% connected to a node
nodeNum = nx*ny*nz;
nodeLattice = zeros(nodeNum, 3);
if periodic
    edgeLattice = zeros(nodeNum*3+ny*nz+nx*nz, 6);
else
    edgeLattice = zeros(nodeNum*3, 6);
end
nodeId = 1;
edgeId = 1;
for z = 1:nz
    for x = 1:nx
        for y = 1:ny
            nodeLattice(nodeId, :) = [x, y, z]; nodeId = nodeId + 1;
            edgeLattice(edgeId, :) = [x, y, z, x+1, y, z]; edgeId = edgeId + 1;
            edgeLattice(edgeId, :) = [x, y, z, x, y+1, z]; edgeId = edgeId + 1;
            edgeLattice(edgeId, :) = [x, y, z, x, y, z+1]; edgeId = edgeId + 1;
        end
    end
end
assert(nodeId==(nodeNum+1))
assert(edgeId==(nodeNum*3+1))

% Introduce periodic boundary
if periodic
    for y = 1:ny
        for z = 1:nz
            edgeLattice(edgeId, :) = [nx, y, z, 1, y, z]; edgeId = edgeId + 1;
        end
    end
    for x = 1:nx
        for z = 1:nz
            edgeLattice(edgeId, :) = [x, ny, z, x, 1, z]; edgeId = edgeId + 1;
        end
    end
    assert(edgeId==(nodeNum*3+ny*nz+nx*nz+1))
end

% Remove excess boundary edges. These are those entries with edgeLattice
% coordinate of 0, nx+1, ny+1, nz+1.
remX = edgeLattice(:,4)>nx;
remY = edgeLattice(:,5)>ny;
remZ = edgeLattice(:,6)>nz;
edgeLattice = edgeLattice((remX + remY + remZ)==0, :);

% Identify the nodes corresponding to each edge based on the lattice
% coordinates.
edgeLatticeFrom = edgeLattice(:, 1:3);
edgeLatticeTo = edgeLattice(:, 4:6);
[bool, edgeIndexFrom] = ismember(edgeLatticeFrom, nodeLattice, 'rows'); assert(all(bool));
[bool, edgeIndexTo] = ismember(edgeLatticeTo, nodeLattice, 'rows'); assert(all(bool));
edgeList = [edgeIndexFrom, edgeIndexTo];
edgeNum = size(edgeList,1);

% Determine the neighbours of each nodes
nodeConnectivity = zeros(nodeNum,1);
fprintf('Generating neighbour lists ');
switch directedOrUndirected
    case 'directed'
        nodeEdgeTo = zeros(nodeNum,3);
        nodeNodeTo = nodeEdgeTo;
        for i = 1:nodeNum
            idxEdges = find(edgeList(:,1)==i);
            idxNodes = edgeList(edgeList(:,1)==i,2)';
            nodeEdgeTo(i,1:length(idxEdges)) = idxEdges;
            nodeNodeTo(i,1:length(idxEdges)) = idxNodes;
            nodeConnectivity(i) = length(idxEdges);
            if i~=1
                fprintf(repmat('\b', 1, 7));
            end
            fprintf('%6.2f%%', i/nodeNum*100);
        end
    case 'undirected'
        nodeEdgeTo = zeros(nodeNum,6);
        nodeNodeTo = nodeEdgeTo;
        for i = 1:nodeNum
            idxEdges = [find(edgeList(:,1)==i);find(edgeList(:,2)==i)];
            idxNodes = [edgeList(edgeList(:,1)==i,2);edgeList(edgeList(:,2)==i,1)];
            nodeEdgeTo(i,1:length(idxEdges)) = idxEdges;
            nodeNodeTo(i,1:length(idxEdges)) = idxNodes;
            nodeConnectivity(i) = length(idxEdges);
            if i~=1
                fprintf(repmat('\b', 1, 7));
            end
            fprintf('%6.2f%%', i/nodeNum*100);
        end
end
fprintf('\n');

% Determine the inlets and the outlets
nodeInlet = zeros(nodeNum,1);
nodeOutlet = nodeInlet;
nodeInlet(1:(nx*ny)) = 1;
nodeOutlet((nz-1)*(nx*ny)+1:end) = 1;
end
function [nodeNum, edgeNum, edgeList, nodeEdgeTo, nodeNodeTo, nodeConnectivity, nodeLattice, nodeInlet, nodeOutlet] = generateSquareNetwork(nx, nz, directedOrUndirected, periodic)
% This function generates a square lattice

% Define the number of nodes, and the list of nodes (with one indexing,
% i.e. node index starts from 1). Loop through and assign a lattice
% coordinate to each node. Determine the lattice coordiantes of all edges
% connected to a node
nodeNum = nx*nz;
nodeLattice = zeros(nodeNum, 3);
if periodic
    edgeLattice = zeros(nodeNum*2+nz, 6);
else
    edgeLattice = zeros(nodeNum*2, 6);
end
nodeId = 1;
edgeId = 1;
for z = 1:nz
    for x = 1:nx
        nodeLattice(nodeId, :) = [x, 1, z]; nodeId = nodeId + 1;
        edgeLattice(edgeId, :) = [x, 1, z, x+1, 1, z]; edgeId = edgeId + 1;
        edgeLattice(edgeId, :) = [x, 1, z, x, 1, z+1]; edgeId = edgeId + 1;
    end
end
assert(nodeId==(nodeNum+1))
assert(edgeId==(nodeNum*2+1))

% Introduce periodic boundary
if periodic
    for z = 1:nz
        edgeLattice(edgeId, :) = [nx, 1, z, 1, 1, z]; edgeId = edgeId + 1;
    end
    assert(edgeId==(nodeNum*2+nz+1))
end

% Remove excess boundary edges. These are those entries with edgeLattice
% coordinate of 0, nx+1, ny+1, nz+1.
remX = edgeLattice(:,4)>nx;
remZ = edgeLattice(:,6)>nz;
edgeLattice = edgeLattice((remX + remZ)==0, :);

% Identify the nodes corresponding to each edge based on the lattice
% coordinates.
edgeLatticeFrom = edgeLattice(:, 1:3);
edgeLatticeTo = edgeLattice(:, 4:6);
[bool, edgeIndexFrom] = ismember(edgeLatticeFrom, nodeLattice, 'rows'); assert(all(bool));
[bool, edgeIndexTo] = ismember(edgeLatticeTo, nodeLattice, 'rows'); assert(all(bool));
edgeList = [edgeIndexFrom, edgeIndexTo];
edgeNum = size(edgeList,1);

% Determine the neighbours of each nodes
nodeConnectivity = zeros(nodeNum,1);
fprintf('Generating neighbour lists ');
switch directedOrUndirected
    case 'directed'
        nodeEdgeTo = zeros(nodeNum,2);
        nodeNodeTo = nodeEdgeTo;
        for i = 1:nodeNum
            idxEdges = find(edgeList(:,1)==i);
            idxNodes = edgeList(edgeList(:,1)==i,2)';
            nodeEdgeTo(i,1:length(idxEdges)) = idxEdges;
            nodeNodeTo(i,1:length(idxEdges)) = idxNodes;
            nodeConnectivity(i) = length(idxEdges);
            if i~=1
                fprintf(repmat('\b', 1, 7));
            end
            fprintf('%6.2f%%', i/nodeNum*100);
        end
    case 'undirected'
        nodeEdgeTo = zeros(nodeNum,4);
        nodeNodeTo = nodeEdgeTo;
        for i = 1:nodeNum
            idxEdges = [find(edgeList(:,1)==i);find(edgeList(:,2)==i)];
            idxNodes = [edgeList(edgeList(:,1)==i,2);edgeList(edgeList(:,2)==i,1)];
            nodeEdgeTo(i,1:length(idxEdges)) = idxEdges;
            nodeNodeTo(i,1:length(idxEdges)) = idxNodes;
            nodeConnectivity(i) = length(idxEdges);
            if i~=1
                fprintf(repmat('\b', 1, 7));
            end
            fprintf('%6.2f%%', i/nodeNum*100);
        end
end
fprintf('\n');

% Determine the inlets and the outlets
nodeInlet = zeros(nodeNum,1);
nodeOutlet = nodeInlet;
nodeInlet(1:nx) = 1;
nodeOutlet((nz-1)*nx+1:end) = 1;
end
function [nodeNum, edgeNum, edgeList, nodeEdgeTo, nodeNodeTo, nodeConnectivity, nodeLattice, nodeInlet, nodeOutlet] = generateBCCNetwork(nx, ny, nz, directedOrUndirected, periodic)
% This function generates a BCC network

% Define the number of nodes, and the list of nodes (with one indexing,
% i.e. node index starts from 1). Loop through and assign a lattice
% coordinate to each node. Determine the lattice coordiantes of all edges
% connected to a node.
nodeNum = nx*ny*nz;
nodeLattice = zeros(nodeNum, 3);
edgeLattice = zeros(nodeNum*4-nx*ny*4, 6);
nodeId = 1;
edgeId = 1;
for z = 1:nz
    if mod(z,2)==1
        for x = 1:nx
            for y = 1:ny
                nodeLattice(nodeId, :) = [2*(x-1)+1, 2*(y-1)+1, z]; nodeId = nodeId + 1;
                if z~=nz
                    edgeLattice(edgeId, :) = [2*(x-1)+1, 2*(y-1)+1, z, 2*(x-1), 2*(y-1), z+1]; edgeId = edgeId + 1;
                    edgeLattice(edgeId, :) = [2*(x-1)+1, 2*(y-1)+1, z, 2*x, 2*(y-1), z+1]; edgeId = edgeId + 1;
                    edgeLattice(edgeId, :) = [2*(x-1)+1, 2*(y-1)+1, z, 2*(x-1), 2*y, z+1]; edgeId = edgeId + 1;
                    edgeLattice(edgeId, :) = [2*(x-1)+1, 2*(y-1)+1, z, 2*x, 2*y, z+1]; edgeId = edgeId + 1;
                end
            end
        end
    else
        for x = 1:nx
            for y = 1:ny
                nodeLattice(nodeId, :) = [2*x, 2*y, z]; nodeId = nodeId + 1;
                if z~=nz
                    edgeLattice(edgeId, :) = [2*x, 2*y, z, 2*(x-1)+1, 2*(y-1)+1, z+1]; edgeId = edgeId + 1;
                    edgeLattice(edgeId, :) = [2*x, 2*y, z, 2*x+1, 2*(y-1)+1, z+1]; edgeId = edgeId + 1;
                    edgeLattice(edgeId, :) = [2*x, 2*y, z, 2*(x-1)+1, 2*y+1, z+1]; edgeId = edgeId + 1;
                    edgeLattice(edgeId, :) = [2*x, 2*y, z, 2*x+1, 2*y+1, z+1]; edgeId = edgeId + 1;
                end
            end
        end
    end
end
assert(nodeId==(nodeNum+1));
assert(edgeId==(nodeNum*4-nx*ny*4+1));

% Introduce periodic boundary
if periodic
    edgeLattice(edgeLattice(:,4)==0,4) = 2*nx;
    edgeLattice(edgeLattice(:,5)==0,5) = 2*ny;
    edgeLattice(edgeLattice(:,4)==2*nx+1,4) = 1;
    edgeLattice(edgeLattice(:,5)==2*ny+1,5) = 1;
else
    edgeLattice = edgeLattice(edgeLattice(:,4)~=0,:);
    edgeLattice = edgeLattice(edgeLattice(:,5)~=0,:);
    edgeLattice = edgeLattice(edgeLattice(:,4)~=2*nx+1,:);
    edgeLattice = edgeLattice(edgeLattice(:,5)~=2*ny+1,:);
end

% Identify the nodes corresponding to each edge based on the lattice
% coordinates.
edgeLatticeFrom = edgeLattice(:, 1:3);
edgeLatticeTo = edgeLattice(:, 4:6);
[bool, edgeIndexFrom] = ismember(edgeLatticeFrom, nodeLattice, 'rows'); assert(all(bool));
[bool, edgeIndexTo] = ismember(edgeLatticeTo, nodeLattice, 'rows'); assert(all(bool));
edgeList = [edgeIndexFrom, edgeIndexTo];
edgeNum = size(edgeList,1);

% Determine the neighbours of each nodes
nodeConnectivity = zeros(nodeNum,1);
fprintf('Generating neighbour lists ');
switch directedOrUndirected
    case 'directed'
        nodeEdgeTo = zeros(nodeNum,4);
        nodeNodeTo = nodeEdgeTo;
        for i = 1:nodeNum
            idxEdges = find(edgeList(:,1)==i);
            idxNodes = edgeList(edgeList(:,1)==i,2)';
            nodeEdgeTo(i,1:length(idxEdges)) = idxEdges;
            nodeNodeTo(i,1:length(idxEdges)) = idxNodes;
            nodeConnectivity(i) = length(idxEdges);
            if i~=1
                fprintf(repmat('\b', 1, 7));
            end
            fprintf('%6.2f%%', i/nodeNum*100);
        end
    case 'undirected'
        nodeEdgeTo = zeros(nodeNum,8);
        nodeNodeTo = nodeEdgeTo;
        for i = 1:nodeNum
            idxEdges = [find(edgeList(:,1)==i);find(edgeList(:,2)==i)];
            idxNodes = [edgeList(edgeList(:,1)==i,2);edgeList(edgeList(:,2)==i,1)];
            nodeEdgeTo(i,1:length(idxEdges)) = idxEdges;
            nodeNodeTo(i,1:length(idxEdges)) = idxNodes;
            nodeConnectivity(i) = length(idxEdges);
            if i~=1
                fprintf(repmat('\b', 1, 7));
            end
            fprintf('%6.2f%%', i/nodeNum*100);
        end
end
fprintf('\n');

% Determine the inlets and the outlets
nodeInlet = zeros(nodeNum,1);
nodeOutlet = nodeInlet;
nodeInlet(1:(nx*ny)) = 1;
nodeOutlet((nz-1)*(nx*ny)+1:end) = 1;
end
function [nodeNum, edgeNum, edgeList, nodeEdgeTo, nodeNodeTo, nodeConnectivity, nodeLattice, nodeInlet, nodeOutlet] = generateFCCNetwork(nx, ny, nz, directedOrUndirected, periodic)
% This function generates a FCC network

% Define the number of nodes, and the list of nodes (with one indexing,
% i.e. node index starts from 1). Loop through and assign a lattice
% coordinate to each node. Determine the lattice coordiantes of all edges
% connected to a node.
nodeNum = nx*ny*nz*2;
nodeLattice = zeros(nodeNum, 3);
edgeLattice = zeros(nodeNum*6-nx*ny*2*4, 6);
nodeId = 1;
edgeId = 1;
for z = 1:nz
    if mod(z,2)==1
        for x = 1:nx
            for y = 1:ny
                nodeLattice(nodeId, :) = [2*x-1, 2*y-1, z]; nodeId = nodeId + 1;
                nodeLattice(nodeId, :) = [2*x, 2*y, z]; nodeId = nodeId + 1;
                edgeLattice(edgeId, :) = [2*x-1, 2*y-1, z, 2*x-2, 2*y, z]; edgeId = edgeId + 1;
                edgeLattice(edgeId, :) = [2*x-1, 2*y-1, z, 2*x, 2*y, z]; edgeId = edgeId + 1;
                edgeLattice(edgeId, :) = [2*x, 2*y, z, 2*x-1, 2*y+1, z]; edgeId = edgeId + 1;
                edgeLattice(edgeId, :) = [2*x, 2*y, z, 2*x+1, 2*y+1, z]; edgeId = edgeId + 1;
                if z~=nz
                    edgeLattice(edgeId, :) = [2*x-1, 2*y-1, z, 2*x-2, 2*y-1, z+1]; edgeId = edgeId + 1;
                    edgeLattice(edgeId, :) = [2*x-1, 2*y-1, z, 2*x-1, 2*y-2, z+1]; edgeId = edgeId + 1;
                    edgeLattice(edgeId, :) = [2*x-1, 2*y-1, z, 2*x, 2*y-1, z+1]; edgeId = edgeId + 1;
                    edgeLattice(edgeId, :) = [2*x-1, 2*y-1, z, 2*x-1, 2*y, z+1]; edgeId = edgeId + 1;
                    edgeLattice(edgeId, :) = [2*x, 2*y, z, 2*x-1, 2*y, z+1]; edgeId = edgeId + 1;
                    edgeLattice(edgeId, :) = [2*x, 2*y, z, 2*x, 2*y-1, z+1]; edgeId = edgeId + 1;
                    edgeLattice(edgeId, :) = [2*x, 2*y, z, 2*x+1, 2*y, z+1]; edgeId = edgeId + 1;
                    edgeLattice(edgeId, :) = [2*x, 2*y, z, 2*x, 2*y+1, z+1]; edgeId = edgeId + 1;
                end
            end
        end
    else
        for x = 1:nx
            for y = 1:ny
                nodeLattice(nodeId, :) = [2*x-1, 2*y, z]; nodeId = nodeId + 1;
                nodeLattice(nodeId, :) = [2*x, 2*y-1, z]; nodeId = nodeId + 1;
                edgeLattice(edgeId, :) = [2*x-1, 2*y, z, 2*x-2, 2*y+1, z]; edgeId = edgeId + 1;
                edgeLattice(edgeId, :) = [2*x-1, 2*y, z, 2*x, 2*y+1, z]; edgeId = edgeId + 1;
                edgeLattice(edgeId, :) = [2*x, 2*y-1, z, 2*x-1, 2*y, z]; edgeId = edgeId + 1;
                edgeLattice(edgeId, :) = [2*x, 2*y-1, z, 2*x+1, 2*y, z]; edgeId = edgeId + 1;
                if z~=nz
                    edgeLattice(edgeId, :) = [2*x-1, 2*y, z, 2*x-2, 2*y, z+1]; edgeId = edgeId + 1;
                    edgeLattice(edgeId, :) = [2*x-1, 2*y, z, 2*x-1, 2*y-1, z+1]; edgeId = edgeId + 1;
                    edgeLattice(edgeId, :) = [2*x-1, 2*y, z, 2*x, 2*y, z+1]; edgeId = edgeId + 1;
                    edgeLattice(edgeId, :) = [2*x-1, 2*y, z, 2*x-1, 2*y+1, z+1]; edgeId = edgeId + 1;
                    edgeLattice(edgeId, :) = [2*x, 2*y-1, z, 2*x-1, 2*y-1, z+1]; edgeId = edgeId + 1;
                    edgeLattice(edgeId, :) = [2*x, 2*y-1, z, 2*x, 2*y-2, z+1]; edgeId = edgeId + 1;
                    edgeLattice(edgeId, :) = [2*x, 2*y-1, z, 2*x+1, 2*y-1, z+1]; edgeId = edgeId + 1;
                    edgeLattice(edgeId, :) = [2*x, 2*y-1, z, 2*x, 2*y, z+1]; edgeId = edgeId + 1;
                end
            end
        end
    end
end
assert(nodeId==(nodeNum+1));
assert(edgeId==(nodeNum*6-nx*ny*2*4+1));

% Introduce periodic boundary
if periodic
    edgeLattice(edgeLattice(:,4)==0,4) = 2*nx;
    edgeLattice(edgeLattice(:,5)==0,5) = 2*ny;
    edgeLattice(edgeLattice(:,4)==2*nx+1,4) = 1;
    edgeLattice(edgeLattice(:,5)==2*ny+1,5) = 1;
else
    edgeLattice = edgeLattice(edgeLattice(:,4)~=0,:);
    edgeLattice = edgeLattice(edgeLattice(:,5)~=0,:);
    edgeLattice = edgeLattice(edgeLattice(:,4)~=2*nx+1,:);
    edgeLattice = edgeLattice(edgeLattice(:,5)~=2*ny+1,:);
end

% Identify the nodes corresponding to each edge based on the lattice
% coordinates.
edgeLatticeFrom = edgeLattice(:, 1:3);
edgeLatticeTo = edgeLattice(:, 4:6);
[bool, edgeIndexFrom] = ismember(edgeLatticeFrom, nodeLattice, 'rows'); assert(all(bool));
[bool, edgeIndexTo] = ismember(edgeLatticeTo, nodeLattice, 'rows'); assert(all(bool));
edgeList = [edgeIndexFrom, edgeIndexTo];
edgeNum = size(edgeList,1);

% Determine the neighbours of each nodes
nodeConnectivity = zeros(nodeNum,1);
fprintf('Generating neighbour lists ');
switch directedOrUndirected
    case 'directed'
        nodeEdgeTo = zeros(nodeNum,6);
        nodeNodeTo = nodeEdgeTo;
        for i = 1:nodeNum
            idxEdges = find(edgeList(:,1)==i);
            idxNodes = edgeList(edgeList(:,1)==i,2)';
            nodeEdgeTo(i,1:length(idxEdges)) = idxEdges;
            nodeNodeTo(i,1:length(idxEdges)) = idxNodes;
            nodeConnectivity(i) = length(idxEdges);
            if i~=1
                fprintf(repmat('\b', 1, 7));
            end
            fprintf('%6.2f%%', i/nodeNum*100);
        end
    case 'undirected'
        nodeEdgeTo = zeros(nodeNum,12);
        nodeNodeTo = nodeEdgeTo;
        for i = 1:nodeNum
            idxEdges = [find(edgeList(:,1)==i);find(edgeList(:,2)==i)];
            idxNodes = [edgeList(edgeList(:,1)==i,2);edgeList(edgeList(:,2)==i,1)];
            nodeEdgeTo(i,1:length(idxEdges)) = idxEdges;
            nodeNodeTo(i,1:length(idxEdges)) = idxNodes;
            nodeConnectivity(i) = length(idxEdges);
            if i~=1
                fprintf(repmat('\b', 1, 7));
            end
            fprintf('%6.2f%%', i/nodeNum*100);
        end
end
fprintf('\n');

% Determine the inlets and the outlets
nodeInlet = zeros(nodeNum,1);
nodeOutlet = nodeInlet;
nodeInlet(1:(nx*ny*2)) = 1;
nodeOutlet((nz-1)*(nx*ny*2)+1:end) = 1;
end
function [nodeNum, edgeNum, edgeList, nodeEdgeTo, nodeNodeTo, nodeConnectivity, nodeLattice, nodeInlet, nodeOutlet] = generateSlantedCubicNetwork(nx, ny, nz, directedOrUndirected, periodic)
% This function generates a slanted cubic lattice

% Define the number of nodes, and the list of nodes (with one indexing,
% i.e. node index starts from 1). Loop through and assign a lattice
% coordinate to each node. Determine the lattice coordiantes of all edges
% connected to a node
nodeNum = nx*ny*nz*4;
nodeLattice = zeros(nodeNum, 3);
edgeLattice = zeros(nodeNum*3-nx*ny*4*3, 6);
nodeId = 1;
edgeId = 1;
for z = 1:nz
    if mod(z,3)==1
        for x = 1:nx
            for y = 1:ny
                nodeLattice(nodeId, :) = [4*x-3, 6*y-4, z]; nodeId = nodeId + 1;
                nodeLattice(nodeId, :) = [4*x-1, 6*y-4, z]; nodeId = nodeId + 1;
                nodeLattice(nodeId, :) = [4*x-2, 6*y-1, z]; nodeId = nodeId + 1;
                nodeLattice(nodeId, :) = [4*x, 6*y-1, z]; nodeId = nodeId + 1;
                if z~=nz
                    edgeLattice(edgeId, :) = [4*x-3, 6*y-4, z, 4*x-3, 6*y-2, z+1]; edgeId = edgeId + 1;
                    edgeLattice(edgeId, :) = [4*x-3, 6*y-4, z, 4*x-4, 6*y-5, z+1]; edgeId = edgeId + 1;
                    edgeLattice(edgeId, :) = [4*x-3, 6*y-4, z, 4*x-2, 6*y-5, z+1]; edgeId = edgeId + 1;
                    edgeLattice(edgeId, :) = [4*x-1, 6*y-4, z, 4*x-1, 6*y-2, z+1]; edgeId = edgeId + 1;
                    edgeLattice(edgeId, :) = [4*x-1, 6*y-4, z, 4*x-2, 6*y-5, z+1]; edgeId = edgeId + 1;
                    edgeLattice(edgeId, :) = [4*x-1, 6*y-4, z, 4*x, 6*y-5, z+1]; edgeId = edgeId + 1;
                    edgeLattice(edgeId, :) = [4*x-2, 6*y-1, z, 4*x-2, 6*y+1, z+1]; edgeId = edgeId + 1;
                    edgeLattice(edgeId, :) = [4*x-2, 6*y-1, z, 4*x-3, 6*y-2, z+1]; edgeId = edgeId + 1;
                    edgeLattice(edgeId, :) = [4*x-2, 6*y-1, z, 4*x-1, 6*y-2, z+1]; edgeId = edgeId + 1;
                    edgeLattice(edgeId, :) = [4*x, 6*y-1, z, 4*x, 6*y+1, z+1]; edgeId = edgeId + 1;
                    edgeLattice(edgeId, :) = [4*x, 6*y-1, z, 4*x-1, 6*y-2, z+1]; edgeId = edgeId + 1;
                    edgeLattice(edgeId, :) = [4*x, 6*y-1, z, 4*x+1, 6*y-2, z+1]; edgeId = edgeId + 1;
                end
            end
        end
    elseif mod(z,3)==2
        for x = 1:nx
            for y = 1:ny
                nodeLattice(nodeId, :) = [4*x-2, 6*y-5, z]; nodeId = nodeId + 1;
                nodeLattice(nodeId, :) = [4*x, 6*y-5, z]; nodeId = nodeId + 1;
                nodeLattice(nodeId, :) = [4*x-3, 6*y-2, z]; nodeId = nodeId + 1;
                nodeLattice(nodeId, :) = [4*x-1, 6*y-2, z]; nodeId = nodeId + 1;
                if z~=nz
                    edgeLattice(edgeId, :) = [4*x-2, 6*y-5, z, 4*x-2, 6*y-3, z+1]; edgeId = edgeId + 1;
                    edgeLattice(edgeId, :) = [4*x-2, 6*y-5, z, 4*x-3, 6*y-6, z+1]; edgeId = edgeId + 1;
                    edgeLattice(edgeId, :) = [4*x-2, 6*y-5, z, 4*x-1, 6*y-6, z+1]; edgeId = edgeId + 1;
                    edgeLattice(edgeId, :) = [4*x, 6*y-5, z, 4*x, 6*y-3, z+1]; edgeId = edgeId + 1;
                    edgeLattice(edgeId, :) = [4*x, 6*y-5, z, 4*x-1, 6*y-6, z+1]; edgeId = edgeId + 1;
                    edgeLattice(edgeId, :) = [4*x, 6*y-5, z, 4*x+1, 6*y-6, z+1]; edgeId = edgeId + 1;
                    edgeLattice(edgeId, :) = [4*x-3, 6*y-2, z, 4*x-3, 6*y, z+1]; edgeId = edgeId + 1;
                    edgeLattice(edgeId, :) = [4*x-3, 6*y-2, z, 4*x-4, 6*y-3, z+1]; edgeId = edgeId + 1;
                    edgeLattice(edgeId, :) = [4*x-3, 6*y-2, z, 4*x-2, 6*y-3, z+1]; edgeId = edgeId + 1;
                    edgeLattice(edgeId, :) = [4*x-1, 6*y-2, z, 4*x-1, 6*y, z+1]; edgeId = edgeId + 1;
                    edgeLattice(edgeId, :) = [4*x-1, 6*y-2, z, 4*x-2, 6*y-3, z+1]; edgeId = edgeId + 1;
                    edgeLattice(edgeId, :) = [4*x-1, 6*y-2, z, 4*x, 6*y-3, z+1]; edgeId = edgeId + 1;
                end
            end
        end
    else
        for x = 1:nx
            for y = 1:ny
                nodeLattice(nodeId, :) = [4*x-2, 6*y-3, z]; nodeId = nodeId + 1;
                nodeLattice(nodeId, :) = [4*x, 6*y-3, z]; nodeId = nodeId + 1;
                nodeLattice(nodeId, :) = [4*x-3, 6*y, z]; nodeId = nodeId + 1;
                nodeLattice(nodeId, :) = [4*x-1, 6*y, z]; nodeId = nodeId + 1;
                if z~=nz
                    edgeLattice(edgeId, :) = [4*x-2, 6*y-3, z, 4*x-2, 6*y-1, z+1]; edgeId = edgeId + 1;
                    edgeLattice(edgeId, :) = [4*x-2, 6*y-3, z, 4*x-3, 6*y-4, z+1]; edgeId = edgeId + 1;
                    edgeLattice(edgeId, :) = [4*x-2, 6*y-3, z, 4*x-1, 6*y-4, z+1]; edgeId = edgeId + 1;
                    edgeLattice(edgeId, :) = [4*x, 6*y-3, z, 4*x, 6*y-1, z+1]; edgeId = edgeId + 1;
                    edgeLattice(edgeId, :) = [4*x, 6*y-3, z, 4*x-1, 6*y-4, z+1]; edgeId = edgeId + 1;
                    edgeLattice(edgeId, :) = [4*x, 6*y-3, z, 4*x+1, 6*y-4, z+1]; edgeId = edgeId + 1;
                    edgeLattice(edgeId, :) = [4*x-3, 6*y, z, 4*x-3, 6*y+2, z+1]; edgeId = edgeId + 1;
                    edgeLattice(edgeId, :) = [4*x-3, 6*y, z, 4*x-4, 6*y-1, z+1]; edgeId = edgeId + 1;
                    edgeLattice(edgeId, :) = [4*x-3, 6*y, z, 4*x-2, 6*y-1, z+1]; edgeId = edgeId + 1;
                    edgeLattice(edgeId, :) = [4*x-1, 6*y, z, 4*x-1, 6*y+2, z+1]; edgeId = edgeId + 1;
                    edgeLattice(edgeId, :) = [4*x-1, 6*y, z, 4*x-2, 6*y-1, z+1]; edgeId = edgeId + 1;
                    edgeLattice(edgeId, :) = [4*x-1, 6*y, z, 4*x, 6*y-1, z+1]; edgeId = edgeId + 1;
                end
            end
        end
    end
end
assert(nodeId==(nodeNum+1))
assert(edgeId==(nodeNum*3-nx*ny*4*3+1))

% Introduce periodic boundary
if periodic
    edgeLattice(edgeLattice(:,4)==0,4) = 4*nx;
    edgeLattice(edgeLattice(:,5)==-1,5) = 6*ny-1;
    edgeLattice(edgeLattice(:,5)==0,5) = 6*ny;
    edgeLattice(edgeLattice(:,4)==4*nx+1,4) = 1;
    edgeLattice(edgeLattice(:,5)==6*ny+1,5) = 1;
    edgeLattice(edgeLattice(:,5)==6*ny+2,5) = 2;
else
    edgeLattice = edgeLattice(edgeLattice(:,4)~=0,:);
    edgeLattice = edgeLattice(edgeLattice(:,5)~=-1,:);
    edgeLattice = edgeLattice(edgeLattice(:,5)~=0,:);
    edgeLattice = edgeLattice(edgeLattice(:,4)~=4*nx+1,:);
    edgeLattice = edgeLattice(edgeLattice(:,5)~=6*ny+1,:);
    edgeLattice = edgeLattice(edgeLattice(:,5)~=6*ny+2,:);
end

% Identify the nodes corresponding to each edge based on the lattice
% coordinates.
edgeLatticeFrom = edgeLattice(:, 1:3);
edgeLatticeTo = edgeLattice(:, 4:6);
[bool, edgeIndexFrom] = ismember(edgeLatticeFrom, nodeLattice, 'rows'); assert(all(bool));
[bool, edgeIndexTo] = ismember(edgeLatticeTo, nodeLattice, 'rows'); assert(all(bool));
edgeList = [edgeIndexFrom, edgeIndexTo];
edgeNum = size(edgeList,1);

% Determine the neighbours of each nodes
nodeConnectivity = zeros(nodeNum,1);
fprintf('Generating neighbour lists ');
switch directedOrUndirected
    case 'directed'
        nodeEdgeTo = zeros(nodeNum,3);
        nodeNodeTo = nodeEdgeTo;
        for i = 1:nodeNum
            idxEdges = find(edgeList(:,1)==i);
            idxNodes = edgeList(edgeList(:,1)==i,2)';
            nodeEdgeTo(i,1:length(idxEdges)) = idxEdges;
            nodeNodeTo(i,1:length(idxEdges)) = idxNodes;
            nodeConnectivity(i) = length(idxEdges);
            if i~=1
                fprintf(repmat('\b', 1, 7));
            end
            fprintf('%6.2f%%', i/nodeNum*100);
        end
    case 'undirected'
        nodeEdgeTo = zeros(nodeNum,6);
        nodeNodeTo = nodeEdgeTo;
        for i = 1:nodeNum
            idxEdges = [find(edgeList(:,1)==i);find(edgeList(:,2)==i)];
            idxNodes = [edgeList(edgeList(:,1)==i,2);edgeList(edgeList(:,2)==i,1)];
            nodeEdgeTo(i,1:length(idxEdges)) = idxEdges;
            nodeNodeTo(i,1:length(idxEdges)) = idxNodes;
            nodeConnectivity(i) = length(idxEdges);
            if i~=1
                fprintf(repmat('\b', 1, 7));
            end
            fprintf('%6.2f%%', i/nodeNum*100);
        end
end
fprintf('\n');

% Determine the inlets and the outlets
nodeInlet = zeros(nodeNum,1);
nodeOutlet = nodeInlet;
nodeInlet(1:(4*nx*ny)) = 1;
nodeOutlet((nz-1)*(4*nx*ny)+1:end) = 1;
end
function [nodeNum, edgeNum, edgeList, nodeEdgeTo, nodeNodeTo, nodeConnectivity, nodeLattice, nodeInlet, nodeOutlet] = generateSlantedSquareNetwork(nx, nz, directedOrUndirected, periodic)
% This function generates a slanted square lattice

% Define the number of nodes, and the list of nodes (with one indexing,
% i.e. node index starts from 1). Loop through and assign a lattice
% coordinate to each node. Determine the lattice coordiantes of all edges
% connected to a node
nodeNum = nx*nz;
nodeLattice = zeros(nodeNum, 3);
edgeLattice = zeros(nodeNum*2-nx*2, 6);
nodeId = 1;
edgeId = 1;
for z = 1:nz
    if mod(z,2)==1
        for x = 1:nx
            nodeLattice(nodeId, :) = [2*x-1, 1, z]; nodeId = nodeId + 1;
            if z~=nz
                edgeLattice(edgeId, :) = [2*x-1, 1, z, 2*x, 1, z+1]; edgeId = edgeId + 1;
                edgeLattice(edgeId, :) = [2*x-1, 1, z, 2*x-2, 1, z+1]; edgeId = edgeId + 1;
            end
        end
    else
        for x = 1:nx
            nodeLattice(nodeId, :) = [2*x, 1, z]; nodeId = nodeId + 1;
            if z~=nz
                edgeLattice(edgeId, :) = [2*x, 1, z, 2*x+1, 1, z+1]; edgeId = edgeId + 1;
                edgeLattice(edgeId, :) = [2*x, 1, z, 2*x-1, 1, z+1]; edgeId = edgeId + 1;
            end
        end
    end
end
assert(nodeId==(nodeNum+1))
assert(edgeId==(nodeNum*2-nx*2+1))

% Introduce periodic boundary
if periodic
    edgeLattice(edgeLattice(:,4)==0,4) = 2*nx;
    edgeLattice(edgeLattice(:,4)==2*nx+1,4) = 1;
else
    edgeLattice = edgeLattice(edgeLattice(:,4)~=0,:);
    edgeLattice = edgeLattice(edgeLattice(:,4)~=2*nx+1,:);
end

% Identify the nodes corresponding to each edge based on the lattice
% coordinates.
edgeLatticeFrom = edgeLattice(:, 1:3);
edgeLatticeTo = edgeLattice(:, 4:6);
[bool, edgeIndexFrom] = ismember(edgeLatticeFrom, nodeLattice, 'rows'); assert(all(bool));
[bool, edgeIndexTo] = ismember(edgeLatticeTo, nodeLattice, 'rows'); assert(all(bool));
edgeList = [edgeIndexFrom, edgeIndexTo];
edgeNum = size(edgeList,1);

% Determine the neighbours of each nodes
nodeConnectivity = zeros(nodeNum,1);
fprintf('Generating neighbour lists ');
switch directedOrUndirected
    case 'directed'
        nodeEdgeTo = zeros(nodeNum,2);
        nodeNodeTo = nodeEdgeTo;
        for i = 1:nodeNum
            idxEdges = find(edgeList(:,1)==i);
            idxNodes = edgeList(edgeList(:,1)==i,2)';
            nodeEdgeTo(i,1:length(idxEdges)) = idxEdges;
            nodeNodeTo(i,1:length(idxEdges)) = idxNodes;
            nodeConnectivity(i) = length(idxEdges);
            if i~=1
                fprintf(repmat('\b', 1, 7));
            end
            fprintf('%6.2f%%', i/nodeNum*100);
        end
    case 'undirected'
        nodeEdgeTo = zeros(nodeNum,4);
        nodeNodeTo = nodeEdgeTo;
        for i = 1:nodeNum
            idxEdges = [find(edgeList(:,1)==i);find(edgeList(:,2)==i)];
            idxNodes = [edgeList(edgeList(:,1)==i,2);edgeList(edgeList(:,2)==i,1)];
            nodeEdgeTo(i,1:length(idxEdges)) = idxEdges;
            nodeNodeTo(i,1:length(idxEdges)) = idxNodes;
            nodeConnectivity(i) = length(idxEdges);
            if i~=1
                fprintf(repmat('\b', 1, 7));
            end
            fprintf('%6.2f%%', i/nodeNum*100);
        end
end
fprintf('\n');

% Determine the inlets and the outlets
nodeInlet = zeros(nodeNum,1);
nodeOutlet = nodeInlet;
nodeInlet(1:nx) = 1;
nodeOutlet((nz-1)*nx+1:end) = 1;
end