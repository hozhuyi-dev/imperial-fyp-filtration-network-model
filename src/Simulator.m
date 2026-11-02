classdef Simulator
    % The top-level class responsible for simulating filtration processes
    % using network models
    % 
    % =====================================================================
    % AVAILABLE OPTIONS: 
    % 
    % introductionMethod (Affect private function introduceFinesParticle)
    % (1) refillInlets (Default) - Refill the inlets to the maximum volume,
    %     i.e. 0.5*nodeVolume by default, every time step until t reaches
    %     maxFinesTimeStep.
    % (2) fixedParticleNum - Refill the inlets to the maximum volume until
    %     the number of particles reaches storeNum. 
    % (3) fixedParticleNum2 (Not used) - Introduce a fixed number of
    %     particles every time step (equal to inletNum*particleInlet) until
    %     t reaches maxFinesTimeStep. If the maximum volume is exceeded,
    %     the particle is store in timeOutside. 
    % 
    % transportLogicType (Affect private function transportLogic)
    % (1) flowbiased (Default) - The particles 'choose' the edge to follow
    %     based on flow-biased probabilities. Particles assumed to be in
    %     suspension with fluid
    % (2) flowbiasedpos - The particles still obey flow-biased
    %     probabilities, but an additional restraint is applied such the
    %     particles cannot flow in negative z-direction. 
    % (3) flowbiasednode - Flow-biased probabilities, and a node is
    %     considered clogged once a particle is retained in the node. 
    % (4) flowbiasedextreme - Flow-biased probabilities with conditions (2)
    %     and (3) as well. 
    % (5) noclogging - Particles follow the edge that (i) meet edgeDiameter 
    %     and nodeVolume constraints and (ii) has the highest flow. 
    % (6) noflow - Fluid flows are not calculated, edges flow based on
    %     area-biased probabilities as used in Shire & O'Sullivan (2017).
    %     Total volume of particles in each node cannot exceed the limit.
    % (7) noflowrandom - Instead of area-biased probabilies, the particles
    %     chooses the exits with equal probabilities. 
    % 
    % storeOption
    % (1) normal (Default)
    % (2) storePath - Store the paths of all particles; memory intensive.
    %     Suggest to limit the number of particles with fixedParticleNum or
    %     use only with effective filters such that less particles are
    %     simulated.
    % (3) more - Store timeNodeRetained (see Statistics class for details), 
    %     very large data set
    % (4) storeEdgePath - Store the number of eroded particles passing
    %     through each edge. 
    % (5) storeNodeClogging - Store nodeClogging (see Statistics class for
    %     details). 
    %
    % flowMethod
    % (1) matlab (Default) - Solve with MATLAB using Gaussion elimination
    %     with sparse matrix
    % (2) python - Solve with OpenPNM (Gostick et. al. 2016); slower than
    %     with matlab
    % =====================================================================
    % PARAMETERS
    %
    % densityLimit (Default is 0.5) - A node is considered full when
    %     the total volume of base particles within the node reaches 
    %     densityLimit*nodeVolume
    % 
    % contentLimit (Default is 1) - Only affects introduction of particles.
    %     Inlets are considered full when
    %     contentLimit*densityLimit*nodeVolume is reached
    % 
    % maxFinesParticle (Default is 5e5) - When the number of simultaneous
    %     reaches 5e5, the introduction of base particles is stopped, until
    %     the total number reduces to below 5e5
    %
    % maxFinesTimeStep (Default is 5*nz for equivalent cubic network) -
    %     Maximum number of time steps to introduce base particles
    %
    % maxTimeStep (Default is 5*maxFinesTimeStep) - Stop simulations
    %     regardless of whether the number of base particles has stabilised
    %     or reduced to zero.
    % 
    % maxEndSamePath (Default is 5) - End simulation if particles alternate
    %     between the same edge for maxEndSamePath number of time steps.
    %
    % maxEndSameNumber (Default is 100) - End simulations if the number of
    %     particles stabilised, i.e. remain the same, for maxEndSameNumber
    %     number of time steps.
    %
    % storeNum - When fixedParticleNum (introductionMethod) is used,
    %     storeNum defines the number of particles to introduce in total. 
    % 
    % particleInlet - When refillInlets or fixedParticleNum2
    %     (introductionMethod) is used, defines the maximum number of
    %     particles to introduce each time step as inletNum*particleInlet
    % =====================================================================
    % FLOW SOLVER PARAMETERS
    %
    % pressureEntry (Default is 1) - Pressure imposed at the inlets. 
    %
    % pressureExit (Default is 0) - Pressure imposed at the outlets.
    % =====================================================================
    
    properties
        % Filter-base combination
        filter
        base = struct;
        
        % Parameters
        simulatorParameters = struct('introductionMethod', [], ...
            'transportLogicType', [], ...
            'storeOption', [], ...
            'densityLimit', [], ...
            'contentLimit', [], ...
            'maxFinesParticle', [], ...
            'maxFinesTimeStep', [], ...
            'maxTimeStep', [], ...
            'maxEndSamePath', [], ...
            'maxEndSameNumber', [], ...
            'storeNum', [], ...
            'particleInlet', []);
        flowSolverParameters = struct('pressureEntry', [], ...
            'pressureExit', [], ...
            'flowMethod', []);
        
        % Results
        statistics = struct;
        
        % Save location (if saved or loaded)
        fileLocation
    end
    
    methods
        % ===== Initialisation ============================================
        % Initialise Simulator object, choose options through this function
        function [obj] = Simulator(varargin)
            % Creates a Simulator object, while defining simulation
            % and flow solver properties
            parser = inputParser;
            expectedSimulationMethod = Constants.introductionMethod;
            addParameter(parser, 'introductionMethod', 'refillInlet', @(x) any(validatestring(x, expectedSimulationMethod)))
            expectedTransportLogicType = Constants.transportLogicType;
            addParameter(parser, 'transportLogicType', 'flowbiased', @(x) any(validatestring(x, expectedTransportLogicType)))
            expectedStoreOption = Constants.storeOption;
            addParameter(parser, 'storeOption', 'normal', @(x) any(validatestring(x, expectedStoreOption)))
            addParameter(parser, 'densityLimit', 0.5)
            addParameter(parser, 'contentLimit', [])
            addParameter(parser, 'maxFinesParticle', 5e5)
            addParameter(parser, 'maxFinesTimeStep', [])
            expectedFlowMethod = {'matlab', 'python'};
            addParameter(parser, 'pressureEntry', 1)
            addParameter(parser, 'pressureExit', 0)
            addParameter(parser, 'flowMethod', 'matlab', @(x) any(validatestring(x, expectedFlowMethod)))
            addParameter(parser, 'maxTimeStep', [])
            addParameter(parser, 'maxEndSamePath', 5)
            addParameter(parser, 'maxEndSameNumber', 100)
            addParameter(parser, 'storeNum', [])
            addParameter(parser, 'particleInlet', [])
            parse(parser, varargin{:})
            obj.simulatorParameters.introductionMethod = parser.Results.introductionMethod;
            obj.simulatorParameters.transportLogicType = parser.Results.transportLogicType;
            obj.simulatorParameters.storeOption = parser.Results.storeOption;
            obj.simulatorParameters.densityLimit = parser.Results.densityLimit;
            obj.simulatorParameters.contentLimit = parser.Results.contentLimit;
            obj.simulatorParameters.maxFinesParticle = parser.Results.maxFinesParticle;
            obj.simulatorParameters.maxFinesTimeStep = parser.Results.maxFinesTimeStep;
            obj.simulatorParameters.maxTimeStep = parser.Results.maxTimeStep;
            obj.simulatorParameters.maxEndSamePath = parser.Results.maxEndSamePath;
            obj.simulatorParameters.maxEndSameNumber = parser.Results.maxEndSameNumber;
            obj.simulatorParameters.storeNum = parser.Results.storeNum;
            obj.simulatorParameters.particleInlet = parser.Results.particleInlet;
            obj.flowSolverParameters.pressureEntry = parser.Results.pressureEntry;
            obj.flowSolverParameters.pressureExit = parser.Results.pressureExit;
            obj.flowSolverParameters.flowMethod = parser.Results.flowMethod;
        end
        
        % Creating associated Filter object
        function [obj] = createFilter(obj, cuCoarse, varargin)
            % Validate input
            expectedCuCoarse = Constants.cuCoarse;
            if ~any(expectedCuCoarse==cuCoarse)
                errorMessage = sprintf('Error: The following cuCoarse does not exist: %.1f\n', cuCoarse);
                fprintf(errorMessage);
                return;
            end
            
            % Input parser
            parser = inputParser;
            addParameter(parser, 'minCoarse', 1)
            expectedDensity = Constants.density;
            addParameter(parser, 'density', 'dense', @(x) any(validatestring(x, expectedDensity)))
            expectedSize = Constants.size;
            addParameter(parser, 'size', [1,1,1], @(x) any(all(expectedSize==x, 2)))
            expectedTriangulationMethod = Constants.triangulationMethod;
            addParameter(parser, 'triangulationMethod', 'Del', @(x) any(validatestring(x, expectedTriangulationMethod)))
            expectedLoadMethod = Constants.loadMethod;
            addParameter(parser, 'loadMethod', 'read', @(x) any(validatestring(x, expectedLoadMethod)))
            expectedNetworkType = Constants.networkType;
            expectedAssignmentMethod = Constants.assignmentMethod;
            addParameter(parser, 'networkType', 'actual', @(x) any(validatestring(x, expectedNetworkType)))
            addParameter(parser, 'assignmentMethod', 'actual', @(x) any(validatestring(x, expectedAssignmentMethod)))
            expectedDirectedOrUndirected = Constants.directedOrUndirected;
            expectedFlowDirection = Constants.flowDirection;
            addParameter(parser, 'directedOrUndirected', 'undirected', @(x) any(validatestring(x, expectedDirectedOrUndirected)))
            addParameter(parser, 'nxyz', [])
            addParameter(parser, 'flowDirection', 'zpos', @(x) any(validatestring(x, expectedFlowDirection)))
            parse(parser, varargin{:})
            minCoarse = parser.Results.minCoarse;
            density = parser.Results.density;
            size = parser.Results.size;
            triangulationMethod = parser.Results.triangulationMethod;
            loadMethod = parser.Results.loadMethod;
            networkType = parser.Results.networkType;
            if strcmp(networkType, 'actual')
                assignmentMethod = parser.Results.assignmentMethod;
            else
                assignmentMethod = 'random';
            end
            directedOrUndirected = parser.Results.directedOrUndirected;
            nxyz = parser.Results.nxyz;
            flowDirection = parser.Results.flowDirection;
            
            % Determine nxyz such the nodeNum is comparable
            if and(isempty(nxyz), ~strcmp(networkType, 'actual'))
                % Determine the number of nodes
                tempFilter = Filter(1,cuCoarse,density,size,triangulationMethod, ...
                    'actual','actual',[],flowDirection, directedOrUndirected);
                tempFilter = tempFilter.loadDTM('read');
                nxyz = Constants.comparableNxyz(tempFilter.nodeNum,networkType,size);
            end
            
            % Create filter
            obj.filter = Filter(minCoarse, cuCoarse, density, size, triangulationMethod, ...
                networkType, assignmentMethod, nxyz, flowDirection, directedOrUndirected);
            obj.filter = obj.filter.loadDTM(loadMethod);
            obj.filter = obj.filter.generateNetwork();
            
            % Number of time steps to introduce fines
            tempNxyz = Constants.comparableNxyz(obj.filter.nodeNum, 'cubic', size);
            switch obj.simulatorParameters.introductionMethod
                case 'refillInlet'
                    if isempty(obj.simulatorParameters.maxFinesTimeStep)
                        obj.simulatorParameters.maxFinesTimeStep = tempNxyz(3)*5;
                    end
                    if isempty(obj.simulatorParameters.particleInlet)
                        obj.simulatorParameters.particleInlet = 100;
                    end
                    if isempty(obj.simulatorParameters.contentLimit)
                        obj.simulatorParameters.contentLimit = 1;
                    end
                case 'fixedParticleNum'
                    if isempty(obj.simulatorParameters.maxFinesTimeStep)
                        obj.simulatorParameters.maxFinesTimeStep = tempNxyz(3)*5;
                    end
                    if isempty(obj.simulatorParameters.storeNum)
                        obj.simulatorParameters.storeNum = 1e4;
                    end
                    if isempty(obj.simulatorParameters.contentLimit)
                        obj.simulatorParameters.contentLimit = 1;
                    end
                case 'fixedParticleNum2'
                    if isempty(obj.simulatorParameters.maxFinesTimeStep)
                        obj.simulatorParameters.maxFinesTimeStep = 100;
                    end
                    if isempty(obj.simulatorParameters.particleInlet)
                        obj.simulatorParameters.particleInlet = 1;
                    end
                    assert(obj.simulatorParameters.particleInlet*obj.simulatorParameters.maxFinesTimeStep*sum(obj.filter.nodeInlet)<obj.simulatorParameters.maxFinesParticle);
            end
            if isempty(obj.simulatorParameters.maxTimeStep)
                obj.simulatorParameters.maxTimeStep = tempNxyz(3)*25;
            end
        end
        
        % Creating Base objects
        function [obj] = createBase(obj, d50, varargin)
            % Input parser
            parser = inputParser;
            addParameter(parser, 'cuFine', 1.2)
            addParameter(parser, 'run', 1)
            parse(parser, varargin{:});
            cuFine = parser.Results.cuFine;
            run = parser.Results.run;
            
            % Create base
            obj.base.(sprintf('base%i',run)) = Base();
            obj.base.(sprintf('base%i',run)) = obj.base.(sprintf('base%i',run)).linearVolBase(d50, cuFine);
        end
        
        % ===== Simulation ================================================
        % Main function for the simulation of filtration
        function [obj] = simulateFiltration(obj, varargin)
            % The function simulates the filtration process; main steps are
            % given here:
            % (1) Base particles are introduced at the inlets
            % (2) Particles flow through the network of the filter
            % (3) Until they are retained or eroded
            
            % Inputs
            parser = inputParser;
            addParameter(parser, 'followingRun', 1);
            addParameter(parser, 'run', 1);
            parse(parser, varargin{:});
            followingRun  = parser.Results.followingRun;
            run = parser.Results.run;
            
            % Unpack variables to reduce the size of variables passed
            % around the different functions
            [flowDirection, networkType] = obj.unpack(obj.filter.networkParameters, 'flowDirection', 'networkType');
            [nodeNum, nodeCoord, nodeInlet, nodeLattice, nodeOutlet, nodeSurfaceArea, nodeVolume] = obj.unpack(obj.filter, 'nodeNum', 'nodeCoord', 'nodeInlet', 'nodeLattice', 'nodeOutlet', 'nodeSurfaceArea', 'nodeVolume');
            [nodeConnectivity, nodeEdgeTo, nodeNodeTo, nodeEdgeToPos, nodeEdgeDirection] = obj.unpack(obj.filter, 'nodeConnectivity', 'nodeEdgeTo', 'nodeNodeTo', 'nodeEdgeToPos', 'nodeEdgeDirection');
            [edgeNum, edgeArea, edgeDiameter, edgeLength, edgeList, edgeReduced] = obj.unpack(obj.filter, 'edgeNum', 'edgeArea', 'edgeDiameter', 'edgeLength', 'edgeList', 'edgeReduced');
            [psdFineNum, psdFineBucket] = obj.unpack(obj.base.(sprintf('base%i',run)), 'psdFineNum', 'psdFineBucket');
            [introductionMethod, densityLimit, contentLimit, maxFinesParticle, maxFinesTimeStep, transportLogicType, maxTimeStep, maxEndSamePath, maxEndSameNumber, storeNum, particleInlet, storeOption] = obj.unpack(obj.simulatorParameters, ...
                'introductionMethod', 'densityLimit', 'contentLimit', 'maxFinesParticle', 'maxFinesTimeStep', 'transportLogicType', 'maxTimeStep', 'maxEndSamePath', 'maxEndSameNumber', 'storeNum', 'particleInlet', 'storeOption');
            [pressureEntry, pressureExit, flowMethod] = obj.unpack(obj.flowSolverParameters, 'pressureEntry', 'pressureExit', 'flowMethod');
            
            % Create a Statistics object to store the output
            obj.statistics.(sprintf('statistics%i', run)) = Statistics(nodeNum, edgeNum, nodeLattice, psdFineBucket, maxTimeStep);
            [edgePath, particlePathEroded, particlePathRetained, timeEroded, timeIntroduced, timeRetained, timeOutside, timeTotalFlow] = obj.unpack(obj.statistics.(sprintf('statistics%i', run)), ...
                'edgePath', 'particlePathEroded', 'particlePathRetained', 'timeEroded', 'timeIntroduced', 'timeRetained', 'timeOutside', 'timeTotalFlow');
            [particleDistanceRetained, particleDistanceEroded, edgePathEroded, edgePathRetained] = obj.unpack(obj.statistics.(sprintf('statistics%i', run)), ...
                'particleDistanceRetained', 'particleDistanceEroded', 'edgePathEroded', 'edgePathRetained');
            [nodeClogging, nodeRetained, edgeFlowing, layerRetained] = obj.unpack(obj.statistics.(sprintf('statistics%i',followingRun)), 'nodeClogging', 'nodeRetained', 'edgeFlowing', 'layerRetained');
            [particleRoute, timeNodeRetained] = obj.unpack(obj.statistics.(sprintf('statistics%i', run)), 'particleRoute', 'timeNodeRetained');
            
            % Initialisation
            finesNum = 1;
            t = 1;
            isFull = 0;
            [nodeActive, nodeClogging, F, idxParticle] = initialiseArrays(nodeNum, nodeClogging);
            
            % Update nodeActive
            [nodeActive] = updateNodeActive(nodeNum, nodeActive, nodeRetained, F);
            
            if strcmp(transportLogicType(1:6), 'noflow')
                % Initialise
                nodePressure = NaN;
                edgeConductance = NaN;
                
                % Main time loop
                while t<=maxTimeStep && (finesNum~=0 || t<=maxFinesTimeStep)
                    % Introduce fines
                    [nodeActive, F, timeIntroduced, timeOutside, particleRoute, idxParticle, maxFinesTimeStep] = introduceFinesParticle(edgeNum, nodeActive, nodeInlet, nodeVolume, F, timeIntroduced, timeOutside, particleRoute, idxParticle, psdFineNum, psdFineBucket, introductionMethod, densityLimit, contentLimit, maxFinesParticle, maxFinesTimeStep, storeOption, storeNum, particleInlet, t);
                    
                    % Update nodeActive
                    [nodeActive] = updateNodeActive(nodeNum, nodeActive, nodeRetained, F);
                    
                    % Perform transport logic
                    [nodeActive, nodeClogging, edgeFlowing, edgePath, particleRoute, F, finesNum] = transportLogic(nodeActive, nodeClogging, nodeConnectivity, nodeEdgeTo, nodeNodeTo, nodeEdgeToPos, nodePressure, nodeOutlet, nodeRetained, nodeVolume, edgeArea, edgeConductance, edgeDiameter, edgeFlowing, edgeLength, edgePath, particleRoute, F, transportLogicType, densityLimit, storeOption, t);
                    
                    % Memory management
                    [F, timeEroded, timeRetained, timeNodeRetained, nodeRetained, layerRetained, particlePathEroded, particlePathRetained, particleDistanceEroded, particleDistanceRetained, particleRoute, idxParticle, edgePathEroded, edgePathRetained] = rearrangeF(nodeNum, edgeNum, nodeConnectivity, nodeEdgeTo, nodeEdgeDirection, nodeInlet, nodeLattice, nodePressure, edgeConductance, edgeFlowing, edgeList, edgePathEroded, edgePathRetained, F, psdFineBucket, timeEroded, timeRetained, timeNodeRetained, nodeRetained, layerRetained, particlePathEroded, particlePathRetained, particleDistanceEroded, particleDistanceRetained, particleRoute, idxParticle, storeOption, t);

                    % Move time counter forward
                    t = t+1;
                end
            else
                % Initialise counter
                counter1 = 0;
                counter2 = 0;
                
                % Main time loop
                while t<=maxTimeStep && ~isFull && (finesNum~=0 || t<=maxFinesTimeStep)
                    % Introduce fines
                    [nodeActive, F, timeIntroduced, timeOutside, particleRoute, idxParticle, maxFinesTimeStep] = introduceFinesParticle(edgeNum, nodeActive, nodeInlet, nodeVolume, F, timeIntroduced, timeOutside, particleRoute, idxParticle, psdFineNum, psdFineBucket, introductionMethod, densityLimit, contentLimit, maxFinesParticle, maxFinesTimeStep, storeOption, storeNum, particleInlet, t);
                    
                    % Update nodeActive
                    [nodeActive] = updateNodeActive(nodeNum, nodeActive, nodeRetained, F);
                    
                    % Flow solver
                    [nodePressure, nodeSource, edgeConductance] = obj.determinePressure(nodeNum, edgeNum, nodeActive, nodeConnectivity, nodeCoord, nodeEdgeTo, nodeInlet, nodeNodeTo, nodeOutlet, nodeRetained, nodeSurfaceArea, nodeVolume, edgeArea, edgeFlowing, edgeLength, edgeList, edgeReduced, pressureEntry, pressureExit, flowMethod, flowDirection, networkType);
                    timeTotalFlow(t) = sum(nodeSource(nodeInlet==1));
                    
                    % Perform transport logic
                    edgePath1 = edgePath;
                    [nodeActive, nodeClogging, edgeFlowing, edgePath, particleRoute, F, finesNum] = transportLogic(nodeActive, nodeClogging, nodeConnectivity, nodeEdgeTo, nodeNodeTo, nodeEdgeToPos, nodePressure, nodeOutlet, nodeRetained, nodeVolume, edgeArea, edgeConductance, edgeDiameter, edgeFlowing, edgeLength, edgePath, particleRoute, F, transportLogicType, densityLimit, storeOption, t);
                    edgePath2 = edgePath;
                    edgePathChange = edgePath2-edgePath1;
                    clear edgepath1 edgePath2;
                    
                    if t>maxFinesTimeStep
                        % Increase counter if particles move through the same edges
                        if all(edgePathChangePrevious==edgePathChange)
                            counter1 = counter1+1;
                            if counter1==maxEndSamePath
                                F(:,5) = 2;
                                finesNum = 0;
                                isFull = 1;
                            end
                        else
                            counter1 = 0;
                        end
                        
                        % Increase counter if number of particles stabilises
                        if finesNumPrevious==finesNum
                            counter2 = counter2+1;
                            if counter2==maxEndSameNumber
                                F(:,5) = 2;
                                finesNum = 0;
                            end
                        else
                            counter2 = 0;
                        end
                    end
                    edgePathChangePrevious = edgePathChange;
                    finesNumPrevious = finesNum;
                    
                    % Memory management
                    [F, timeEroded, timeRetained, timeNodeRetained, nodeRetained, layerRetained, particlePathEroded, particlePathRetained, particleDistanceEroded, particleDistanceRetained, particleRoute, idxParticle, edgePathEroded, edgePathRetained] = rearrangeF(nodeNum, edgeNum, nodeConnectivity, nodeEdgeTo, nodeEdgeDirection, nodeInlet, nodeLattice, nodePressure, edgeConductance, edgeFlowing, edgeList, edgePathEroded, edgePathRetained, F, psdFineBucket, timeEroded, timeRetained, timeNodeRetained, nodeRetained, layerRetained, particlePathEroded, particlePathRetained, particleDistanceEroded, particleDistanceRetained, particleRoute, idxParticle, storeOption, t);
                    
                    % Move time counter forward
                    t = t+1;
                end
            end
            
            % Finalise arrays
            [timeRetained, timeEroded, timeIntroduced, timeOutside, timeTotalFlow] = finaliseArrays(timeRetained, timeEroded, timeIntroduced, timeOutside, timeTotalFlow, maxFinesTimeStep, t);
            
            % Pack into struct
            obj.statistics.(sprintf('statistics%i', run)).nodeClogging = nodeClogging;
            obj.statistics.(sprintf('statistics%i', run)).nodeRetained = nodeRetained;
            obj.statistics.(sprintf('statistics%i', run)).edgeFlowing = edgeFlowing;
            obj.statistics.(sprintf('statistics%i', run)).edgePath = edgePath;
            obj.statistics.(sprintf('statistics%i', run)).edgePathEroded = edgePathEroded;
            obj.statistics.(sprintf('statistics%i', run)).edgePathRetained = edgePathRetained;
            obj.statistics.(sprintf('statistics%i', run)).particlePathEroded = particlePathEroded;
            obj.statistics.(sprintf('statistics%i', run)).particlePathRetained = particlePathRetained;
            obj.statistics.(sprintf('statistics%i', run)).particleDistanceEroded = particleDistanceEroded;
            obj.statistics.(sprintf('statistics%i', run)).particleDistanceRetained = particleDistanceRetained;
            obj.statistics.(sprintf('statistics%i', run)).particleRoute = particleRoute;
            obj.statistics.(sprintf('statistics%i', run)).layerRetained = layerRetained;
            obj.statistics.(sprintf('statistics%i', run)).timeRetained = timeRetained;
            obj.statistics.(sprintf('statistics%i', run)).timeEroded = timeEroded;
            obj.statistics.(sprintf('statistics%i', run)).timeIntroduced = timeIntroduced;
            obj.statistics.(sprintf('statistics%i', run)).timeOutside = timeOutside;
            obj.statistics.(sprintf('statistics%i', run)).timeTotalFlow = timeTotalFlow;
            obj.statistics.(sprintf('statistics%i', run)).timeNodeRetained = timeNodeRetained;
        end
        
        % ===== Packaging =================================================
        % Packaging of monodisperse simulations
        function [obj] = simulateMonodisperse(obj, d, varargin)
            % Input parser
            parser = inputParser;
            expectedParticleNum = {'normal','less'};
            addParameter(parser, 'particleNum', 'normal', @(x) any(validatestring(x, expectedParticleNum)))
            addParameter(parser, 'continuous', 0)
            addParameter(parser, 'endWith', '')
            addParameter(parser, 'fixedParticleNum', 0)
            addParameter(parser, 'introduceNum', [])
            parse(parser, varargin{:})
            particleNum = parser.Results.particleNum;
            continuous = parser.Results.continuous;
            endWith = parser.Results.endWith;
            fixedParticleNum = parser.Results.fixedParticleNum;
            introduceNum = parser.Results.introduceNum;
            
            % Number of bases
            baseNum = length(d);
            
            % Number of bases created
            baseCreatedNum = numel(fieldnames(obj.base));
            
            % Number of simulations performed
            simulatedNum = numel(fieldnames(obj.statistics));
            
            for i = simulatedNum+1:baseNum
                % Assert
                if baseCreatedNum-simulatedNum==1
                    assert(obj.base.(sprintf('base%i',i)).d50 == d(i));
                end
                
                % Change number of particles to introduce
                if fixedParticleNum
                    obj.simulatorParameters.storeNum = introduceNum(i);
                end
                
                % Create base
                obj = obj.createBase(d(i),'cuFine',1,'run',i);
                
                % Simulate filtration
                if or(~continuous, i==1)
                    x = 0;
                else
                    x = 1;
                end
                obj = obj.simulateFiltration('run',i,'followingRun',i-x);
                
                % Save every time a simulation is done
                switch particleNum
                    case 'normal'
                        if continuous
                            obj.saveResults(FileLocation.simulatorMonodisperse(obj.simulatorParameters.transportLogicType, ...
                                obj.filter.networkParameters.networkType, obj.filter.networkParameters.assignmentMethod, ...
                                obj.filter.filterParameters.cuCoarse, obj.filter.filterParameters.density, ...
                                obj.filter.filterParameters.size, obj.filter.filterParameters.triangulationMethod, ...
                                'cont'));
                        else
                            obj.saveResults(FileLocation.simulatorMonodisperse(obj.simulatorParameters.transportLogicType, ...
                                obj.filter.networkParameters.networkType, obj.filter.networkParameters.assignmentMethod, ...
                                obj.filter.filterParameters.cuCoarse, obj.filter.filterParameters.density, ...
                                obj.filter.filterParameters.size, obj.filter.filterParameters.triangulationMethod, ...
                                endWith));
                        end
                    case 'less'
                        if continuous
                            obj.saveResults(FileLocation.simulatorMonodisperseLess(obj.simulatorParameters.transportLogicType, ...
                                obj.filter.networkParameters.networkType, obj.filter.networkParameters.assignmentMethod, ...
                                obj.filter.filterParameters.cuCoarse, obj.filter.filterParameters.density, ...
                                obj.filter.filterParameters.size, obj.filter.filterParameters.triangulationMethod, ...
                                'cont'));
                        else
                            obj.saveResults(FileLocation.simulatorMonodisperseLess(obj.simulatorParameters.transportLogicType, ...
                                obj.filter.networkParameters.networkType, obj.filter.networkParameters.assignmentMethod, ...
                                obj.filter.filterParameters.cuCoarse, obj.filter.filterParameters.density, ...
                                obj.filter.filterParameters.size, obj.filter.filterParameters.triangulationMethod, ...
                                endWith));
                        end
                end
            end
        end
        function [obj] = runnerMonodisperse(obj, cuCoarse, d, varargin)
            % Input parser
            parser = inputParser;
            addParameter(parser, 'density', 'dense')
            addParameter(parser, 'size', [1,1,1])
            addParameter(parser, 'triangulationMethod', 'Del')
            addParameter(parser, 'networkType', 'actual')
            addParameter(parser, 'assignmentMethod', 'actual')
            expectedLoadMethod = Constants.loadMethod;
            addParameter(parser, 'loadMethod', 'read', @(x) any(validatestring(x, expectedLoadMethod)))
            expectedParticleNum = {'normal','less'};
            addParameter(parser, 'particleNum', 'normal', @(x) any(validatestring(x, expectedParticleNum)))
            addParameter(parser, 'continuous', 0)
            addParameter(parser, 'fixedParticleNum', 0)
            addParameter(parser, 'introduceNum', [])
            parse(parser, varargin{:})
            density = parser.Results.density;
            size = parser.Results.size;
            triangulationMethod = parser.Results.triangulationMethod;
            networkType = parser.Results.networkType;
            if strcmp(networkType, 'actual')
                assignmentMethod = parser.Results.assignmentMethod;
            else
                assignmentMethod = 'random';
            end
            loadMethod = parser.Results.loadMethod;
            particleNum = parser.Results.particleNum;
            continuous = parser.Results.continuous;
            fixedParticleNum = parser.Results.fixedParticleNum;
            introduceNum = parser.Results.introduceNum;
            if fixedParticleNum
                assert(length(introduceNum)==length(d))
                obj.simulatorParameters.introductionMethod = 'fixedParticleNum';
            end
            
            % File location
            transportLogicType = obj.simulatorParameters.transportLogicType;
            if continuous
                endWith = 'cont';
            elseif fixedParticleNum
                endWith = 'fixed';
            else
                endWith = '';
            end
            switch particleNum
                case 'normal'
                    s = FileLocation.simulatorMonodisperse(transportLogicType, ...
                        networkType, assignmentMethod, cuCoarse, density, size, triangulationMethod, endWith);
                case 'less'
                    s = FileLocation.simulatorMonodisperseLess(transportLogicType, ...
                        networkType, assignmentMethod, cuCoarse, density, size, triangulationMethod, endWith);
            end
            
            % Load or create filter
            if and(strcmp(loadMethod, 'read'), isfile(s))
                % Load data if available
                obj = Simulator.loadResults(s);
            else
                % Create filter
                obj = obj.createFilter(cuCoarse, 'density', density, 'size', size, ...
                    'triangulationMethod', triangulationMethod, 'networkType', networkType, ...
                    'assignmentMethod', assignmentMethod);
            end
            
            % Perform simulations and save
            obj = obj.simulateMonodisperse(d, 'particleNum', particleNum, 'continuous', continuous, ...
                'fixedParticleNum', fixedParticleNum, 'introduceNum', introduceNum, 'endWith', endWith);
        end
        function [obj] = repeatMonodisperse(obj, cuCoarse, d, varargin)
            % Input parser
            parser = inputParser;
            addParameter(parser, 'density', 'dense')
            addParameter(parser, 'size', [1,1,1])
            addParameter(parser, 'triangulationMethod', 'Del')
            addParameter(parser, 'networkType', 'actual')
            addParameter(parser, 'assignmentMethod', 'actual')
            expectedLoadMethod = Constants.loadMethod;
            addParameter(parser, 'loadMethod', 'read', @(x) any(validatestring(x, expectedLoadMethod)))
            expectedParticleNum = {'normal','less'};
            addParameter(parser, 'particleNum', 'normal', @(x) any(validatestring(x, expectedParticleNum)))
            addParameter(parser, 'continuous', 0)
            addParameter(parser, 'runNum', 10)
            parse(parser, varargin{:})
            density = parser.Results.density;
            size = parser.Results.size;
            triangulationMethod = parser.Results.triangulationMethod;
            networkType = parser.Results.networkType;
            if strcmp(networkType, 'actual')
                assignmentMethod = parser.Results.assignmentMethod;
            else
                assignmentMethod = 'random';
            end
            loadMethod = parser.Results.loadMethod;
            particleNum = parser.Results.particleNum;
            continuous = parser.Results.continuous;
            runNum = parser.Results.runNum;
            
            % Perform simulation for runNum times
            for i = 1:runNum
                % File names end with
                endWith = sprintf('%i',i);
                
                % File location
                switch particleNum
                    case 'normal'
                        if continuous
                            s = FileLocation.simulatorMonodisperse('flowbiased', ...
                                networkType, assignmentMethod, cuCoarse, density, size, triangulationMethod, 'cont');
                        else
                            s = FileLocation.simulatorMonodisperse('flowbiased', ...
                                networkType, assignmentMethod, cuCoarse, density, size, triangulationMethod, endWith);
                        end
                    case 'less'
                        if continuous
                            s = FileLocation.simulatorMonodisperseLess('flowbiased', ...
                                networkType, assignmentMethod, cuCoarse, density, size, triangulationMethod, 'cont');
                        else
                            s = FileLocation.simulatorMonodisperseLess('flowbiased', ...
                                networkType, assignmentMethod, cuCoarse, density, size, triangulationMethod, endWith);
                        end
                end
                
                % Load or create filter
                if and(strcmp(loadMethod, 'read'), isfile(s))
                    % Load data if available
                    obj = Simulator.loadResults(s);
                else
                    % Clear 
                    obj.base = struct;
                    obj.statistics = struct;
                    
                    % Create filter
                    obj = obj.createFilter(cuCoarse, 'density', density, 'size', size, ...
                        'triangulationMethod', triangulationMethod, 'networkType', networkType, ...
                        'assignmentMethod', assignmentMethod);
                end
                
                % Perform simulations and save
                obj = obj.simulateMonodisperse(d, 'particleNum', particleNum, 'continuous', continuous, 'endWith', endWith);
            end
        end
        
        % Packaging for monodisperse simulations in parts
        function [obj] = simulateMonodisperseParts(obj, d, varargin)
            % Input parser
            parser = inputParser;
            addParameter(parser, 'partNum', 10)
            addParameter(parser, 'divide', 0)
            parse(parser, varargin{:})
            partNum = parser.Results.partNum;
            divide = parser.Results.divide;
            
            if divide
                % Change the maxFinesTimeStep
                obj.simulatorParameters.maxFinesTimeStep = obj.simulatorParameters.maxFinesTimeStep/partNum;
            end
            
            for i = 1:partNum
                % Create base
                obj = obj.createBase(d,'cuFine',1,'run',i);
                
                % Simulate filtration
                if i==1
                    x = 0;
                else
                    x = 1;
                end
                obj = obj.simulateFiltration('run',i,'followingRun',i-x);
                
                % Save
                obj.saveResults(FileLocation.simulatorMonodisperseParts(obj.simulatorParameters.transportLogicType, ...
                    obj.filter.networkParameters.networkType, obj.filter.networkParameters.assignmentMethod, ...
                    obj.filter.filterParameters.cuCoarse, obj.filter.filterParameters.density, ...
                    obj.filter.filterParameters.size, obj.filter.filterParameters.triangulationMethod, d, divide));
            end
        end
        
        % Packaging for monodisperse simulations in parts and store paths
        function [obj] = simulateMonodisperseStorePaths(obj, d, storeNum)
            % Number of simulations
            simulationNum = length(storeNum);
            
            for i = 1:simulationNum
                % Update the number of particles to introduce
                obj.simulatorParameters.storeNum = storeNum(i);
                
                % Create base
                obj = obj.createBase(d,'cuFine',1,'run',i);
                
                % Perform simulation
                if i==1
                    x = 0;
                else
                    x = 1;
                end
                obj = obj.simulateFiltration('run',i,'followingRun',i-x);
                
                % Save
                obj.saveResults(FileLocation.simulatorMonodisperseStorePaths(obj.simulatorParameters.transportLogicType, ...
                    obj.filter.networkParameters.networkType, obj.filter.networkParameters.assignmentMethod, ...
                    obj.filter.filterParameters.cuCoarse, obj.filter.filterParameters.density, ...
                    obj.filter.filterParameters.size, obj.filter.filterParameters.triangulationMethod, d));
            end
        end
        
        % Packaging for polydisperse simulations
        function [obj] = simulatePolydisperse(obj, r, cuFine)
            % Number of simulations
            simulationNum = length(r);
            
            % Number of simulations performed
            simulatedNum = numel(fieldnames(obj.statistics));
            
            % Determine D15 for filter
            D15 = obj.filter.percentileParticleDiameter(15);
            
            for i = simulatedNum+1:simulationNum
                % Determine corresponding d50 for a specific ratio
                d50 = Base.givenRatio(r(i),D15,cuFine);
                
                % Create base
                obj = obj.createBase(d50,'cuFine',cuFine,'run',i);
                
                % Simulate filtration
                obj = obj.simulateFiltration('run',i,'followingRun',i);
                
                % Save every time a simulation is done
                obj.saveResults(FileLocation.simulatorPolydisperse(obj.simulatorParameters.transportLogicType, ...
                    obj.filter.networkParameters.networkType, obj.filter.networkParameters.assignmentMethod, ...
                    obj.filter.filterParameters.cuCoarse, obj.filter.filterParameters.density, ...
                    obj.filter.filterParameters.size, obj.filter.filterParameters.triangulationMethod, cuFine, r(i)));
            end
        end
        function [obj] = runnerPolydisperse(obj, cuCoarse, r, cuFine, varargin)
            % Input parser
            parser = inputParser;
            addParameter(parser, 'density', 'dense')
            addParameter(parser, 'size', [1,1,1])
            addParameter(parser, 'triangulationMethod', 'Del')
            addParameter(parser, 'networkType', 'actual')
            addParameter(parser, 'assignmentMethod', 'actual')
            expectedLoadMethod = Constants.loadMethod;
            addParameter(parser, 'loadMethod', 'read', @(x) any(validatestring(x, expectedLoadMethod)))
            parse(parser, varargin{:})
            density = parser.Results.density;
            size = parser.Results.size;
            triangulationMethod = parser.Results.triangulationMethod;
            networkType = parser.Results.networkType;
            if strcmp(networkType, 'actual')
                assignmentMethod = parser.Results.assignmentMethod;
            else
                assignmentMethod = 'random';
            end
            loadMethod = parser.Results.loadMethod;
            
            % File location
            transportLogicType = obj.simulatorParameters.transportLogicType;
            s = FileLocation.simulatorPolydisperse(transportLogicType, ...
                networkType, assignmentMethod, cuCoarse, density, size, triangulationMethod, cuFine);
            
            % Load or create filter
            if and(strcmp(loadMethod, 'read'), isfile(s))
                % Load data if available
                obj = Simulator.loadResults(s);
            else
                % Create filter
                obj = obj.createFilter(cuCoarse, 'density', density, 'size', size, ...
                    'triangulationMethod', triangulationMethod, 'networkType', networkType, ...
                    'assignmentMethod', assignmentMethod);
            end
            
            % Perform simulations and save
            obj = obj.simulatePolydisperse(r, cuFine);
        end
        function [obj] = repeatPolydisperse(obj, cuCoarse, r, cuFine, varargin)
            % Input parser
            parser = inputParser;
            addParameter(parser, 'density', 'dense')
            addParameter(parser, 'size', [1,1,1])
            addParameter(parser, 'triangulationMethod', 'Del')
            addParameter(parser, 'networkType', 'actual')
            addParameter(parser, 'assignmentMethod', 'actual')
            expectedLoadMethod = Constants.loadMethod;
            addParameter(parser, 'loadMethod', 'read', @(x) any(validatestring(x, expectedLoadMethod)))
            addParameter(parser, 'partNum', 10)
            addParameter(parser, 'run', 1)
            parse(parser, varargin{:})
            density = parser.Results.density;
            size = parser.Results.size;
            triangulationMethod = parser.Results.triangulationMethod;
            networkType = parser.Results.networkType;
            if strcmp(networkType, 'actual')
                assignmentMethod = parser.Results.assignmentMethod;
            else
                assignmentMethod = 'random';
            end
            loadMethod = parser.Results.loadMethod;
            partNum = parser.Results.partNum;
            run = parser.Results.run;
            
            % File location
            transportLogicType = obj.simulatorParameters.transportLogicType;
            
            for i = 1:length(r)
                % File location
                s = FileLocation.simulatorPolydisperseRepeat(transportLogicType, ...
                    networkType, assignmentMethod, cuCoarse, density, size, triangulationMethod, cuFine, r(i), sprintf('continuous%i',run));
                
                % Load or create filter
                if and(strcmp(loadMethod, 'read'), isfile(s))
                    % Load data if available
                    obj = Simulator.loadResults(s);
                else
                    % Create filter
                    obj = obj.createFilter(cuCoarse, 'density', density, 'size', size, ...
                        'triangulationMethod', triangulationMethod, 'networkType', networkType, ...
                        'assignmentMethod', assignmentMethod);
                    
                    % Clear data
                    obj.base = struct;
                    obj.statistics = struct;
                end
                
                % Determine D15 for filter
                D15 = obj.filter.percentileParticleDiameterVol(15);
                
                % Number of simulations performed
                simulatedNum = numel(fieldnames(obj.statistics));
                
                % Determine corresponding d50 for a specific ratio
                d50 = Base.givenRatio(r(i),D15,cuFine);
                
                for j = simulatedNum+1:partNum
                    
                    % Create base
                    obj = obj.createBase(d50,'cuFine',cuFine,'run',j);
                    
                    if j==1
                        following = j;
                    else
                        following = j-1;
                    end
                    
                    % Simulate filtration
                    obj = obj.simulateFiltration('run',j,'followingRun',following);
                    
                    % Save every time a simulation is done
                    obj.saveResults(FileLocation.simulatorPolydisperseRepeat(obj.simulatorParameters.transportLogicType, ...
                        obj.filter.networkParameters.networkType, obj.filter.networkParameters.assignmentMethod, ...
                        obj.filter.filterParameters.cuCoarse, obj.filter.filterParameters.density, ...
                        obj.filter.filterParameters.size, obj.filter.filterParameters.triangulationMethod, cuFine, r(i), sprintf('continuous%i',run)));
                end
            end
        end
        
        % ===== Statistics ================================================
        % Number of particles introduced
        function [n] = introducedNum(obj)
            % Number of simulations
            simulationNum = numel(fieldnames(obj.statistics));
            
            % Initialise array(s)
            n = zeros(simulationNum,1);
            
            % Determine number
            for i = 1:simulationNum
                n(i) = obj.statistics.(sprintf('statistics%i',i)).timeIntroduced(end,1);
            end
        end
        
        % Number of particles retained in each layer
        function [n] = layerRetainedNum(obj, layerNum)
            % Number of simulations
            simulationNum = numel(fieldnames(obj.statistics));
            
            % Initialise array(s)
            n = zeros(layerNum,simulationNum);
            
            % For each simulation
            for i = 1:simulationNum
                % Determine the number of particles retained in each layer
                n(:,i) = obj.statistics.(sprintf('statistics%i',i)).layerRetainedNum(obj.filter, ...
                    layerNum, obj.filter.networkParameters.flowDirection);
            end
        end
        
        % Percentage particle eroded, either by number or volume. Overall
        % is used when consecutive simulations, such that the percentage
        % is determined from all Statistics object under the property
        % 'statistics', i.e. statistics1, statistics2 ...
        function [p] = percentageErodedNum(obj, varargin)
            % Input parser
            parser = inputParser;
            addParameter(parser, 'continuous', 0)
            parse(parser, varargin{:})
            continuous = parser.Results.continuous;
            
            % Number of simulations
            simulationNum = numel(fieldnames(obj.statistics));
            
            % Initialise array(s)
            p = zeros(simulationNum,1);
            
            % Determine percentages
            for i = 1:simulationNum
                if or(~continuous, i==1)
                    p(i) = obj.statistics.(sprintf('statistics%i',i)).percentageErodedNum(obj.filter.nodeInlet);
                else
                    p(i) = obj.statistics.(sprintf('statistics%i',i)).percentageErodedNum(obj.filter.nodeInlet, ...
                        'previousStatistics',obj.statistics.(sprintf('statistics%i',i-1)));
                end
            end
        end
        function [p] = percentageErodedNumOverall(obj)
            % Number of simulations
            simulationNum = numel(fieldnames(obj.statistics));
            
            % Initialise array(s)
            nEroded = zeros(simulationNum,1);
            nIntroduced = nEroded;
            
            % Determine number of particles eroded and introduced
            for i = 1:simulationNum
                [~, nEroded(i), nIntroduced(i)] = obj.statistics.(sprintf('statistics%i',i)).percentageErodedNum(obj.filter.nodeInlet, 'overall', 1);
            end
            
            % Determine percentage
            p = sum(nEroded)/sum(nIntroduced)*100;
        end
        function [p] = percentageErodedVol(obj)
            % Number of simulations
            simulationNum = numel(fieldnames(obj.statistics));
            
            % Initialise array(s)
            p = zeros(simulationNum,1);
            
            % Determine percentages
            for i = 1:simulationNum
                p(i) = obj.statistics.(sprintf('statistics%i',i)).percentageErodedVol;
            end
        end
        
        % Total volume of base particles eroded. Area is for the volume
        % that is normalised by the cross-sectional area of the filter
        % sample
        function [v] = totalErodedVol(obj)
            % Number of simulations
            simulationNum = numel(fieldnames(obj.statistics));
            
            % Total volume of void space
            totalV = sum(obj.filter.nodeVolume)*obj.simulatorParameters.densityLimit;
            
            % Initialise array(s)
            v = zeros(simulationNum,1);
            
            % Determine volumes
            for i = 1:simulationNum
                v(i) = obj.statistics.(sprintf('statistics%i',i)).timeEroded(end,2);
            end
            
            % Total volume
            v = sum(v)/totalV;
        end
        function [v] = totalErodedVolArea(obj)
            % Number of simulations
            simulationNum = numel(fieldnames(obj.statistics));
            
            % Area of inlet
            area = (obj.filter.boundary(4)-obj.filter.boundary(3))*(obj.filter.boundary(2)-obj.filter.boundary(1));
            
            % Initialise array(s)
            v = zeros(simulationNum,1);
            
            % Determine volumes
            for i = 1:simulationNum
                v(i) = obj.statistics.(sprintf('statistics%i',i)).timeEroded(end,2);
            end
            
            % Total volume
            v = sum(v)/area;
        end
        
        % Percentage of edge that are clogged
        function [p] = percentageEdgeClogged(obj, mechanism)
            % Number of simulations
            simulationNum = numel(fieldnames(obj.statistics));
            
            % Initialise array(s)
            p = zeros(simulationNum,1);
            
            % Determine percentages
            for i = 1:simulationNum
                p(i) = obj.statistics.(sprintf('statistics%i',i)).percentageEdgeClogged(mechanism);
            end
        end
        
        % Distance travelled by particles (options available for eroded or
        % retained particles)
        % directness = distance travelled normalised by shortest path 
        %              (a measure of path tortuosity)
        % distanceTravelled = not normalised
        function [r] = directness(obj, varargin)
            % Input parser
            parser = inputParser;
            expectedMechanism = {'eroded','retained','retained1','retained2','retained3'};
            addParameter(parser, 'mechanism', 'eroded',  @(x) any(validatestring(x, expectedMechanism)))
            addParameter(parser, 'd', [])
            parse(parser, varargin{:})
            mechanism = parser.Results.mechanism;
            d = parser.Results.d;
            
            % Determine the shortest path from inlet to outlet
            if isempty(d)
                [~,d] = obj.filter.shortestPath();
            end
            
            % Number of simulations
            simulationNum = numel(fieldnames(obj.statistics));
            
            % Initialise array(s)
            r = zeros(simulationNum,1);
            
            % Determine average distance
            for i = 1:simulationNum
                r(i) = obj.statistics.(sprintf('statistics%i',i)).averageDistance('mechanism',mechanism);
            end
            
            % Determine ratio
            r = r/d;
        end
        function [d] = distanceTravelled(obj, varargin)
            % Input parser
            parser = inputParser;
            expectedMechanism = {'eroded','retained','retained1','retained2','retained3'};
            addParameter(parser, 'mechanism', 'eroded',  @(x) any(validatestring(x, expectedMechanism)));
            parse(parser, varargin{:});
            mechanism = parser.Results.mechanism;
            
            % Number of simulations
            simulationNum = numel(fieldnames(obj.statistics));
            
            % Initialise array(s)
            d = zeros(simulationNum,1);
            
            % Determine average distance
            for i = 1:simulationNum
                d(i) = obj.statistics.(sprintf('statistics%i',i)).averageDistance('mechanism',mechanism);
            end
        end
        
        % Average number of constrictions met by particles (options 
        % available for eroded or retained particles)
        function [n] = constrictionNum(obj, varargin)
            % Input parser
            parser = inputParser;
            expectedMechanism = {'eroded','retained','retained1','retained2','retained3'};
            addParameter(parser, 'mechanism', 'eroded',  @(x) any(validatestring(x, expectedMechanism)));
            parse(parser, varargin{:});
            mechanism = parser.Results.mechanism;
            
            % Number of simulations
            simulationNum = numel(fieldnames(obj.statistics));
            
            % Initialise array(s)
            n = zeros(simulationNum,1);
            
            % Determine average distance
            for i = 1:simulationNum
                n(i) = obj.statistics.(sprintf('statistics%i',i)).constrictionNum('mechanism',mechanism);
            end
        end
        
        % Extract the diameters of base used in simulations
        function [d] = monodisperseBaseDiameter(obj)
            % Number of bases
            baseNum = numel(fieldnames(obj.base));
            
            % Initialise array(s)
            d = zeros(baseNum,1);
            
            % Determine percentages
            for i = 1:baseNum
                d(i) = obj.base.(sprintf('base%i',i)).d50;
            end
        end
        
        % Firstly, determine the number of exits that are smaller than a 
        %          specified d for all nodes
        % Secondly, determine the percentage of nodes that are affected by 
        %           premature node clogging, i.e. nodes becoming full
        %           before all exits are clogged
        function [m] = numberExitsSmaller(obj, d)
            % Unpack variables
            [flowDirection, networkType] = obj.unpack(obj.filter.networkParameters, 'flowDirection', 'networkType');
            [nodeNum, nodeCoord, nodeInlet, nodeLattice, nodeOutlet, nodeSurfaceArea, nodeVolume] = obj.unpack(obj.filter, 'nodeNum', 'nodeCoord', 'nodeInlet', 'nodeLattice', 'nodeOutlet', 'nodeSurfaceArea', 'nodeVolume');
            [nodeConnectivity, nodeEdgeTo, nodeNodeTo, ~, nodeEdgeDirection] = obj.unpack(obj.filter, 'nodeConnectivity', 'nodeEdgeTo', 'nodeNodeTo', 'nodeEdgeToPos', 'nodeEdgeDirection');
            [edgeNum, edgeArea, edgeDiameter, edgeLength, edgeList, edgeReduced] = obj.unpack(obj.filter, 'edgeNum', 'edgeArea', 'edgeDiameter', 'edgeLength', 'edgeList', 'edgeReduced');
            [pressureEntry, pressureExit, flowMethod] = obj.unpack(obj.flowSolverParameters, 'pressureEntry', 'pressureExit', 'flowMethod');
            
            % Initialise array(s)
            statisticsTemp = Statistics(nodeNum, edgeNum, nodeLattice, 1, 1);
            edgeFlowing = statisticsTemp.edgeFlowing;
            nodeClogging = statisticsTemp.nodeClogging;
            nodeRetained = statisticsTemp.nodeRetained;
            [nodeActive, ~, ~, ~] = initialiseArrays(nodeNum, nodeClogging);
            
            % Determine flow
            [nodePressure, ~, edgeConductance] = obj.determinePressure(nodeNum, edgeNum, nodeActive, nodeConnectivity, nodeCoord, nodeEdgeTo, nodeInlet, nodeNodeTo, nodeOutlet, nodeRetained, nodeSurfaceArea, nodeVolume, edgeArea, edgeFlowing, edgeLength, edgeList, edgeReduced, pressureEntry, pressureExit, flowMethod, flowDirection, networkType);
            [edgeFlow] = obj.determineEdgeFlow(edgeNum, nodePressure, edgeConductance, edgeList);
            
            % Initialise m
            m = zeros(nodeNum,length(d));
            
            for i = 1:nodeNum
                % Determine edges connected to the node and outflow
                n = nodeConnectivity(i);
                edgeTo = nodeEdgeTo(i,1:n)';
                outflow = edgeFlow(edgeTo).*nodeEdgeDirection(i,1:n)';
                idxEntrances = outflow<=0;
                
                for j = 1:length(d)
                    % Exits smaller than base particles
                    edgeToDiameter = edgeDiameter(edgeTo)<d(j);
                    
                    % Determine number of open exits
                    m(i,j) = sum(and(~idxEntrances, edgeToDiameter));
                end
            end
        end
        function [p] = percentageNodeAffected(obj)
            % Adopted tranport logic type
            transportLogicType = obj.simulatorParameters.transportLogicType;
            
            % Determine the base particle diameter
            d = obj.monodisperseBaseDiameter();
            
            % Determine the mean number of exits
            m = obj.numberExitsSmaller(d);
            
            % Indices of inlets or outlet
            idxIgnore = or(obj.filter.nodeInlet, obj.filter.nodeOutlet);
            
            if or(strcmp(transportLogicType, 'flowbiasednode'),strcmp(transportLogicType, 'flowbiasedextreme'))
                % Determine the percentage of nodes affected, ignore inlets and
                % outlets
                p = zeros(length(d),1);
                for i = 1:length(d)
                    p(i) = sum(and(m(:,i)>0, ~idxIgnore))/(obj.filter.nodeNum-sum(idxIgnore))*100;
                end
            else
                % Density limit
                f = obj.simulatorParameters.densityLimit;
                
                % Determine the percentage of nodes affected, ignore inlets and
                % outlets
                p = zeros(length(d),1);
                for i = 1:length(d)
                    p(i) = sum(and(obj.filter.nodeVolume*f<(m(:,i)+1)*(pi/6)*(d(i)^3), ~idxIgnore))/(obj.filter.nodeNum-sum(idxIgnore))*100;
                end
            end
        end
        
        % Mean size of clusters of fully-filled nodes
        function [m] = meanClusterSize(obj)
            % Determine the base particle diameter
            d = obj.monodisperseBaseDiameter();
            
            % Initialise m
            m = zeros(length(d),1);
            
            for i = 1:length(d)
                % Determine the cluster sizes
                [~, binsize] = obj.statistics.(sprintf('statistics%i',i)).monoNodeFullClusterSize(d(i), ...
                    obj.filter, obj.simulatorParameters.densityLimit);
                m(i) = mean(binsize);
            end
        end
        
        % Mean size of clusters of nodes with the largest edge diameter
        % smaller than the simulated monodisperse base particle diameter d
        function [m] = meanClusterEdgeDiameter(obj)
            % Determine the base particle diameter
            d = obj.monodisperseBaseDiameter();
            
            % Initialise m
            m = zeros(length(d),1);
            
            for i = 1:length(d)
                % Obtain the size of clusters of nodes with largest edge
                % smaller than d
                [s] = obj.filter.nodeSmaller(d(i));
                m(i) = mean(s);
            end
        end
        
        % Percentage of full nodes or edge clogged in each layer
        function [p] = percentageLayerNodeFull(obj, layerNum)
            % Number of simulations
            simulationNum = numel(fieldnames(obj.statistics));
            
            % Initialise p
            p = zeros(layerNum, simulationNum);
            
            % For each simulation
            for i = 1:simulationNum
                % Obtain the percentage clogged in each layer
                p(:,i) = obj.statistics.(sprintf('statistics%i',i)).percentageLayerFull(obj.base.(sprintf('base%i',i)).d50, obj.filter, ...
                    obj.simulatorParameters.densityLimit, layerNum);
            end
        end
        function [p] = percentageLayerEdgeClogged(obj, layerNum, mechanism)
            % Number of simulations
            simulationNum = numel(fieldnames(obj.statistics));
            
            % Initialise p
            p = zeros(layerNum, simulationNum);
            
            % For each simulation
            for i = 1:simulationNum
                % Obtain the percentage clogged in each layer
                p(:,i) = obj.statistics.(sprintf('statistics%i',i)).percentageLayerClogged(obj.filter, ...
                    layerNum, mechanism);
            end
        end
        
        % Plot the variation of number of particles introduced with time
        % for all simulations
        function [] = plotTimeNum(obj, varargin)
            % Input parser
            parser = inputParser;
            addParameter(parser, 'continuous', 0)
            parse(parser, varargin{:})
            continuous = parser.Results.continuous;
            
            % Number of simulations
            simulationNum = numel(fieldnames(obj.statistics));
            
            % For each simulation
            for i = 1:simulationNum
                if or(~continuous, i==1)
                    obj.statistics.(sprintf('statistics%i',i)).plotTimeNum(obj.filter.nodeInlet);
                else
                    obj.statistics.(sprintf('statistics%i',i)).plotTimeNum(obj.filter.nodeInlet, ...
                        'previousStatistics',obj.statistics.(sprintf('statistics%i',i-1)));
                end
                
            end
        end
        
        % Determine flow, either with a 'clean' filter, i.e. without any
        % retained particles, or following a simulation (input the
        % simulation number)
        function [edgeFlow, nodePressure] = flowFilter(obj, varargin)
            % Input parser
            parser = inputParser;
            addParameter(parser, 'followingRun', 0)
            parse(parser, varargin{:})
            followingRun = parser.Results.followingRun;
            
            % Unpack variables
            [flowDirection, networkType] = obj.unpack(obj.filter.networkParameters, 'flowDirection', 'networkType');
            [nodeNum, nodeCoord, nodeInlet, nodeLattice, nodeOutlet, nodeSurfaceArea, nodeVolume] = obj.unpack(obj.filter, 'nodeNum', 'nodeCoord', 'nodeInlet', 'nodeLattice', 'nodeOutlet', 'nodeSurfaceArea', 'nodeVolume');
            [nodeConnectivity, nodeEdgeTo, nodeNodeTo, ~, ~] = obj.unpack(obj.filter, 'nodeConnectivity', 'nodeEdgeTo', 'nodeNodeTo', 'nodeEdgeToPos', 'nodeEdgeDirection');
            [edgeNum, edgeArea, ~, edgeLength, edgeList, edgeReduced] = obj.unpack(obj.filter, 'edgeNum', 'edgeArea', 'edgeDiameter', 'edgeLength', 'edgeList', 'edgeReduced');
            [pressureEntry, pressureExit, flowMethod] = obj.unpack(obj.flowSolverParameters, 'pressureEntry', 'pressureExit', 'flowMethod');
            
            % Initialise array(s)
            [nodeActive, F, ~, ~] = initialiseArrays(nodeNum, 0);
            if followingRun==0
                statisticsTemp = Statistics(nodeNum, edgeNum, nodeLattice, 1, 1);
                edgeFlowing = statisticsTemp.edgeFlowing;
                
                statisticsTemp = Statistics(nodeNum, edgeNum, nodeLattice, 0, 0);
                nodeRetained = statisticsTemp.nodeRetained;
            else
                [nodeRetained, edgeFlowing] = obj.unpack(obj.statistics.(sprintf('statistics%i',followingRun)), 'nodeRetained', 'edgeFlowing');
                
                % Update nodeActive
                [nodeActive] = updateNodeActive(nodeNum, nodeActive, nodeRetained, F);
            end
            
            % Determine flow
            [nodePressure, ~, edgeConductance] = obj.determinePressure(nodeNum, edgeNum, nodeActive, nodeConnectivity, nodeCoord, nodeEdgeTo, nodeInlet, nodeNodeTo, nodeOutlet, nodeRetained, nodeSurfaceArea, nodeVolume, edgeArea, edgeFlowing, edgeLength, edgeList, edgeReduced, pressureEntry, pressureExit, flowMethod, flowDirection, networkType);
            [edgeFlow] = obj.determineEdgeFlow(edgeNum, nodePressure, edgeConductance, edgeList);
        end
        
        % Percolation threshold using a directed graph
        function [p, d] = percolationThresholdAssignedDirected(obj, varargin)
            % Input parser
            parser = inputParser;
            addParameter(parser, 'minP', 0)
            addParameter(parser, 'display', 'off')
            addParameter(parser, 'pos', 0)
            parse(parser, varargin{:});
            minP = parser.Results.minP;
            display = parser.Results.display;
            pos = parser.Results.pos;
            
            % Calculate flow
            [edgeFlow,~] = obj.flowFilter;
            
            % Determine the direction of flow in each edge
            idxPositive = edgeFlow>0;
            idxNegative = edgeFlow<0;
            
            % Restrict to positive z-direction
            if pos
                % Initialise edgeDirection
                edgeDirection = -1*ones(obj.filter.edgeNum,1);
                
                % Determine the direction for each edge
                for i = 1:obj.filter.edgeNum
                    if nodeCoord(edgeList(i,2),3)-nodeCoord(edgeList(i,1),3)>0
                        edgeDirection(i) = 1;
                    end
                end
                
                % Update idxPositive and idxNegative
                idxPositive = and(idxPositive,edgeDirection==1);
                idxNegative = and(idxNegative,edgeDirection==-1);
            end
            
            % Determine percolation threshold
            [p, d] = obj.filter.percolationThresholdAssignedDirected(idxPositive, idxNegative, ...
                'minP', minP, 'display', display);
        end
        
        % Similarity of edgePathEroded for all simulations. Two methods
        % were implemented. The method 'similarityIndexEdgePathEroded'
        % compares the edge indices as proposed in Sufian et. al. (2018)
        function [r] = similarityBoolEdgePathEroded(obj)
            % Number of simulations
            simulationNum = numel(fieldnames(obj.statistics));
            
            % List of combinations
            combList = combnk(1:simulationNum,2);
            combNum = size(combList,1);
            
            % Initialise array(s)
            r = zeros(combNum,1);
            
            % For each combination
            for i = 1:combNum
                % Determine the similarity
                r(i) = Statistics.similarityBool(obj.statistics.(sprintf('statistics%i',combList(i,1))).edgePathEroded>1, ...
                    obj.statistics.(sprintf('statistics%i',combList(i,2))).edgePathEroded>1);
            end
            
            % Average
            r = mean(r);
        end
        function [r] = similarityIndexEdgePathEroded(obj)
            % Number of simulations
            simulationNum = numel(fieldnames(obj.statistics));
            
            % List of combinations
            combList = combnk(1:simulationNum,2);
            combNum = size(combList,1);
            
            % Initialise array(s)
            r = zeros(combNum,1);
            
            % For each combination
            for i = 1:combNum
                % Determine the similarity
                r(i) = Statistics.similarityIndex(obj.statistics.(sprintf('statistics%i',combList(i,1))).edgePathEroded>1, ...
                    obj.statistics.(sprintf('statistics%i',combList(i,2))).edgePathEroded>1);
            end
            
            % Average
            r = mean(r);
        end
        
        % Determine ratios of particle sizes
        function [r] = terzaghiRatio(obj)
            % Number of simulations
            simulationNum = numel(fieldnames(obj.statistics));
            
            % Initialise r
            r = zeros(simulationNum,1);
            
            for i = 1:simulationNum
                % psdFineVol and psdCourseVol
                psdFineVol = obj.base.(sprintf('base%i',i)).psdFineVol;
                
                % Determine Terzaghi's ratio
                d85 = interp1(psdFineVol(:,1),psdFineVol(:,2),85);
                D15 = obj.filter.percentileParticleDiameterVol(15);
                r(i) = D15/d85;
            end
        end
        function [r] = fosterFellRatio(obj)
            % Number of simulations
            simulationNum = numel(fieldnames(obj.statistics));
            
            % Initialise r
            r = zeros(simulationNum,1);
            
            for i = 1:simulationNum
                % psdFineVol and psdCourseVol
                psdFineVol = obj.base.(sprintf('base%i',i)).psdFineVol;
                
                % Determine Terzaghi's ratio
                d95 = interp1(psdFineVol(:,1),psdFineVol(:,2),95);
                D15 = obj.filter.percentileParticleDiameterVol(15);
                r(i) = D15/d95;
            end
        end
        function [r] = testRatio1(obj)
            % Number of simulations
            simulationNum = numel(fieldnames(obj.statistics));
            
            % Initialise r
            r = zeros(simulationNum,1);
            
            for i = 1:simulationNum
                % psdFineVol and psdCourseVol
                psdFineArea = Statistics.convertPSDVolToArea(obj.base.(sprintf('base%i',i)).psdFineVol);
                
                % Determine Terzaghi's ratio
                d85 = interp1(psdFineArea(:,1),psdFineArea(:,2),85);
                D15 = obj.filter.percentileParticleDiameterVol(15);
                r(i) = D15/d85;
            end
        end
        
        % Determine the ratios (constriction based method)
        function [r] = testRatio2(obj)
            % Number of simulations
            simulationNum = numel(fieldnames(obj.statistics));
            
            % Initialise r
            r = zeros(simulationNum,1);
            
            for i = 1:simulationNum
                % psdFineVol and psdCourseVol
                psdFineVol = obj.base.(sprintf('base%i',i)).psdFineVol;
                csd = obj.filter.csdNum;
                
                % Determine Terzaghi's ratio
                d85 = interp1(psdFineVol(:,1),psdFineVol(:,2),85);
                D35c = interp1(csd(:,1),csd(:,2),35);
                r(i) = D35c/d85;
            end
        end
        function [r] = testRatio3(obj)
            % Number of simulations
            simulationNum = numel(fieldnames(obj.statistics));
            
            % Initialise r
            r = zeros(simulationNum,1);
            
            for i = 1:simulationNum
                % psdFineVol and psdCourseVol
                psdFineVol = obj.base.(sprintf('base%i',i)).psdFineVol;
                csd = obj.filter.csdNum;
                
                % Determine Terzaghi's ratio
                d95 = interp1(psdFineVol(:,1),psdFineVol(:,2),95);
                D35c = interp1(csd(:,1),csd(:,2),35);
                r(i) = D35c/d95;
            end
        end
        function [r] = indraratnaRautRatio(obj)
            % Number of simulations
            simulationNum = numel(fieldnames(obj.statistics));
            
            % Initialise r
            r = zeros(simulationNum,1);
            
            for i = 1:simulationNum
                % psdFineArea and csd
                psdFineArea = Statistics.convertPSDVolToArea(obj.base.(sprintf('base%i',i)).psdFineVol);
                % psdFineNum = obj.base.(sprintf('base%i',i)).psdFineNum;
                csd = obj.filter.csdNum;
                
                % Determine Terzaghi's ratio
                d85 = interp1(psdFineArea(:,1),psdFineArea(:,2),85);
                % d85 = interp1(psdFineNum(:,1),psdFineNum(:,2),90);
                D35c = interp1(csd(:,1),csd(:,2),35);
                r(i) = D35c/d85;
            end
        end
        
        % Whether the filter is seal at the end of simulations
        function [bool] = filterSeal(obj)
            % Number of simulations
            simulationNum = numel(fieldnames(obj.statistics));
            
            % Initialise bool
            bool = zeros(simulationNum,1);
            
            for i = 1:simulationNum
                bool(i) = obj.statistics.(sprintf('statistics%i',i)).inletFull(obj.filter,obj.simulatorParameters.densityLimit);
            end
        end
        function [category] = filterSealTime(obj)
            % Categories defined as: 
            % 1: Seal with little erosion i.e. seal after first simulation
            % run
            % 2: Seal with signigicant erosion i.e. not seal after first,
            % but seal after last simulation
            % 3: No seal i.e. not seal even after 10 consecutive
            % simulations
            
            % Number of simulations
            simulationNum = numel(fieldnames(obj.statistics));
            
            % Initialise bool
            bool = zeros(2,1);
            
            j = 0;
            for i = [1,simulationNum]
                j = j+1;
                bool(j) = obj.statistics.(sprintf('statistics%i',i)).inletFull(obj.filter,obj.simulatorParameters.densityLimit);
            end
            
            if all(bool)
                category = 1;
            elseif bool(2)==1
                category = 2;
            else
                category = 3;
            end
        end
        
        % Largest overall particle eroded (95 implies using the 95th
        % percentile instead, which was considered in Moraci et. al. (2011)
        function [d] = largestParticleEroded(obj)
            % Number of simulations
            simulationNum = numel(fieldnames(obj.statistics));
            
            % Initialise d
            d = zeros(simulationNum,1);
            
            for i = 1:simulationNum
                d(i) = obj.statistics.(sprintf('statistics%i',i)).largestParticle;
            end
            
            d = max(d);
        end
        function [d] = largestParticleEroded95(obj)
            % Number of simulations
            simulationNum = numel(fieldnames(obj.statistics));
            
            % Initialise d
            d = zeros(simulationNum,1);
            
            for i = 1:simulationNum
                d(i) = obj.statistics.(sprintf('statistics%i',i)).largestParticle95(obj.base.(sprintf('base%i',i)).psdFineBucket);
            end
            
            d = max(d);
        end
        
        % ===== Saving data ===============================================
        % Reduce the size and save
        function [] = saveResults(obj, saveLocation)
            % Reduce size of Filter object
            obj.filter = obj.filter.removeProperties();
            
            % Reduce size of Base object(s)
            baseNum = numel(fieldnames(obj.base));
            for i = 1:baseNum
                obj.base.(sprintf('base%i',i)) = obj.base.(sprintf('base%i',i)).removePSDs();
            end
            
            % Save
            simulator = obj;
            save(saveLocation, 'simulator', '-v7.3');
        end
    end
    
    methods (Static)
        % ===== Flow Solver ===============================================
        % Calculate the conductance of the edges and the pressure at the
        % nodes; main steps:
        % (1) Determine edge conductance based on Chareyre et. al. (2012)
        % (2) Based on conservation of mass, determine the pressures at the
        %     nodes, and the source into the nodes at the boundaries
        function [nodePressure, nodeSource, edgeConductance] = determinePressure(nodeNum, edgeNum, nodeActive, nodeConnectivity, nodeCoord, nodeEdgeTo, nodeInlet, nodeNodeTo, nodeOutlet, nodeRetained, nodeSurfaceArea, nodeVolume, edgeArea, edgeFlowing, edgeLength, edgeList, edgeReduced, pressureEntry, pressureExit, flowMethod, flowDirection, networkType) %#ok<INUSL>
            % The function solves for nodePressure, nodeSource and edgeConductance
            
            % The matrix nodeActive stores:
            % 1: Total number of particles active in the node
            % 2: Total surface area of particles active or retained in the node
            % 3: Total volume of particles active or retained in the node
            if 1
                nodeFinesArea = nodeActive(:,2);
                nodeFinesVolume = nodeActive(:,3);
            else
                % IMPROVEMENT: total flow more meaningful if active
                % particles are not considered
                nodeFinesArea = nodeRetained(:,2); %#ok<UNRCH>
                nodeFinesVolume = nodeRetained(:,3);
            end
            
            [edgeConductance] = setNetworkConductance(edgeFlowing(:,1), edgeList, nodeVolume, nodeSurfaceArea, edgeArea, edgeLength, nodeFinesVolume, nodeFinesArea, edgeReduced);
            [nodePressure, nodeSource] = calculatePressure(flowMethod, nodeNum, edgeNum, nodeConnectivity, nodeInlet, nodeOutlet, nodeEdgeTo, nodeNodeTo, edgeConductance, pressureEntry, pressureExit, edgeList, nodeCoord, flowDirection, networkType);
        end
        
        % Calculate the flow in each edge based on the conductance and
        % pressure calculated with Static method 'determinePressure'
        function [edgeFlow] = determineEdgeFlow(edgeNum, nodePressure, edgeConductance, edgeList)
            % The function determines the flow in each edge
            
            % Initialise edgeFlow
            edgeFlow = zeros(edgeNum,1);
            
            for i = 1:edgeNum
                % Determine flow
                edgeFlow(i) = edgeConductance(i)*(nodePressure(edgeList(i,1))-nodePressure(edgeList(i,2)));
            end
        end
        
        % Calculate the outflows from a specified node, i.e. positice for 
        % outflow negative for inflow
        function [outflow] = determineOutflow(nodeFrom, nodeTo, edgeTo, edgeConductance, nodePressure)
            % Determine the outflows from a node
            outflow = zeros(size(nodeTo));
            for i = 1:length(nodeTo)
                outflow(i) = edgeConductance(edgeTo(i))*(nodePressure(nodeFrom)-nodePressure(nodeTo(i)));
            end
        end
        
        % ===== Unpack elements ===========================================
        function [varargout] = unpack(obj, varargin)
            % The function unpacks the struct
            
            varargout = cell(1,length(varargin));
            for i = 1:length(varargin)
                varargout{i} = obj.(sprintf('%s', varargin{i}));
            end
        end
        
        % ===== Loading data ==============================================
        function [simulator] = loadResults(saveLocation)
            % Load
            simulator = load(saveLocation);
            simulator = simulator.simulator;
            
            % Reload other data
            simulator.filter = simulator.filter.regenerateProperties('read');
            baseNum = numel(fieldnames(simulator.base));
            for i = 1:baseNum
                simulator.base.(sprintf('base%i',i)) = simulator.base.(sprintf('base%i',i)).generatePSDs();
            end
            
            % Save saveLocation
            simulator.fileLocation = saveLocation;
        end
    end
end

% PRIVATE FUNCTIONS
% ===== Initialise Arrays =================================================
function [nodeActive, nodeClogging, F, idxParticle] = initialiseArrays(nodeNum, nodeClogging)
% The function creates empty arrays to store data

% The matrix nodeActive stores:
% 1: Total number of particles active in the node
% 2: Total surface area of particles active or retained in the node
% 3: Total volume of particles active or retained in the node
nodeActive = zeros(nodeNum, 3);

% The matrix F stores:
% 1: Node index
% 2: Particle diameter
% 3: Particle surface area
% 4: Particle volume
% 5: Particle status (-1:eroded, 0:active, 1:clogged, 2:retained, 3:nospace)
% 6: Number of constrictions met
% 7: Distance travelled
F = [];

% Store the index of particles
idxParticle = [];

% Set the time steps at which node is clogged to 0 if it is clogged in
% previous simulations
nodeClogging(~isnan(nodeClogging(:,1)),1) = 0;
end

% ===== Introducing Fines =================================================
function [nodeActive, F, timeIntroduced, timeOutside, particleRoute, idxParticle, maxFinesTimeStep] = introduceFinesParticle(edgeNum, nodeActive, nodeInlet, nodeVolume, F, timeIntroduced, timeOutside, particleRoute, idxParticle, psdFineNum, psdFineBucket, introductionMethod, densityLimit, contentLimit, maxFinesParticle, maxFinesTimeStep, storeOption, storeNum, particleInlet, t)
% This function generates fines from a PSD by number and at the entry face

% The matrix nodeActive stores:
% 1: Total number of particles active in the node
% 2: Total surface area of particles active in the node
% 3: Total volume of particles active in the node

% The matrix F stores:
% 1: Node index
% 2: Particle diameter
% 3: Particle surface area
% 4: Particle volume
% 5: Particle status (-1:eroded, 0:active, 1:clogged, 2:retained, 3:nospace)
% 6: Number of constrictions
% 7: Distance travelled

% PSD by number of base particles
PSD = psdFineNum(:,2);

% Introduce
switch introductionMethod
    case 'fixedParticleNum'
        if strcmp(storeOption, 'storePath')
            % storeNum of particles
            numFinesParticle = storeNum-length(idxParticle);
            
            % Place fines at the face where lattice z = 1
            if t==1
                % Generate fines of random diameter based on PSD
                F = zeros(numFinesParticle,7);
                F(:,2) = PSD(ceil(length(PSD)*rand(numFinesParticle,1)));
                F(:,3) = pi*(F(:,2).^2);
                F(:,4) = 1/6*pi*(F(:,2).^3);
                
                % Introduce a particles in each node
                idxInlet = find(nodeInlet==1);
                isFull = zeros(size(idxInlet));
                particleNum = 0;
                b = 1;
                while sum(isFull==0)
                    % Randomly chooses a node that is not full
                    idxInletCurrent = idxInlet(isFull==0);
                    nodeTo = idxInletCurrent(ceil(rand(1)*length(idxInletCurrent)));
                    
                    if nodeActive(nodeTo,3)<(nodeVolume(nodeTo)*densityLimit*contentLimit)
                        % If there is space, add particle
                        F(b,1) = nodeTo;
                        nodeActive(nodeTo,1) = nodeActive(nodeTo,1)+1;
                        nodeActive(nodeTo,3) = nodeActive(nodeTo,3)+F(b,4);
                        timeIntroduced(t,1) = timeIntroduced(t,1)+1;
                        timeIntroduced(t,2) = timeIntroduced(t,2)+F(b,4);
                        timeIntroduced(t,3) = max([timeIntroduced(t,3); F(b,2)]);
                        idx = find(F(b,2)<=psdFineBucket, 1, 'first');
                        timeIntroduced(t,idx+3) = timeIntroduced(t,idx+3)+1;
                        particleNum = particleNum+1;
                        if particleNum==numFinesParticle
                            break;
                        end
                        b = b+1;
                    else
                        % If there is no space, it becomes full
                        isFull(idxInlet==nodeTo) = 1;
                    end
                end
                if particleNum~=numFinesParticle
                    F = F(1:b-1,:);
                end
                finesNum = size(F,1);
                fprintf('Introduced %i/%i, ', particleNum, finesNum);
                
                % Store the nodes to which the base particles were introduced
                for i = 1:finesNum
                    particleRoute.(sprintf('particle%i',i)) = F(i,1);
                end
                idxParticle = (1:finesNum)';
            else
                % Generate fines of random diameter based on PSD
                currentNum = size(F, 1);
                
                % Initialise timeIntroduced
                timeIntroduced(t,:) = timeIntroduced(t-1,:);
                
                if t>maxFinesTimeStep
                    % Exceed maximumm time step to introduce fines
                    fprintf('Max t %i/%i, ', maxFinesTimeStep, currentNum);
                elseif currentNum>maxFinesParticle
                    % Exceed maximum number of particles
                    fprintf('Reached max %i/%i, ', maxFinesParticle, currentNum);
                else
                    % Generate fines of random diameter based on PSD
                    F = [F; zeros(numFinesParticle,7)];
                    F(currentNum+1:end,2) = PSD(ceil(10001*rand(numFinesParticle,1)));
                    F(currentNum+1:end,3) = pi*(F(currentNum+1:end,2).^2);
                    F(currentNum+1:end,4) = 1/6*pi*(F(currentNum+1:end,2).^3);
                    
                    % Introduce a maximum of n particles in each node
                    idxInlet = find(nodeInlet==1);
                    particleNum = 0;
                    isFull = zeros(size(idxInlet));
                    b = currentNum+1;
                    while sum(isFull==0)
                        % Randomly chooses a node that is not full
                        idxInletCurrent = idxInlet(isFull==0);
                        nodeTo = idxInletCurrent(ceil(rand(1)*length(idxInletCurrent)));
                        
                        if nodeActive(nodeTo,3)<(nodeVolume(nodeTo)*densityLimit*contentLimit)
                            % If there is space, add particle
                            F(b,1) = nodeTo;
                            nodeActive(nodeTo,1) = nodeActive(nodeTo,1)+1;
                            nodeActive(nodeTo,3) = nodeActive(nodeTo,3)+F(b,4);
                            timeIntroduced(t,1) = timeIntroduced(t,1)+1;
                            timeIntroduced(t,2) = timeIntroduced(t,2)+F(b,4);
                            timeIntroduced(t,3) = max([timeIntroduced(t,3); F(b,2)]);
                            idx = find(F(b,2)<=psdFineBucket, 1, 'first');
                            timeIntroduced(t,idx+3) = timeIntroduced(t,idx+3)+1;
                            particleNum = particleNum+1;
                            if particleNum==numFinesParticle
                                break;
                            end
                            b = b+1;
                        else
                            % If there is no space, it becomes full
                            isFull(idxInlet==nodeTo) = 1;
                        end
                    end
                    if particleNum~=numFinesParticle
                        F = F(1:b-1,:);
                    end
                    finesNum = size(F,1);
                    fprintf('Introduced %i/%i, ', particleNum, finesNum);
                    
                    % Update idxParticle
                    elementNum = numel(fieldnames(particleRoute));
                    for i = 1:finesNum-currentNum
                        particleRoute.(sprintf('particle%i',elementNum+i)) = F(currentNum+i,1);
                    end
                    idxParticle = [idxParticle; (currentNum+1:finesNum)'];
                end
            end
            
            % Stop when idxParticle reaches storeNum
            if t<maxFinesTimeStep
                if length(idxParticle)==storeNum
                    maxFinesTimeStep = t;
                end
            end
        else
            % storeNum of particles
            if isempty(idxParticle)
                numFinesParticle = storeNum;
            else
                numFinesParticle = storeNum-idxParticle;
            end
            
            % Place fines at the face where lattice z = 1
            if t==1
                % Generate fines of random diameter based on PSD
                F = zeros(numFinesParticle,7);
                F(:,2) = PSD(ceil(length(PSD)*rand(numFinesParticle,1)));
                F(:,3) = pi*(F(:,2).^2);
                F(:,4) = 1/6*pi*(F(:,2).^3);
                
                % Introduce a particles in each node
                idxInlet = find(nodeInlet==1);
                isFull = zeros(size(idxInlet));
                particleNum = 0;
                b = 1;
                while sum(isFull==0)
                    % Randomly chooses a node that is not full
                    idxInletCurrent = idxInlet(isFull==0);
                    nodeTo = idxInletCurrent(ceil(rand(1)*length(idxInletCurrent)));
                    
                    if nodeActive(nodeTo,3)<(nodeVolume(nodeTo)*densityLimit*contentLimit)
                        % If there is space, add particle
                        F(b,1) = nodeTo;
                        nodeActive(nodeTo,1) = nodeActive(nodeTo,1)+1;
                        nodeActive(nodeTo,3) = nodeActive(nodeTo,3)+F(b,4);
                        timeIntroduced(t,1) = timeIntroduced(t,1)+1;
                        timeIntroduced(t,2) = timeIntroduced(t,2)+F(b,4);
                        timeIntroduced(t,3) = max([timeIntroduced(t,3); F(b,2)]);
                        idx = find(F(b,2)<=psdFineBucket, 1, 'first');
                        timeIntroduced(t,idx+3) = timeIntroduced(t,idx+3)+1;
                        particleNum = particleNum+1;
                        if particleNum==numFinesParticle
                            break;
                        end
                        b = b+1;
                    else
                        % If there is no space, it becomes full
                        isFull(idxInlet==nodeTo) = 1;
                    end
                end
                if particleNum~=numFinesParticle
                    F = F(1:b-1,:);
                end
                finesNum = size(F,1);
                fprintf('Introduced %i/%i, ', particleNum, finesNum);
                idxParticle = finesNum;
            else
                % Generate fines of random diameter based on PSD
                currentNum = size(F,1);
                
                % Initialise timeIntroduced
                timeIntroduced(t,:) = timeIntroduced(t-1,:);
                
                if t>maxFinesTimeStep
                    % Exceed maximumm time step to introduce fines
                    fprintf('Max t %i/%i, ', maxFinesTimeStep, currentNum);
                elseif currentNum>maxFinesParticle
                    % Exceed maximum number of particles
                    fprintf('Reached max %i/%i, ', maxFinesParticle, currentNum);
                else
                    % Generate fines of random diameter based on PSD
                    F = [F; zeros(numFinesParticle,7)];
                    F(currentNum+1:end,2) = PSD(ceil(10001*rand(numFinesParticle,1)));
                    F(currentNum+1:end,3) = pi*(F(currentNum+1:end,2).^2);
                    F(currentNum+1:end,4) = 1/6*pi*(F(currentNum+1:end,2).^3);
                    
                    % Introduce a maximum of n particles in each node
                    idxInlet = find(nodeInlet==1);
                    particleNum = 0;
                    isFull = zeros(size(idxInlet));
                    b = currentNum+1;
                    while sum(isFull==0)
                        % Randomly chooses a node that is not full
                        idxInletCurrent = idxInlet(isFull==0);
                        nodeTo = idxInletCurrent(ceil(rand(1)*length(idxInletCurrent)));
                        
                        if nodeActive(nodeTo,3)<(nodeVolume(nodeTo)*densityLimit*contentLimit)
                            % If there is space, add particle
                            F(b,1) = nodeTo;
                            nodeActive(nodeTo,1) = nodeActive(nodeTo,1)+1;
                            nodeActive(nodeTo,3) = nodeActive(nodeTo,3)+F(b,4);
                            timeIntroduced(t,1) = timeIntroduced(t,1)+1;
                            timeIntroduced(t,2) = timeIntroduced(t,2)+F(b,4);
                            timeIntroduced(t,3) = max([timeIntroduced(t,3); F(b,2)]);
                            idx = find(F(b,2)<=psdFineBucket, 1, 'first');
                            timeIntroduced(t,idx+3) = timeIntroduced(t,idx+3)+1;
                            particleNum = particleNum+1;
                            if particleNum==numFinesParticle
                                break;
                            end
                            b = b+1;
                        else
                            % If there is no space, it becomes full
                            isFull(idxInlet==nodeTo) = 1;
                        end
                    end
                    if particleNum~=numFinesParticle
                        F = F(1:b-1,:);
                    end
                    finesNum = size(F,1);
                    fprintf('Introduced %i/%i, ', particleNum, finesNum);
                    idxParticle = idxParticle+particleNum;
                end
            end
            
            % Stop when idxParticle reaches storeNum
            if t<maxFinesTimeStep
                if idxParticle==storeNum
                    maxFinesTimeStep = t;
                end
            end
        end
    case 'refillInlet'
        if strcmp(storeOption, 'storeEdgePath')
            % A maximum of particleInlet particles at each node
            numFinesParticle = sum(nodeInlet)*particleInlet;
            
            % Place fines at the face where lattice z = 1
            if t==1
                % Generate fines of random diameter based on PSD
                F = zeros(numFinesParticle,7);
                F(:,2) = PSD(ceil(length(PSD)*rand(numFinesParticle,1)));
                F(:,3) = pi*(F(:,2).^2);
                F(:,4) = 1/6*pi*(F(:,2).^3);
                
                % Introduce particles in each node
                idxInlet = find(nodeInlet==1);
                isFull = zeros(size(idxInlet));
                particleNum = 0;
                b = 1;
                while sum(isFull==0)
                    % Randomly chooses a node that is not full
                    idxInletCurrent = idxInlet(isFull==0);
                    nodeTo = idxInletCurrent(ceil(rand(1)*length(idxInletCurrent)));
                    
                    if nodeActive(nodeTo,3)<(nodeVolume(nodeTo)*densityLimit*contentLimit)
                        % If there is space, add particle
                        F(b,1) = nodeTo;
                        nodeActive(nodeTo,1) = nodeActive(nodeTo,1)+1;
                        nodeActive(nodeTo,3) = nodeActive(nodeTo,3)+F(b,4);
                        timeIntroduced(t,1) = timeIntroduced(t,1)+1;
                        timeIntroduced(t,2) = timeIntroduced(t,2)+F(b,4);
                        timeIntroduced(t,3) = max([timeIntroduced(t,3); F(b,2)]);
                        idx = find(F(b,2)<=psdFineBucket, 1, 'first');
                        timeIntroduced(t,idx+3) = timeIntroduced(t,idx+3)+1;
                        particleNum = particleNum+1;
                        if particleNum==numFinesParticle
                            break;
                        end
                        b = b+1;
                    else
                        % If there is no space, it becomes full
                        isFull(idxInlet==nodeTo) = 1;
                    end
                end
                if particleNum~=numFinesParticle
                    F = F(1:b-1,:);
                end
                finesNum = size(F,1);
                
                % Store the nodes to which the base particles were introduced
                for i = 1:finesNum
                    particleRoute.(sprintf('idx%i',i)) = zeros(edgeNum,1);
                end
                
                % Display
                fprintf('Introduced %i/%i, ', particleNum, finesNum);
            else
                % Current number of active fines
                currentNum = size(F,1);
                
                % Initialise timeIntroduced
                timeIntroduced(t,:) = timeIntroduced(t-1,:);
                
                if t>maxFinesTimeStep
                    % Exceed maximumm time step to introduce fines
                    fprintf('Max t %i/%i, ', maxFinesTimeStep, currentNum);
                elseif currentNum>maxFinesParticle
                    % Exceed maximum number of particles
                    fprintf('Reached max %i/%i, ', maxFinesParticle, currentNum);
                else
                    % Generate fines of random diameter based on PSD
                    F = [F; zeros(numFinesParticle,7)];
                    F(currentNum+1:end,2) = PSD(ceil(10001*rand(numFinesParticle,1)));
                    F(currentNum+1:end,3) = pi*(F(currentNum+1:end,2).^2);
                    F(currentNum+1:end,4) = 1/6*pi*(F(currentNum+1:end,2).^3);
                    
                    % Introduce a maximum of n particles in each node
                    idxInlet = find(nodeInlet==1);
                    particleNum = 0;
                    isFull = zeros(size(idxInlet));
                    b = currentNum+1;
                    while sum(isFull==0)
                        % Randomly chooses a node that is not full
                        idxInletCurrent = idxInlet(isFull==0);
                        nodeTo = idxInletCurrent(ceil(rand(1)*length(idxInletCurrent)));
                        
                        if nodeActive(nodeTo,3)<(nodeVolume(nodeTo)*densityLimit*contentLimit)
                            % If there is space, add particle
                            F(b,1) = nodeTo;
                            nodeActive(nodeTo,1) = nodeActive(nodeTo,1)+1;
                            nodeActive(nodeTo,3) = nodeActive(nodeTo,3)+F(b,4);
                            timeIntroduced(t,1) = timeIntroduced(t,1)+1;
                            timeIntroduced(t,2) = timeIntroduced(t,2)+F(b,4);
                            timeIntroduced(t,3) = max([timeIntroduced(t,3); F(b,2)]);
                            idx = find(F(b,2)<=psdFineBucket, 1, 'first');
                            timeIntroduced(t,idx+3) = timeIntroduced(t,idx+3)+1;
                            particleNum = particleNum+1;
                            if particleNum==numFinesParticle
                                break;
                            end
                            b = b+1;
                        else
                            % If there is no space, it becomes full
                            isFull(idxInlet==nodeTo) = 1;
                        end
                    end
                    if particleNum~=numFinesParticle
                        F = F(1:b-1,:);
                    end
                    finesNum = size(F,1);
                    
                    % Store the nodes to which the base particles were introduced
                    for i = currentNum+1:finesNum
                        particleRoute.(sprintf('idx%i',i)) = zeros(edgeNum,1);
                    end
                    
                    % Display
                    fprintf('Introduced %i/%i, ', particleNum, finesNum);
                end
            end
        else
            % A maximum of particleInlet particles at each node
            numFinesParticle = sum(nodeInlet)*particleInlet;
            
            % Place fines at the face where lattice z = 1
            if t==1
                % Generate fines of random diameter based on PSD
                F = zeros(numFinesParticle,7);
                F(:,2) = PSD(ceil(length(PSD)*rand(numFinesParticle,1)));
                F(:,3) = pi*(F(:,2).^2);
                F(:,4) = 1/6*pi*(F(:,2).^3);
                
                % Introduce particles in each node
                idxInlet = find(nodeInlet==1);
                isFull = zeros(size(idxInlet));
                particleNum = 0;
                b = 1;
                while sum(isFull==0)
                    % Randomly chooses a node that is not full
                    idxInletCurrent = idxInlet(isFull==0);
                    nodeTo = idxInletCurrent(ceil(rand(1)*length(idxInletCurrent)));
                    
                    if nodeActive(nodeTo,3)<(nodeVolume(nodeTo)*densityLimit*contentLimit)
                        % If there is space, add particle
                        F(b,1) = nodeTo;
                        nodeActive(nodeTo,1) = nodeActive(nodeTo,1)+1;
                        nodeActive(nodeTo,3) = nodeActive(nodeTo,3)+F(b,4);
                        timeIntroduced(t,1) = timeIntroduced(t,1)+1;
                        timeIntroduced(t,2) = timeIntroduced(t,2)+F(b,4);
                        timeIntroduced(t,3) = max([timeIntroduced(t,3); F(b,2)]);
                        idx = find(F(b,2)<=psdFineBucket, 1, 'first');
                        timeIntroduced(t,idx+3) = timeIntroduced(t,idx+3)+1;
                        particleNum = particleNum+1;
                        if particleNum==numFinesParticle
                            break;
                        end
                        b = b+1;
                    else
                        % If there is no space, it becomes full
                        isFull(idxInlet==nodeTo) = 1;
                    end
                end
                if particleNum~=numFinesParticle
                    F = F(1:b-1,:);
                end
                finesNum = size(F,1);
                fprintf('Introduced %i/%i, ', particleNum, finesNum);
            else
                % Generate fines of random diameter based on PSD
                currentNum = size(F, 1);
                
                % Initialise timeIntroduced
                timeIntroduced(t,:) = timeIntroduced(t-1,:);
                
                if t>maxFinesTimeStep
                    % Exceed maximumm time step to introduce fines
                    fprintf('Max t %i/%i, ', maxFinesTimeStep, currentNum);
                elseif currentNum>maxFinesParticle
                    % Exceed maximum number of particles
                    fprintf('Reached max %i/%i, ', maxFinesParticle, currentNum);
                else
                    % Generate fines of random diameter based on PSD
                    F = [F; zeros(numFinesParticle,7)];
                    F(currentNum+1:end,2) = PSD(ceil(10001*rand(numFinesParticle,1)));
                    F(currentNum+1:end,3) = pi*(F(currentNum+1:end,2).^2);
                    F(currentNum+1:end,4) = 1/6*pi*(F(currentNum+1:end,2).^3);
                    
                    % Introduce a maximum of n particles in each node
                    idxInlet = find(nodeInlet==1);
                    particleNum = 0;
                    isFull = zeros(size(idxInlet));
                    b = currentNum+1;
                    while sum(isFull==0)
                        % Randomly chooses a node that is not full
                        idxInletCurrent = idxInlet(isFull==0);
                        nodeTo = idxInletCurrent(ceil(rand(1)*length(idxInletCurrent)));
                        
                        if nodeActive(nodeTo,3)<(nodeVolume(nodeTo)*densityLimit*contentLimit)
                            % If there is space, add particle
                            F(b,1) = nodeTo;
                            nodeActive(nodeTo,1) = nodeActive(nodeTo,1)+1;
                            nodeActive(nodeTo,3) = nodeActive(nodeTo,3)+F(b,4);
                            timeIntroduced(t,1) = timeIntroduced(t,1)+1;
                            timeIntroduced(t,2) = timeIntroduced(t,2)+F(b,4);
                            timeIntroduced(t,3) = max([timeIntroduced(t,3); F(b,2)]);
                            idx = find(F(b,2)<=psdFineBucket, 1, 'first');
                            timeIntroduced(t,idx+3) = timeIntroduced(t,idx+3)+1;
                            particleNum = particleNum+1;
                            if particleNum==numFinesParticle
                                break;
                            end
                            b = b+1;
                        else
                            % If there is no space, it becomes full
                            isFull(idxInlet==nodeTo) = 1;
                        end
                    end
                    if particleNum~=numFinesParticle
                        F = F(1:b-1,:);
                    end
                    finesNum = size(F,1);
                    fprintf('Introduced %i/%i, ', particleNum, finesNum);
                end
            end
        end
    case 'fixedParticleNum2'
        % A maximum of particleInlet particles at each node
        numFinesParticle = sum(nodeInlet)*particleInlet;
        
        % Generate fines of random diameter based on PSD
        currentNum = size(F, 1);
        
        if t<=maxFinesTimeStep
            % Place fines at the face where lattice z = 1
            if t==1
                % Generate fines of random diameter based on PSD
                F = zeros(numFinesParticle,7);
                F(:,2) = PSD(ceil(length(PSD)*rand(numFinesParticle,1)));
                F(:,3) = pi*(F(:,2).^2);
                F(:,4) = 1/6*pi*(F(:,2).^3);
            else
                % Initialise timeIntroduced
                timeIntroduced(t,:) = timeIntroduced(t-1,:);
                timeOutside(t,:) = timeOutside(t-1,:);
                
                % Generate fines of random diameter based on PSD
                F = [F; zeros(numFinesParticle,7)];
                F(currentNum+1:end,2) = PSD(ceil(10001*rand(numFinesParticle,1)));
                F(currentNum+1:end,3) = pi*(F(currentNum+1:end,2).^2);
                F(currentNum+1:end,4) = 1/6*pi*(F(currentNum+1:end,2).^3);
            end
            
            % Introduce particles in each node
            idxInlet = find(nodeInlet==1);
            inletNum = length(idxInlet);
            for i = currentNum+1:currentNum+numFinesParticle
                % Randomly chooses a node
                nodeTo = idxInlet(ceil(rand(1)*inletNum));
                
                if F(i,4)<(nodeVolume(nodeTo)*densityLimit)-nodeActive(nodeTo,3)
                    % If there is space, add particle
                    F(i,1) = nodeTo;
                    nodeActive(nodeTo,1) = nodeActive(nodeTo,1)+1;
                    nodeActive(nodeTo,3) = nodeActive(nodeTo,3)+F(i,4);
                    timeIntroduced(t,1) = timeIntroduced(t,1)+1;
                    timeIntroduced(t,2) = timeIntroduced(t,2)+F(i,4);
                    timeIntroduced(t,3) = max([timeIntroduced(t,3); F(i,2)]);
                    idx = find(F(i,2)<=psdFineBucket, 1, 'first');
                    timeIntroduced(t,idx+3) = timeIntroduced(t,idx+3)+1;
                else
                    % If there is no space, store it in timeOutside
                    F(i,5) = -2;
                    timeOutside(t,1) = timeOutside(t,1)+1;
                    timeOutside(t,2) = timeOutside(t,2)+F(i,4);
                    timeOutside(t,3) = max([timeOutside(t,3); F(i,2)]);
                    idx = find(F(i,2)<=psdFineBucket, 1, 'first');
                    timeOutside(t,idx+3) = timeOutside(t,idx+3)+1;
                end
            end
            
            % Remove particles outside
            F = F(F(:,5)==0,:);
            
            % Display
            finesNum = size(F,1);
            fprintf('Introduced %i/%i, ', finesNum-currentNum, finesNum);
        else
            % Exceed maximumm time step to introduce fines
            fprintf('Max t %i/%i, ', maxFinesTimeStep, currentNum);
        end
end

end

% ===== Flow Solver =======================================================
function [edgeConductance] = setNetworkConductance(edgeFlowing, edgeList, nodeVolume, nodeSurfaceArea, edgeThroatArea, edgeLength, nodeFinesVolume, nodeFinesArea, edgeReduced)
% This function defines the local conductance of the network edge element
% based on the equivalent hydraulic radius model in Chareyre et al. (2012)
% - see reference below. This model relies on the node volume, node
% surface area and edge throat area, all of which are calculated in the
% function 'setNetworkGeometry'. We will also incorporate the presence of
% fines by deducting the volume of fines from the pore volume.

% Ref: Chareyre et al. (2012), Pore-scale Modeling of Viscous Flow and
% Induced Forces in Dense Sphere Packings, Transport in Porous Media, Vol.
% 92(2), p. 64-75

% Set parameters required for the local conductance model, including the
% dynamic viscosity and the shape factor of the constriction. A value of
% 0.5 represents a circular cross-section.
dynamicViscosity = Constants.dynamicViscosity;
shapeFactor = Constants.shapeFactor;

if isempty(edgeLength)
    edgeConductance = edgeThroatArea;
else
    % Considers the effects of base particles
    edgeVolume = nodeVolume-nodeFinesVolume;
    edgeVolume = edgeVolume(edgeList);
    edgeSurfaceArea = nodeSurfaceArea+nodeFinesArea;
    edgeSurfaceArea = edgeSurfaceArea(edgeList);
    edgeHydraulicRadius = sum(edgeVolume, 2)./sum(edgeSurfaceArea, 2);
    edgeConductance = (shapeFactor*edgeThroatArea.*(edgeHydraulicRadius.^2))./(dynamicViscosity*edgeLength);
    
    % Base particles do not affect the flow
    % edgeRadius = sqrt(edgeThroatArea/pi);
    % edgeConductance = (pi*edgeRadius.^4)./(8*dynamicViscosity*edgeLength);
end
edgeConductance(edgeFlowing==0) = edgeConductance(edgeFlowing==0).*edgeReduced(edgeFlowing==0);
end
function [nodePressure, nodeSource] = calculatePressure(flowMethod, nodeNum, edgeNum, nodeConnectivity, nodeInlet, nodeOutlet, nodeEdgeTo, nodeNodeTo, edgeConductance, pressureEntry, pressureExit, edgeList, nodeCoord, flowDirection, networkType)
switch flowMethod
    case 'matlab'
        [nodePressure, nodeSource] = calculatePressureMatlab(nodeNum, nodeCoord, nodeConnectivity, nodeInlet, nodeOutlet, nodeEdgeTo, nodeNodeTo, edgeConductance, pressureEntry, pressureExit, flowDirection, networkType);
    case 'python'
        generateCSV(nodeNum, edgeNum, edgeList, edgeConductance, nodeCoord, nodeInlet, nodeOutlet);
        calculatePressureOpenPNM(FileLocation.flowSolverTXT(), FileLocation.flowSolverCSV(), pressureEntry, pressureExit, min(nodeCoord(:,3)), max(nodeCoord(:,3)), 'z')
        [nodePressure] = readTXT(FileLocation.flowSolverTXT());
        nodeSource = NaN;
end
end
function [nodePressure, nodeSource] = calculatePressureMatlab(nodeNum, nodeCoord, nodeConnectivity, nodeInlet, nodeOutlet, nodeEdgeTo, nodeNodeTo, edgeConductance, pressureEntry, pressureExit, flowDirection, networkType)
% The function solves for the local pressures in all nodes

conductanceNum = sum(nodeConnectivity+1);
idx = zeros(conductanceNum, 2);
conductance = zeros(conductanceNum, 1);

% Assemble the conductance matrix
conductanceID = 1;
for i = 1:nodeNum
    n = nodeConnectivity(i);
    conductanceCurrent = edgeConductance(nodeEdgeTo(i,1:n));
    conductance(conductanceID:conductanceID+n) = [sum(conductanceCurrent);-conductanceCurrent];
    idx(conductanceID:conductanceID+n,1) = i;
    idx(conductanceID:conductanceID+n,2) = [i,nodeNodeTo(i,1:n)];
    conductanceID = conductanceID+n+1;
end
C = sparse(idx(:,1), idx(:,2), conductance);

% Solve for pressure vector
P = zeros(nodeNum,1);
if strcmp(networkType, 'actual')
    switch flowDirection(1)
        case 'x'
            i = 1;
        case 'y'
            i = 2;
        case 'z'
            i = 3;
    end
    minCoord = min(nodeCoord(:,i));
    maxCoord = max(nodeCoord(:,i));
    nodeFraction = (nodeCoord(:,i)-minCoord)/(maxCoord-minCoord);
    if strcmp(flowDirection(2:4),'neg')
        nodeFraction = 1-nodeFraction;
    end
    P(nodeInlet==1) = pressureEntry-(pressureEntry-pressureExit)*nodeFraction(nodeInlet==1);
    P(nodeOutlet==1) = pressureEntry-(pressureEntry-pressureExit)*nodeFraction(nodeOutlet==1);
else
    P(nodeInlet==1) = pressureEntry;
    P(nodeOutlet==1) = pressureExit;
end
A = -C((nodeInlet+nodeOutlet)==0,nodeInlet==1)*P(nodeInlet==1)-C((nodeInlet+nodeOutlet)==0,nodeOutlet==1)*P(nodeOutlet==1);
P((nodeInlet+nodeOutlet)==0) = C((nodeInlet+nodeOutlet)==0,(nodeInlet+nodeOutlet)==0)\A;

% Solve for source vector
S = zeros(nodeNum,1);
S(nodeInlet==1) = C(nodeInlet==1,:)*P;
S(nodeOutlet==1) = C(nodeOutlet==1,:)*P;

% Output
nodePressure = P;
nodeSource = S;
end
function [] = generateCSV(nodeNum, edgeNum, edgeList, edgeConductance, nodeCoord, nodeInlet, nodeOutlet)
% The function creates a CSV file to be used as input to the Stokes flow
% solver with openPNM

edgeList = edgeList-1;

% Create header
data = cell(edgeNum+1,6);
header = {'throat.conns', 'throat.equivalent_conductance', 'pore.id', 'pore.coords', 'pore.inlets', 'pore.outlets'};
for i = 1:6
    data{1,i} = header{i};
end

% Assemble for the edge connectivity and edge conductance
for i = 1:edgeNum
    data{i+1,1} = sprintf('%i %i',edgeList(i,1), edgeList(i,2));
    data{i+1,2} = sprintf('%5f', edgeConductance(i));
end

% Assemble for the node id, coord, inlets and outlets
for i = 1:nodeNum
    data{i+1,3} = i;
    data{i+1,4} = sprintf('%.4f %.4f %.4f', nodeCoord(i,1), nodeCoord(i,2), nodeCoord(i,3));
    data{i+1,5} = nodeInlet(i);
    data{i+1,6} = nodeOutlet(i);
end

% Write to CSV file
T = table(data);
writetable(T, FileLocation.flowSolverCSV(), 'WriteVariableNames', false);
end
function [] = calculatePressureOpenPNM(file_name, csv_path, inletPressure, outletPressure, inletLocation, outletLocation, flowDirection)
% The function calls for a python method

% Add the current folder to the Python search path
if count(py.sys.path,['..', filesep, 'src']) == 0
    insert(py.sys.path,int32(0),['..', filesep, 'src']);
end

% Run the file
py.stokes_flow.stokes_flow(file_name, csv_path, 'delaunay', inletPressure, outletPressure, inletLocation, outletLocation, flowDirection, 1e-3, 298);
end
function [nodePressure] = readTXT(file_name)
% The function reads the node pressures calculated with openPNM
fileID = fopen(file_name, 'r');
nodePressure = fscanf(fileID, '%f');
fclose(fileID);
end

% ===== Update Volume and Surface Area in Nodes ===========================
function [nodeActive] = updateNodeActive(nodeNum, nodeActive, nodeRetained, F)
% Update the surface area and the volume of active fines in each node

% The matrix nodeActive stores:
% 1: Total number of particles active in the node
% 2: Total surface area of particles active or retained in the node
% 3: Total volume of particles active or retained in the node

% The matrix nodeRetained stores:
% 1: Total number of particles retained in the node
% 2: Total surface area of particles retained in the node
% 3: Total volume of particles retained in the node
% 4: Total number of particles that clog through Mechanism 1
% 5: Total number of particles that clog through Mechanism 3
% 6: If the node is clogged (for 'flowbiasedextreme')

if ~isempty(F)
    % If there is active particles
    nodeActive(:,2) = determineFinesArea(F, nodeNum);
    nodeActive(:,3) = determineFinesVolume(F, nodeNum);
    nodeActive(:,2) = nodeActive(:,2)+nodeRetained(:,2);
    nodeActive(:,3) = nodeActive(:,3)+nodeRetained(:,3);
else
    % If there is no active particle
    nodeActive(:,2) = nodeRetained(:,2);
    nodeActive(:,3) = nodeRetained(:,3);
end
end
function [nodeFinesVolume] = determineFinesVolume(F, nodeNum)
% Determine the volume of fines in each node
finesVolume = zeros(nodeNum, 1);
idxActive = find(F(:,5)>=0);
for i = 1:length(idxActive)
    j = idxActive(i);
    finesVolume(F(j,1)) = finesVolume(F(j,1))+ F(j,4);
end
nodeFinesVolume = finesVolume;
end
function [nodeFinesArea] = determineFinesArea(F, nodeNum)
% Determine the area of fines in each node
finesArea = zeros(nodeNum, 1);
idxActive = find(F(:,5)>=0);
for i = 1:length(idxActive)
    j = idxActive(i);
    finesArea(F(j,1)) = finesArea(F(j,1))+ F(j,3);
end
nodeFinesArea = finesArea;
end

% ===== Transport Logic ===================================================
function [nodeActive, nodeClogging, edgeFlowing, edgePath, particleRoute, F, finesNum] = transportLogic(nodeActive, nodeClogging, nodeConnectivity, nodeEdgeTo, nodeNodeTo, nodeEdgeToPos, nodePressure, nodeOutlet, nodeRetained, nodeVolume, edgeArea, edgeConductance, edgeDiameter, edgeFlowing, edgeLength, edgePath, particleRoute, F, transportLogicType, densityLimit, storeOption, t)
% The function performs transport logic

% The matrix nodeActive stores:
% 1: Total number of particles active in the node
% 2: Total surface area of particles active or retained in the node
% 3: Total volume of particles active or retained in the node

% The matrix F stores:
% 1: Node index
% 2: Particle diameter
% 3: Particle surface area
% 4: Particle volume
% 5: Particle status (-1:eroded, 0:active, 1:clogged, 2:retained, 3:nospace)

% The number of active fines
finesNum = size(F,1);

storeEdgePath = strcmp(storeOption, 'storeEdgePath');
storeNodeClogging = strcmp(storeOption, 'storeNodeClogging');

switch transportLogicType
    case 'flowbiased'
        % There are a set of rules governing the transport of fines:
        % 1. Every fine particle gets transported per time step.
        % 2. Fine travels throught the edge based on the concept of flow-biased
        % probability
        % 3. If throat area is smaller, the edge becomes clogged
        % 4. If there is no space on the other end, and the particles in
        % the opposite node are inactive
        % 5. If both of the constraints above are met, move
        
        % Loop through the active fine particles
        edgeFlowingCurrent = edgeFlowing;
        for i = 1:finesNum
            nodeFrom = F(i,1);
            
            % Find the edges and nodes connected to the reference node, as well
            % as the outflows
            n = nodeConnectivity(nodeFrom);
            edgeTo = nodeEdgeTo(nodeFrom,1:n)';
            nodeTo = nodeNodeTo(nodeFrom,1:n)';
            outflow = Simulator.determineOutflow(nodeFrom, nodeTo, edgeTo, edgeConductance, nodePressure);
            
            % Determine whether the node is flowing
            isFlowing = (outflow.*edgeFlowing(edgeTo,1))>0;
            isFlowingCurrent = (outflow.*edgeFlowingCurrent(edgeTo,1))>0;
            
            % Check if the node is flowing before this time step
            if sum(isFlowing)
                % Check if the node is still flowing within this time step
                if sum(isFlowingCurrent)
                    % Apply flow-biased probability
                    outflow = outflow(isFlowing==1);
                    nodeTo = nodeTo(isFlowing==1);
                    edgeTo = edgeTo(isFlowing==1);
                    outflowRatio = zeros(size(outflow));
                    outflowSum = sum(outflow);
                    for j = 1:length(outflow)
                        outflowRatio(j) = sum(outflow(1:j))/outflowSum;
                    end
                    idxChosen = length(outflow)+1-sum(outflowRatio>=rand(1));
                    nodeTo = nodeTo(idxChosen);
                    edgeTo = edgeTo(idxChosen);
                    
                    if edgeFlowingCurrent(edgeTo,1)==1
                        if edgeDiameter(edgeTo)<F(i,2)
                            % MECHANISM 1
                            % If diameter is smaller, the edge is clogged
                            F(i,5) = 1;
                            edgeFlowingCurrent(edgeTo,:) = [0, 1, nodeFrom, F(i,2), t];
                            nodeActive(nodeFrom,1) = nodeActive(nodeFrom,1)-1;
                            F(i,6) = F(i,6)+1;
                        elseif (nodeVolume(nodeTo)*densityLimit-nodeActive(nodeTo,3))>F(i,4)
                            % If there is space on the other node
                            edgePath(edgeTo) = edgePath(edgeTo)+1;
                            if storeEdgePath
                                particleRoute.(sprintf('idx%i',i))(edgeTo) = particleRoute.(sprintf('idx%i',i))(edgeTo)+1;
                            end
                            nodeActive(nodeFrom,3) = nodeActive(nodeFrom,3)-F(i,4);
                            nodeActive(nodeFrom,1) = nodeActive(nodeFrom,1)-1;
                            F(i,1) = nodeTo;
                            F(i,6) = F(i,6)+1;
                            F(i,7) = F(i,7)+edgeLength(edgeTo);
                            if nodeOutlet(nodeTo)
                                % ERODED
                                % Particle is eroded if it reaches outlet
                                F(i,5) = -1;
                            else
                                % MOVEMENT
                                nodeActive(nodeTo,3) = nodeActive(nodeTo,3)+F(i,4);
                                nodeActive(nodeTo,1) = nodeActive(nodeTo,1)+1;
                            end
                        elseif ~nodeActive(nodeTo,1)
                            % MECHANISM 3
                            % If opposite node is fully filled with
                            % inactive particles
                            F(i,5) = 3;
                            edgeFlowingCurrent(edgeTo,:) = [0, 3, nodeFrom, F(i,2), t];
                            nodeActive(nodeFrom,1) = nodeActive(nodeFrom,1)-1;
                            F(i,6) = F(i,6)+1;
                            
                            % Store the time step at which a node is full,
                            % the number of exits, number of exits clogged
                            % due to particle straining, number of exits
                            % clogged due to node clogging
                            if storeNodeClogging
                                if isnan(nodeClogging(nodeTo,1))
                                    nodeClogging(nodeTo,1) = t;
                                    
                                    n = nodeConnectivity(nodeTo);
                                    edgeToNodeTo = nodeEdgeTo(nodeTo,1:n)';
                                    nodeToNodeTo = nodeNodeTo(nodeTo,1:n)';
                                    outflowNodeTo = Simulator.determineOutflow(nodeTo, nodeToNodeTo, edgeToNodeTo, edgeConductance, nodePressure);
                                    idxExits = edgeToNodeTo(outflowNodeTo>0);
                                    nodeClogging(nodeTo,2) = length(idxExits);
                                    
                                    nodeClogging(nodeTo,3) = sum(edgeFlowing(idxExits,2)==1);
                                    nodeClogging(nodeTo,4) = sum(edgeFlowing(idxExits,2)==3);
                                end
                            end
                        end
                    end
                end
            else
                % MECHANISM 2
                % The particle is retained in a clogged node
                F(i,5) = 2;
                nodeActive(nodeFrom,1) = nodeActive(nodeFrom,1)-1;
            end
        end
        
        % Update clogged edge
        edgeFlowing = edgeFlowingCurrent;
    case 'flowbiasedpos'
        % There are a set of rules governing the transport of fines:
        % 1. Every fine particle gets transported per time step.
        % 2. Fine travels throught the edge based on the concept of flow-biased
        % probability (only the edges in the flow direction)
        % 3. If throat area is smaller, the edge becomes clogged
        % 4. If there is no space on the other end, and the particles in
        % the opposite node are inactive
        % 5. If both of the constraints above are met, move
        
        % Loop through the active fine particles
        edgeFlowingCurrent = edgeFlowing;
        for i = 1:finesNum
            nodeFrom = F(i,1);
            
            % Find the edges and nodes connected to the reference node, as well
            % as the outflows
            n = nodeConnectivity(nodeFrom);
            edgeTo = nodeEdgeTo(nodeFrom,1:n)';
            nodeTo = nodeNodeTo(nodeFrom,1:n)';
            edgeToPos = nodeEdgeToPos(nodeFrom,1:n)';
            outflow = Simulator.determineOutflow(nodeFrom, nodeTo, edgeTo, edgeConductance, nodePressure);
            
            % Determine whether the node is flowing
            isFlowing = (outflow.*edgeFlowing(edgeTo,1).*edgeToPos)>0;
            isFlowingCurrent = (outflow.*edgeFlowingCurrent(edgeTo,1).*edgeToPos)>0;
            
            % Check if the node is flowing before this time step
            if sum(isFlowing)
                % Check if the node is still flowing within this time step
                if sum(isFlowingCurrent)
                    % Apply flow-biased probability
                    outflow = outflow(isFlowing==1);
                    nodeTo = nodeTo(isFlowing==1);
                    edgeTo = edgeTo(isFlowing==1);
                    outflowRatio = zeros(size(outflow));
                    outflowSum = sum(outflow);
                    for j = 1:length(outflow)
                        outflowRatio(j) = sum(outflow(1:j))/outflowSum;
                    end
                    idxChosen = length(outflow)+1-sum(outflowRatio>=rand(1));
                    nodeTo = nodeTo(idxChosen);
                    edgeTo = edgeTo(idxChosen);
                    
                    if edgeFlowingCurrent(edgeTo,1)==1
                        if edgeDiameter(edgeTo)<F(i,2)
                            % MECHANISM 1
                            % If diameter is smaller, the edge is clogged
                            F(i,5) = 1;
                            edgeFlowingCurrent(edgeTo,:) = [0, 1, nodeFrom, F(i,2), t];
                            nodeActive(nodeFrom,1) = nodeActive(nodeFrom,1)-1;
                            F(i,6) = F(i,6)+1;
                        elseif (nodeVolume(nodeTo)*densityLimit-nodeActive(nodeTo,3))>F(i,4)
                            % If there is space on the other node
                            edgePath(edgeTo) = edgePath(edgeTo)+1;
                            if storeEdgePath
                                particleRoute.(sprintf('idx%i',i))(edgeTo) = particleRoute.(sprintf('idx%i',i))(edgeTo)+1;
                            end
                            nodeActive(nodeFrom,3) = nodeActive(nodeFrom,3)-F(i,4);
                            nodeActive(nodeFrom,1) = nodeActive(nodeFrom,1)-1;
                            F(i,1) = nodeTo;
                            F(i,6) = F(i,6)+1;
                            F(i,7) = F(i,7)+edgeLength(edgeTo);
                            if nodeOutlet(nodeTo)
                                % ERODED
                                % Particle is eroded if it reaches outlet
                                F(i,5) = -1;
                            else
                                % MOVEMENT
                                nodeActive(nodeTo,3) = nodeActive(nodeTo,3)+F(i,4);
                                nodeActive(nodeTo,1) = nodeActive(nodeTo,1)+1;
                            end
                        elseif ~nodeActive(nodeTo,1)
                            % MECHANISM 3
                            % If opposite node is fully filled
                            F(i,5) = 3;
                            edgeFlowingCurrent(edgeTo,:) = [0, 3, nodeFrom, F(i,2), t];
                            nodeActive(nodeFrom,1) = nodeActive(nodeFrom,1)-1;
                            F(i,6) = F(i,6)+1;
                            
                            % Store the time step at which a node is full,
                            % the number of exits, number of exits clogged
                            % due to particle straining, number of exits
                            % clogged due to node clogging
                            if storeNodeClogging
                                if isnan(nodeClogging(nodeTo,1))
                                    nodeClogging(nodeTo,1) = t;
                                    
                                    n = nodeConnectivity(nodeTo);
                                    edgeToNodeTo = nodeEdgeTo(nodeTo,1:n)';
                                    nodeToNodeTo = nodeNodeTo(nodeTo,1:n)';
                                    outflowNodeTo = Simulator.determineOutflow(nodeTo, nodeToNodeTo, edgeToNodeTo, edgeConductance, nodePressure);
                                    idxExits = edgeToNodeTo(outflowNodeTo>0);
                                    nodeClogging(nodeTo,2) = length(idxExits);
                                    
                                    nodeClogging(nodeTo,3) = sum(edgeFlowing(idxExits,2)==1);
                                    nodeClogging(nodeTo,4) = sum(edgeFlowing(idxExits,2)==3);
                                end
                            end
                        end
                    end
                end
            else
                % MECHANISM 2
                % The particle is retained in a clogged node
                F(i,5) = 2;
                nodeActive(nodeFrom,1) = nodeActive(nodeFrom,1)-1;
            end
        end
        
        % Update clogged edge
        edgeFlowing = edgeFlowingCurrent;
    case 'flowbiasednode'
        % There are a set of rules governing the transport of fines:
        % 1. Every fine particle gets transported per time step.
        % 2. Fine travels throught the edge based on the concept of flow-biased
        % probability
        % 3. If throat area is smaller, BOTH the edge and the node becomes clogged
        % 4. If there is no space on the other end, and the particles in
        % the opposite node are inactive
        % 5. If both of the constraints above are met, move
        
        % Loop through the active fine particles
        edgeFlowingCurrent = edgeFlowing;
        for i = 1:finesNum
            nodeFrom = F(i,1);
            
            % Find the edges and nodes connected to the reference node, as well
            % as the outflows
            n = nodeConnectivity(nodeFrom);
            edgeTo = nodeEdgeTo(nodeFrom,1:n)';
            nodeTo = nodeNodeTo(nodeFrom,1:n)';
            outflow = Simulator.determineOutflow(nodeFrom, nodeTo, edgeTo, edgeConductance, nodePressure);
            
            % Determine whether the node is flowing
            isFlowing = (outflow.*edgeFlowing(edgeTo,1))>0;
            isFlowingCurrent = (outflow.*edgeFlowingCurrent(edgeTo,1))>0;
            
            % Check if the node is clogged
            if ~nodeRetained(nodeFrom,6)
                % Check if the node is flowing before this time step
                if sum(isFlowing)
                    % Check if the node is still flowing within this time step
                    if sum(isFlowingCurrent)
                        % Apply flow-biased probability
                        outflow = outflow(isFlowing==1);
                        nodeTo = nodeTo(isFlowing==1);
                        edgeTo = edgeTo(isFlowing==1);
                        outflowRatio = zeros(size(outflow));
                        outflowSum = sum(outflow);
                        for j = 1:length(outflow)
                            outflowRatio(j) = sum(outflow(1:j))/outflowSum;
                        end
                        idxChosen = length(outflow)+1-sum(outflowRatio>=rand(1));
                        nodeTo = nodeTo(idxChosen);
                        edgeTo = edgeTo(idxChosen);
                        
                        if edgeFlowingCurrent(edgeTo,1)==1
                            if edgeDiameter(edgeTo)<F(i,2)
                                % MECHANISM 1
                                % If diameter is smaller, the edge is clogged
                                F(i,5) = 1;
                                edgeFlowingCurrent(edgeTo,:) = [0, 1, nodeFrom, F(i,2), t];
                                nodeActive(nodeFrom,1) = nodeActive(nodeFrom,1)-1;
                                nodeRetained(nodeFrom,6) = 1;
                                F(i,6) = F(i,6)+1;
                            elseif (nodeVolume(nodeTo)*densityLimit-nodeActive(nodeTo,3))>F(i,4)
                                % If there is space on the other node
                                edgePath(edgeTo) = edgePath(edgeTo)+1;
                                if storeEdgePath
                                    particleRoute.(sprintf('idx%i',i))(edgeTo) = particleRoute.(sprintf('idx%i',i))(edgeTo)+1;
                                end
                                nodeActive(nodeFrom,3) = nodeActive(nodeFrom,3)-F(i,4);
                                nodeActive(nodeFrom,1) = nodeActive(nodeFrom,1)-1;
                                F(i,1) = nodeTo;
                                F(i,6) = F(i,6)+1;
                                F(i,7) = F(i,7)+edgeLength(edgeTo);
                                if nodeOutlet(nodeTo)
                                    % ERODED
                                    % Particle is eroded if it reaches outlet
                                    F(i,5) = -1;
                                else
                                    % MOVEMENT
                                    nodeActive(nodeTo,3) = nodeActive(nodeTo,3)+F(i,4);
                                    nodeActive(nodeTo,1) = nodeActive(nodeTo,1)+1;
                                end
                            elseif ~nodeActive(nodeTo,1)
                                % MECHANISM 3
                                % If opposite node is fully filled
                                F(i,5) = 3;
                                edgeFlowingCurrent(edgeTo,:) = [0, 3, nodeFrom, F(i,2), t];
                                nodeActive(nodeFrom,1) = nodeActive(nodeFrom,1)-1;
                                F(i,6) = F(i,6)+1;
                                
                                % Store the time step at which a node is full,
                                % the number of exits, number of exits clogged
                                % due to particle straining, number of exits
                                % clogged due to node clogging
                                if storeNodeClogging
                                    if isnan(nodeClogging(nodeTo,1))
                                        nodeClogging(nodeTo,1) = t;
                                        
                                        n = nodeConnectivity(nodeTo);
                                        edgeToNodeTo = nodeEdgeTo(nodeTo,1:n)';
                                        nodeToNodeTo = nodeNodeTo(nodeTo,1:n)';
                                        outflowNodeTo = Simulator.determineOutflow(nodeTo, nodeToNodeTo, edgeToNodeTo, edgeConductance, nodePressure);
                                        idxExits = edgeToNodeTo(outflowNodeTo>0);
                                        nodeClogging(nodeTo,2) = length(idxExits);
                                        
                                        nodeClogging(nodeTo,3) = sum(edgeFlowing(idxExits,2)==1);
                                        nodeClogging(nodeTo,4) = sum(edgeFlowing(idxExits,2)==3);
                                    end
                                end
                            end
                        end
                    end
                else
                    % MECHANISM 2
                    % The particle is retained in a clogged node
                    F(i,5) = 2;
                    nodeActive(nodeFrom,1) = nodeActive(nodeFrom,1)-1;
                end
            else
                % MECHANISM 2
                % The particle is retained in a clogged node
                F(i,5) = 2;
                nodeActive(nodeFrom,1) = nodeActive(nodeFrom,1)-1;
            end
        end
        
        % Update clogged edge
        edgeFlowing = edgeFlowingCurrent;
    case 'flowbiasedextreme'
        % Combination of flowbiasednode and flowbiasedpos
        
        % Loop through the active fine particles
        edgeFlowingCurrent = edgeFlowing;
        for i = 1:finesNum
            nodeFrom = F(i,1);
            
            % Find the edges and nodes connected to the reference node, as well
            % as the outflows
            n = nodeConnectivity(nodeFrom);
            edgeTo = nodeEdgeTo(nodeFrom,1:n)';
            nodeTo = nodeNodeTo(nodeFrom,1:n)';
            edgeToPos = nodeEdgeToPos(nodeFrom,1:n)';
            outflow = Simulator.determineOutflow(nodeFrom, nodeTo, edgeTo, edgeConductance, nodePressure);
            
            % Determine whether the node is flowing
            isFlowing = (outflow.*edgeFlowing(edgeTo,1).*edgeToPos)>0;
            isFlowingCurrent = (outflow.*edgeFlowingCurrent(edgeTo,1).*edgeToPos)>0;
            
            % Check if the node is clogged
            if ~nodeRetained(nodeFrom,6)
                % Check if the node is flowing before this time step
                if sum(isFlowing)
                    % Check if the node is still flowing within this time step
                    if sum(isFlowingCurrent)
                        % Apply flow-biased probability
                        outflow = outflow(isFlowing==1);
                        nodeTo = nodeTo(isFlowing==1);
                        edgeTo = edgeTo(isFlowing==1);
                        outflowRatio = zeros(size(outflow));
                        outflowSum = sum(outflow);
                        for j = 1:length(outflow)
                            outflowRatio(j) = sum(outflow(1:j))/outflowSum;
                        end
                        idxChosen = length(outflow)+1-sum(outflowRatio>=rand(1));
                        nodeTo = nodeTo(idxChosen);
                        edgeTo = edgeTo(idxChosen);
                        
                        if edgeFlowingCurrent(edgeTo,1)==1
                            if edgeDiameter(edgeTo)<F(i,2)
                                % MECHANISM 1
                                % If diameter is smaller, the edge is clogged
                                F(i,5) = 1;
                                edgeFlowingCurrent(edgeTo,:) = [0, 1, nodeFrom, F(i,2), t];
                                nodeActive(nodeFrom,1) = nodeActive(nodeFrom,1)-1;
                                nodeRetained(nodeFrom,6) = 1;
                                F(i,6) = F(i,6)+1;
                            elseif (nodeVolume(nodeTo)*densityLimit-nodeActive(nodeTo,3))>F(i,4)
                                % If there is space on the other node
                                edgePath(edgeTo) = edgePath(edgeTo)+1;
                                if storeEdgePath
                                    particleRoute.(sprintf('idx%i',i))(edgeTo) = particleRoute.(sprintf('idx%i',i))(edgeTo)+1;
                                end
                                nodeActive(nodeFrom,3) = nodeActive(nodeFrom,3)-F(i,4);
                                nodeActive(nodeFrom,1) = nodeActive(nodeFrom,1)-1;
                                F(i,1) = nodeTo;
                                F(i,6) = F(i,6)+1;
                                F(i,7) = F(i,7)+edgeLength(edgeTo);
                                if nodeOutlet(nodeTo)
                                    % ERODED
                                    % Particle is eroded if it reaches outlet
                                    F(i,5) = -1;
                                else
                                    % MOVEMENT
                                    nodeActive(nodeTo,3) = nodeActive(nodeTo,3)+F(i,4);
                                    nodeActive(nodeTo,1) = nodeActive(nodeTo,1)+1;
                                end
                            elseif ~nodeActive(nodeTo,1)
                                % MECHANISM 3
                                % If opposite node is fully filled
                                F(i,5) = 3;
                                edgeFlowingCurrent(edgeTo,:) = [0, 3, nodeFrom, F(i,2), t];
                                nodeActive(nodeFrom,1) = nodeActive(nodeFrom,1)-1;
                                F(i,6) = F(i,6)+1;
                                
                                % Store the time step at which a node is full,
                                % the number of exits, number of exits clogged
                                % due to particle straining, number of exits
                                % clogged due to node clogging
                                if storeNodeClogging
                                    if isnan(nodeClogging(nodeTo,1))
                                        nodeClogging(nodeTo,1) = t;
                                        
                                        n = nodeConnectivity(nodeTo);
                                        edgeToNodeTo = nodeEdgeTo(nodeTo,1:n)';
                                        nodeToNodeTo = nodeNodeTo(nodeTo,1:n)';
                                        outflowNodeTo = Simulator.determineOutflow(nodeTo, nodeToNodeTo, edgeToNodeTo, edgeConductance, nodePressure);
                                        idxExits = edgeToNodeTo(outflowNodeTo>0);
                                        nodeClogging(nodeTo,2) = length(idxExits);
                                        
                                        nodeClogging(nodeTo,3) = sum(edgeFlowing(idxExits,2)==1);
                                        nodeClogging(nodeTo,4) = sum(edgeFlowing(idxExits,2)==3);
                                    end
                                end
                            end
                        end
                    end
                else
                    % MECHANISM 2
                    % The particle is retained in a clogged node
                    F(i,5) = 2;
                    nodeActive(nodeFrom,1) = nodeActive(nodeFrom,1)-1;
                end
            else
                % MECHANISM 2
                % The particle is retained in a clogged node
                F(i,5) = 2;
                nodeActive(nodeFrom,1) = nodeActive(nodeFrom,1)-1;
            end
        end
        
        % Update clogged edge
        edgeFlowing = edgeFlowingCurrent;
    case 'noclogging'
        % There are a set of rules governing the transport of fines:
        % 1. Every fine particle gets transported per time step.
        % 2. Of the exit edges connecting a node, the fine travels
        % through the edge with the highest outflow.
        % 3. The particle must satisfy the throat area and pore volume
        % constraint.
        
        % Loop through the active fine particles
        for i = 1:finesNum
            nodeFrom = F(i,1);
            
            % Find the edges and nodes connected to the reference node, as well
            % as the outflows
            n = nodeConnectivity(nodeFrom);
            edgeTo = nodeEdgeTo(nodeFrom,1:n)';
            nodeTo = nodeNodeTo(nodeFrom,1:n)';
            outflow = Simulator.determineOutflow(nodeFrom, nodeTo, edgeTo, edgeConductance, nodePressure);
            
            % Create a data array containing the node id (of node
            % connected to reference node, nodeFrom), the corresponding
            % edge conductance and a boolean stating whether the throat
            % area and volume constraints have been satisfied.
            data = [nodeTo, ...
                outflow, ...
                and(and(edgeDiameter(edgeTo)>F(i,2), (densityLimit*nodeVolume(nodeTo) - nodeActive(nodeTo,3))>F(i,4)), outflow>0)];
            
            % Find the first (i.e. highest flow) edge that
            % satisfies the constraint conditions. If there are none, then
            % skip, otherwise move the particle
            if sum(data(:, 3))
                % Determine the index of the highest flow
                maxFlow = max(data(data(:,3)==1,2));
                idxMaxFlow = find(data(:,2)==maxFlow,1);
                nodeTo = data(idxMaxFlow,1);
                edgeTo = edgeTo(idxMaxFlow,1);
                
                % Move the particle
                edgePath(edgeTo) = edgePath(edgeTo)+1;
                if storeEdgePath
                    particleRoute.(sprintf('idx%i',i))(edgeTo) = particleRoute.(sprintf('idx%i',i))(edgeTo)+1;
                end
                nodeActive(nodeFrom,3) = nodeActive(nodeFrom,3)-F(i,4);
                nodeActive(nodeFrom,1) = nodeActive(nodeFrom,1)-1;
                F(i,1) = nodeTo;
                F(i,6) = F(i,6)+1;
                F(i,7) = F(i,7)+edgeLength(edgeTo);
                if nodeOutlet(nodeTo)
                    % Reached outlet
                    F(i,5) = -1;
                else
                    nodeActive(nodeTo,3) = nodeActive(nodeTo,3)+F(i,4);
                    nodeActive(nodeTo,1) = nodeActive(nodeTo,1)+1;
                end
            end
        end
    case 'noflow'
        % There are a set of rules governing the transport of fines:
        % 1. Every fine particle gets transported per time step.
        % 2. Fine travels throught the edge based on the concept of area-biased
        % probability
        % 3. If throat area is smaller, the edge becomes clogged
        % 4. Otherwise, move
        
        % Loop through the active fine particles
        for i = 1:finesNum
            nodeFrom = F(i,1);
            
            % Find the edges and nodes connected to the reference node, as well
            % as the outflows
            n = nodeConnectivity(nodeFrom);
            edgeTo = nodeEdgeTo(nodeFrom,1:n)';
            nodeTo = nodeNodeTo(nodeFrom,1:n)';
            throatArea = edgeArea(edgeTo);
            
            % Determine whether the node is flowing
            isFlowing = edgeFlowing(edgeTo,1)==1;
            
            % Skip if the node is not flowing
            if sum(isFlowing)
                % Apply area-biased probability
                throatArea = throatArea(isFlowing==1);
                nodeTo = nodeTo(isFlowing==1);
                edgeTo = edgeTo(isFlowing==1);
                randomRatio = zeros(size(throatArea));
                randomSum = sum(throatArea);
                for j = 1:length(throatArea)
                    randomRatio(j) = sum(throatArea(1:j))/randomSum;
                end
                idxChosen = length(throatArea)+1-sum(randomRatio>rand(1));
                nodeTo = nodeTo(idxChosen);
                edgeTo = edgeTo(idxChosen);
                
                if edgeDiameter(edgeTo)<F(i,2)
                    % If diameter is smaller, the edge is clogged
                    F(i,5) = 1;
                    edgeFlowing(edgeTo,:) = [0, 1, nodeFrom, F(i,2), t];
                    nodeActive(nodeFrom,1) = nodeActive(nodeFrom,1)-1;
                    F(i,6) = F(i,6)+1;
                elseif (nodeVolume(nodeTo)*densityLimit-nodeActive(nodeTo,3))>F(i,4)
                    % If there is space on the other node
                    edgePath(edgeTo) = edgePath(edgeTo)+1;
                    if storeEdgePath
                        particleRoute.(sprintf('idx%i',i))(edgeTo) = particleRoute.(sprintf('idx%i',i))(edgeTo)+1;
                    end
                    nodeActive(nodeFrom,3) = nodeActive(nodeFrom,3)-F(i,4);
                    nodeActive(nodeFrom,1) = nodeActive(nodeFrom,1)-1;
                    F(i,1) = nodeTo;
                    F(i,6) = F(i,6)+1;
                    F(i,7) = F(i,7)+edgeLength(edgeTo);
                    if nodeOutlet(nodeTo)
                        % ERODED
                        % Particle is eroded if it reaches outlet
                        F(i,5) = -1;
                    else
                        % MOVEMENT
                        nodeActive(nodeTo,3) = nodeActive(nodeTo,3)+F(i,4);
                        nodeActive(nodeTo,1) = nodeActive(nodeTo,1)+1;
                    end
                elseif ~nodeActive(nodeTo,1)
                    % MECHANISM 3
                    % If opposite node is fully filled
                    F(i,5) = 3;
                    edgeFlowing(edgeTo,:) = [0, 3, nodeFrom, F(i,2), t];
                    nodeActive(nodeFrom,1) = nodeActive(nodeFrom,1)-1;
                    F(i,6) = F(i,6)+1;
                    
                    % Store the time step at which a node is full,
                    % the number of exits, number of exits clogged
                    % due to particle straining, number of exits
                    % clogged due to node clogging
                    if storeNodeClogging
                        if isnan(nodeClogging(nodeTo,1))
                            nodeClogging(nodeTo,1) = t;
                            
                            n = nodeConnectivity(nodeTo);
                            idxExits = nodeEdgeTo(nodeTo,1:n)';
                            nodeClogging(nodeTo,2) = n;
                            
                            nodeClogging(nodeTo,3) = sum(edgeFlowing(idxExits,2)==1);
                            nodeClogging(nodeTo,4) = sum(edgeFlowing(idxExits,2)==3);
                        end
                    end
                end
            else
                % The particle is retained in a clogged node
                F(i,5) = 2;
                nodeActive(nodeFrom,1) = nodeActive(nodeFrom,1)-1;
            end
        end
    case 'noflowrandom'
        % There are a set of rules governing the transport of fines:
        % 1. Every fine particle gets transported per time step.
        % 2. Fine randomly chooses an edge to flow through
        % 3. If throat area is smaller, the edge becomes clogged
        % 4. Otherwise, move
        
        % Loop through the active fine particles
        for i = 1:finesNum
            nodeFrom = F(i,1);
            
            % Find the edges and nodes connected to the reference node, as well
            % as the outflows
            n = nodeConnectivity(nodeFrom);
            edgeTo = nodeEdgeTo(nodeFrom,1:n)';
            nodeTo = nodeNodeTo(nodeFrom,1:n)';
            
            % Determine whether the node is flowing
            isFlowing = edgeFlowing(edgeTo,1)==1;
            
            % Skip if the node is not flowing
            if sum(isFlowing)
                % Apply random-walk
                nodeTo = nodeTo(isFlowing==1);
                edgeTo = edgeTo(isFlowing==1);
                random = ones(size(edgeTo));
                randomRatio = zeros(size(random));
                randomSum = sum(random);
                for j = 1:length(random)
                    randomRatio(j) = sum(random(1:j))/randomSum;
                end
                idxChosen = length(random)+1-sum(randomRatio>rand(1));
                nodeTo = nodeTo(idxChosen);
                edgeTo = edgeTo(idxChosen);
                
                if edgeDiameter(edgeTo)<F(i,2)
                    % If diameter is smaller, the edge is clogged
                    F(i,5) = 1;
                    edgeFlowing(edgeTo,:) = [0, 1, nodeFrom, F(i,2), t];
                    nodeActive(nodeFrom,1) = nodeActive(nodeFrom,1)-1;
                    F(i,6) = F(i,6)+1;
                elseif (nodeVolume(nodeTo)*densityLimit-nodeActive(nodeTo,3))>F(i,4)
                    % If there is space on the other node
                    edgePath(edgeTo) = edgePath(edgeTo)+1;
                    if storeEdgePath
                        particleRoute.(sprintf('idx%i',i))(edgeTo) = particleRoute.(sprintf('idx%i',i))(edgeTo)+1;
                    end
                    nodeActive(nodeFrom,3) = nodeActive(nodeFrom,3)-F(i,4);
                    nodeActive(nodeFrom,1) = nodeActive(nodeFrom,1)-1;
                    F(i,1) = nodeTo;
                    F(i,6) = F(i,6)+1;
                    F(i,7) = F(i,7)+edgeLength(edgeTo);
                    if nodeOutlet(nodeTo)
                        % ERODED
                        % Particle is eroded if it reaches outlet
                        F(i,5) = -1;
                    else
                        % MOVEMENT
                        nodeActive(nodeTo,3) = nodeActive(nodeTo,3)+F(i,4);
                        nodeActive(nodeTo,1) = nodeActive(nodeTo,1)+1;
                    end
                elseif ~nodeActive(nodeTo,1)
                    % MECHANISM 3
                    % If opposite node is fully filled
                    F(i,5) = 3;
                    edgeFlowing(edgeTo,:) = [0, 3, nodeFrom, F(i,2), t];
                    nodeActive(nodeFrom,1) = nodeActive(nodeFrom,1)-1;
                    F(i,6) = F(i,6)+1;
                    
                    % Store the time step at which a node is full,
                    % the number of exits, number of exits clogged
                    % due to particle straining, number of exits
                    % clogged due to node clogging
                    if storeNodeClogging
                        if isnan(nodeClogging(nodeTo,1))
                            nodeClogging(nodeTo,1) = t;
                            
                            n = nodeConnectivity(nodeTo);
                            idxExits = nodeEdgeTo(nodeTo,1:n)';
                            nodeClogging(nodeTo,2) = n;
                            
                            nodeClogging(nodeTo,3) = sum(edgeFlowing(idxExits,2)==1);
                            nodeClogging(nodeTo,4) = sum(edgeFlowing(idxExits,2)==3);
                        end
                    end
                end
            else
                % The particle is retained in a clogged node
                F(i,5) = 2;
                nodeActive(nodeFrom,1) = nodeActive(nodeFrom,1)-1;
            end
        end
end

% Display progress
idxActive = find(F(:,5)==0);
finesNum = length(idxActive);
fprintf('%i active particles, t = %i\n', finesNum, t);
end

% ===== Finalise Arrays ===================================================
function [F, timeEroded, timeRetained, timeNodeRetained, nodeRetained, layerRetained, particlePathEroded, particlePathRetained, particleDistanceEroded, particleDistanceRetained, particleRoute, idxParticle, edgePathEroded, edgePathRetained] = rearrangeF(nodeNum, edgeNum, nodeConnectivity, nodeEdgeTo, nodeEdgeDirection, nodeInlet, nodeLattice, nodePressure, edgeConductance, edgeFlowing, edgeList, edgePathEroded, edgePathRetained, F, psdFineBucket, timeEroded, timeRetained, timeNodeRetained, nodeRetained, layerRetained, particlePathEroded, particlePathRetained, particleDistanceEroded, particleDistanceRetained, particleRoute, idxParticle, storeOption, t)
% The function records the data of the eroded, captured or retained
% base particles

% idxInlet
idxInlet = find(nodeInlet==1);

% Divide F
FEroded = F(F(:,5)==-1,:);
FRetained = F(F(:,5)>0,:);

if t~=1
    timeEroded(t,:) = timeEroded(t-1,:);
    timeRetained(t,:) = timeRetained(t-1,:);
end

% Record the data of eroded particles (timeEroded)
timeEroded(t,1) = timeEroded(t,1)+size(FEroded,1);
timeEroded(t,2) = timeEroded(t,2)+sum(FEroded(:,4));
timeEroded(t,3) = max([timeEroded(t,3); FEroded(:,2)]);
for i = 1:size(FEroded,1)
    idx = find(FEroded(i,2)<=psdFineBucket, 1, 'first');
    timeEroded(t,idx+3) = timeEroded(t,idx+3)+1;
    
    % Determine the number of constrictions met by eroded particles
    % (particlePathEroded)
    idx = FEroded(i,6)+1;
    if size(particlePathEroded,1)<idx
        particlePathEroded(idx) = 0;
    end
    particlePathEroded(idx) = particlePathEroded(idx)+1;
end

% Record the data of captured particles
for i = 1:size(FRetained,1)
    % The matrix (nodeRetained) stores:
    % 1: Total number of particles retained in the node
    % 2: Total surface area of particles retained in the node
    % 3: Total volume of particles retained in the node
    % 4: Total number of particles that clog through Mechanism 1
    % 5: Total number of particles that clog through Mechanism 3
    node = FRetained(i,1);
    nodeRetained(node,1) = nodeRetained(node,1)+1;
    nodeRetained(node,2) = nodeRetained(node,2)+FRetained(i,3);
    nodeRetained(node,3) = nodeRetained(node,3)+FRetained(i,4);
    if FRetained(i,5)==1
        nodeRetained(node,4) = nodeRetained(node,4)+1;
    end
    if FRetained(i,5)==3
        nodeRetained(node,5) = nodeRetained(node,5)+1;
    end
    
    if ~ismember(node, idxInlet)
        % (timeRetained)
        timeRetained(t,1) = timeRetained(t,1)+1;
        timeRetained(t,2) = timeRetained(t,2)+FRetained(i,4);
        timeRetained(t,3) = max([timeRetained(t,3); FRetained(i,2)]);
        idx = find(FRetained(i,2)<=psdFineBucket, 1, 'first');
        timeRetained(t,idx+3) = timeRetained(t,idx+3)+1;
        
        % Determine the PSD of captured particles in each layer
        % (layerRetainedPSD)
        layer = nodeLattice(FRetained(i,1),end);
        layerRetained(layer,idx) = layerRetained(layer,idx)+1;
        
        % Determine the number of constrictions met by retained particles
        % (particlePathRetained)
        idx = FRetained(i,6)+1;
        if size(particlePathRetained,1)<idx
            particlePathRetained(idx,:) = [0,0,0];
        end
        if FRetained(i,5)==1
            particlePathRetained(idx,1) = particlePathRetained(idx,1)+1;
        elseif FRetained(i,5)==2
            particlePathRetained(idx,2) = particlePathRetained(idx,2)+1;
        elseif FRetained(i,5)==3
            particlePathRetained(idx,3) = particlePathRetained(idx,3)+1;
        end
    end
end

% Store the distance travelled
idxNotInlet = ~ismember(FRetained(:,1), idxInlet);
particleDistanceEroded = [particleDistanceEroded; FEroded(:,7)];
particleDistanceRetained.mechanism1 = [particleDistanceRetained.mechanism1; FRetained(FRetained(idxNotInlet,5)==1,7)];
particleDistanceRetained.mechanism2 = [particleDistanceRetained.mechanism2; FRetained(FRetained(idxNotInlet,5)==2,7)];
particleDistanceRetained.mechanism3 = [particleDistanceRetained.mechanism3; FRetained(FRetained(idxNotInlet,5)==3,7)];

% More store options
if strcmp(storeOption, 'storePath')
    % Store the path taken by the particles
    particleNum = length(idxParticle);
    j = 0;
    for i = 1:particleNum
        particle = idxParticle(i);
        if particle>0
            if F(particle,5)==0
                % For the active particles
                if particleRoute.(sprintf('particle%i',i))(end)~=F(particle,1)
                    particleRoute.(sprintf('particle%i',i)) = [particleRoute.(sprintf('particle%i',i));F(particle,1)];
                end
            elseif F(particle,5)==-1
                % For the eroded particles
                particleRoute.(sprintf('particle%i',i)) = [particleRoute.(sprintf('particle%i',i));F(particle,1)];
                idxParticle(i) = 0;
            else
                % For the retained particles
                idxParticle(i) = 0;
            end
        end
        
        % Update idxParticle
        particle = idxParticle(i);
        if particle>0
            j = j+1;
            idxParticle(i) = j;
        end
    end
elseif strcmp(storeOption, 'more')
    % Initialise timeNodeRetained
    timeNodeRetained.(sprintf('t%i',t)) = zeros(nodeNum,4);
    
    % Determine edgeFlow
    [edgeFlow] = Simulator.determineEdgeFlow(edgeNum, nodePressure, edgeConductance, edgeList);
    
    for i = 1:size(FRetained,1)
        % Store the number of retained particles in this time step
        node = FRetained(i,1);
        timeNodeRetained.(sprintf('t%i',t))(node,1) = timeNodeRetained.(sprintf('t%i',t))(node,1)+1;
    end
    
    for i = 1:nodeNum
        % Determine edges connected to the node and outflow
        n = nodeConnectivity(i);
        edgeTo = nodeEdgeTo(i,1:n)';
        outflow = edgeFlow(edgeTo).*nodeEdgeDirection(i,1:n)';
        edgeToFlowing = edgeFlowing(edgeTo,1);
        idxEntrances = outflow<=0;
        
        % Determine number of open entrances
        timeNodeRetained.(sprintf('t%i',t))(i,2) = sum(and(idxEntrances, edgeToFlowing));
        
        % Determine number of open exits
        timeNodeRetained.(sprintf('t%i',t))(i,3) = sum(and(~idxEntrances, edgeToFlowing));
        
        % Determine number of clogged entrances
        timeNodeRetained.(sprintf('t%i',t))(i,4) = sum(and(idxEntrances, ~edgeToFlowing));
        
        % Determine number of clogged exits
        timeNodeRetained.(sprintf('t%i',t))(i,5) = sum(and(~idxEntrances, ~edgeToFlowing));
    end
elseif strcmp(storeOption, 'storeEdgePath')
    % Eroded particles
    idxEroded = find(F(:,5)==-1);
    for i = 1:length(idxEroded)
        idx = idxEroded(i);
        edgePathEroded = edgePathEroded+particleRoute.(sprintf('idx%i',idx));
    end
    
    % Retained particles
    idxRetained = find(F(:,5)>0);
    for i = 1:length(idxRetained)
        idx = idxRetained(i);
        edgePathRetained = edgePathRetained+particleRoute.(sprintf('idx%i',idx));
    end
    
    % Reindex particleRoute
    idxRemaining = find(F(:,5)==0);
    for i = 1:length(idxRemaining)
        idx = idxRemaining(i);
        particleRoute.(sprintf('idx%i',i)) = particleRoute.(sprintf('idx%i',idx));
    end
    
    % Remove saved data
    finesNum = size(F,1);
    for i = length(idxRemaining)+1:finesNum
        particleRoute = rmfield(particleRoute,(sprintf('idx%i',i)));
    end
end

% Remaining active particles
F = F(F(:,5)==0,:);
end
function [timeRetained, timeEroded, timeIntroduced, timeOutside, timeTotalFlow] = finaliseArrays(timeRetained, timeEroded, timeIntroduced, timeOutside, timeTotalFlow, maxFinesTimeStep, t)
timeRetained = timeRetained(1:t-1,:);
timeEroded = timeEroded(1:t-1,:);
timeIntroduced = timeIntroduced(1:maxFinesTimeStep,:);
timeOutside = timeOutside(1:maxFinesTimeStep,:);
timeTotalFlow = timeTotalFlow(1:t-1,:);
end