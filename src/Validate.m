classdef Validate
    % Basic validations for the overall implementation of the network
    % model
    
    methods (Static)
        function [bool, test] = allTests()
            % Initialise test
            test = NaN(5,4);
            
            % Test 1.1 Node connectivity must be 4
            [test(1,1),~] = Validate.nodeConnectivity();
            
            % Test 4.1 Compare with Python using a 10x10x10 network
            test(4,1) = Validate.compareWithPythonRegular();
            
            % Test 4.2 Compare with Python using actual network cu1pt2 Mod
            test(4,2) = Validate.compareWithPythonActual();
            
            % Test 4.3 Create a regular cubic network 40x40x60 with
            % constant conductance of 1, check how pressure varies
            test(4,3) = Validate.flowConstantConductance();
            
            % Test 4.4 Create a regular cubic network 40x40x60 with
            % random geometric properties, check for mass conservation
            test(4,4) = Validate.massConservation();
            
            % Test 5.1 Particle movements in cubic network with constant
            % conductance
            test(5,1) = Validate.constantConductance();
            
            % Output
            bool = test(:);
            bool = bool(~isnan(bool));
            bool = all(bool);
        end
        
        % Section 1 ===== Generation of network
        function [] = generateSmallNetwork(networkType, directedOrUndirected)
            switch networkType
                case 'square'
                    nxyz = [3,1,3];
                case 'slantedsquare'
                    nxyz = [3,1,3];
                case 'cubic'
                    nxyz = [3,3,3];
                case 'slantedcubic'
                    nxyz = [2,2,3];
                case 'bcc'
                    nxyz = [3,3,3];
                case 'fcc'
                    nxyz = [2,2,3];
                case 'squarenp'
                    nxyz = [3,1,3];
                case 'slantedsquarenp'
                    nxyz = [3,1,3];
                case 'cubicnp'
                    nxyz = [3,3,3];
                case 'slantedcubicnp'
                    nxyz = [2,2,3];
                case 'bccnp'
                    nxyz = [3,3,3];
                case 'fccnp'
                    nxyz = [2,2,3];
            end
            
            % Generate network
            filter = Filter([],[],[],[],[],networkType,[],nxyz,[],directedOrUndirected);
            filter = filter.generateNetworkTopologyRegular();
            
            % (1) ===== nodeNum, edgeNum
            % Display the nodeNum and edgeNum
            fprintf('Validation: The numbers of nodes and edges are %i and %i respectively.\n', filter.nodeNum, filter.edgeNum);
            
            % (2) ===== nodeCoord, nodeLattice, edgeList
            % Generate a plot
            figure;
            set(gcf, 'WindowState', 'maximized');
            for i = 1:filter.edgeNum
                plot3(filter.nodeLattice(filter.edgeList(i,:),1), filter.nodeLattice(filter.edgeList(i,:),2), filter.nodeLattice(filter.edgeList(i,:),3), 'Color', 0.3*ones(3,1));
                hold on;
            end
            scatter3(filter.nodeLattice(:,1), filter.nodeLattice(:,2), filter.nodeLattice(:,3))
            hold off;
            set(gca,'FontSize',18);
            set(gca, 'DataAspectRatio', [1,1,1])
            title(['Validation: ', networkType, sprintf(' lattice: %i %i %i, transport logic: ', nxyz(1), nxyz(2), nxyz(3)), directedOrUndirected], 'FontSize', 24);
            
            % (3) ===== nodeConnectivity, nodeInlet, nodeOutlet, nodeEdgeTo, nodeNodeTo
            % Count the number of inlets and outlets
            fprintf('Validation: The numbers of inlets and outlets are %i and %i respectively.\n', sum(filter.nodeInlet), sum(filter.nodeOutlet));
            
            % Connectivity of nodes
            fprintf('Validation: Inlets have on average %i edges, internal nodes have %i, outlets have %i.\n', ...
                sum(filter.nodeConnectivity(filter.nodeInlet==1))/sum(filter.nodeInlet), ...
                sum(filter.nodeConnectivity((filter.nodeInlet+filter.nodeOutlet)==0))/sum((filter.nodeInlet+filter.nodeOutlet)==0), ...
                sum(filter.nodeConnectivity(filter.nodeOutlet==1))/sum(filter.nodeOutlet));
        end
        function [bool, test] = nodeConnectivity()
            % Data available
            cuCoarseList = Constants.cuCoarse;
            densityList = Constants.density;
            sizeList = Constants.size;
            sizeList = sizeList(2:end,:);
            % sizeList = [1,1,4];
            triangulationMethod = 'Del';
            
            % Initialise test
            test = ones(length(cuCoarseList),length(densityList),size(sizeList,1));
            
            for i = 1:length(cuCoarseList)
                if cuCoarseList(i)==3
                    n = size(sizeList,1)-4;
                else
                    n = size(sizeList,1);
                end
                for j = 1:length(densityList)
                    for k = 1:n
                        % Create filter
                        simulator = Simulator();
                        simulator = simulator.createFilter(cuCoarseList(i), 'density', densityList{j}, ...
                            'size', sizeList(k,:), 'triangulationMethod', triangulationMethod);
                        
                        % Determine nodeConnectivity
                        nodeConnectivity = simulator.filter.nodeConnectivity+max(simulator.filter.nodeBoundary,[],2);
                        
                        % All nodeConnectivity must be 4
                        test(i,j,k) = all(nodeConnectivity==4);
                    end
                end
            end
            
            % Bool
            bool = all(all(all(test)));
        end
        
        % Section 2 ===== Assignment of geometric properties
        function [] = plotNetworkProperties(simulator)
            % Unpack
            cuCoarse = simulator.filter.filterParameters.cuCoarse;
            density = simulator.filter.filterParameters.density;
            size = simulator.filter.filterParameters.size;
            triangulationMethod = simulator.filter.filterParameters.triangulationMethod;
            
            % Generate an actual network
            filterActual = Filter(1,cuCoarse,density,size,triangulationMethod,'actual','actual',[],'zpos','undirected');
            filterActual = filterActual.loadDTM('read');
            filterActual = filterActual.generateNetwork();
            
            % Current filter under consideration
            filter = simulator.filter;
            
            % Plot nodeSurfaceArea, nodeVolume, edgeArea, edgeDiameter,
            % edgeLength, edgeReduced
            plotList = {'nodeSurfaceArea', 'nodeVolume', 'edgeArea', 'edgeDiameter', 'edgeLength', 'edgeReduced'};
            for i = 1:length(plotList)
                % Create CDF
                CDF = Statistics.createCDF(filter.(sprintf('%s',plotList{i})),20);
                
                % Plot
                figure;
                set(gcf, 'WindowState', 'maximized');
                h = cdfplot(filterActual.(sprintf('%s',plotList{i})));
                h.LineWidth = 2.5;
                hold on;
                plot(CDF(:,2),CDF(:,1)/100,'LineStyle','none','LineWidth',2.5,'Marker','x');
                hold off;
                set(gca, 'FontSize', 14);
                legend('Expected CDF', 'Generated CDF', 'FontSize', 14);
                xlabel(sprintf('Geometric property: %s', plotList{i}), 'FontSize', 18);
                ylabel('Percentage smaller by number', 'FontSize', 18);
                title(sprintf('Validation: Generated Geometric Property: %s', plotList{i}), 'FontSize', 24);
            end
        end
        function [] = plotCompareShire()
            % Data available
            cuCoarseList = Constants.cuCoarse;
            triangulationMethod = 'Mod';
            densityList = {'dense'};
            size = [0,0,0];
            
            % Initialise struct
            filter = struct;
            
            % Load data
            k = 0;
            for i = 1:length(cuCoarseList)
                for j = 1:length(densityList)
                    % Current iteration
                    k = k+1;
                    
                    % Create filter object
                    simulator = Simulator;
                    simulator = simulator.createFilter(cuCoarseList(i),'density',densityList{j},'size',size,'triangulationMethod',triangulationMethod);
                    filter.(sprintf('filter%i',k)) = simulator.filter;
                end
            end
            
            % Load Shire's data
            csdNum = struct;
            for i = 1:length(cuCoarseList)
                [csdNum.(sprintf('csdNum%i',i)), ~, ~, ~, ~, ~, ~] = Validate.generateInputShire(1.2, 0.5, ...
                    sprintf('Cu%ipt%i',floor(cuCoarseList(i)),round((cuCoarseList(i)-floor(cuCoarseList(i)))*10)), 'dense');
            end
            
            % Plot properties
            colourList = Constants.colourList;
            
            % Legend list
            legendList = cell(length(cuCoarseList)*(1+length(densityList)),1);
            k = 0;
            for i = 1:length(cuCoarseList)
                for j = 1:length(densityList)
                    k = k+1;
                    legendList{k} = sprintf('C_u = %.1f Density = %s', cuCoarseList(i), densityList{j});
                end
                legendList{i+(length(cuCoarseList)*length(densityList))} = sprintf('C_u = %.1f Shire (2014)', cuCoarseList(i));
            end
            
            % Plot
            figure;
            set(gcf,'WindowState','maximized');
            k = 0;
            for i = 1:length(cuCoarseList)
                for j = 1:length(densityList)
                    % Current iteration
                    k = k+1;
                    
                    % Create CDF
                    CDF = Statistics.createCDF(filter.(sprintf('filter%i',k)).edgeDiameter,1e3);
                    
                    % Plot
                    semilogx(CDF(:,2),CDF(:,1),'LineWidth',5,'LineStyle','-.','Color',colourList(i,:));
                    hold on;
                end
            end
            for i = 1:length(cuCoarseList)
                % Create CDF
                CDF = csdNum.(sprintf('csdNum%i',i));
                
                % Plot
                semilogx(CDF(:,2),CDF(:,1)*100,'LineWidth',5,'LineStyle','-','Color',colourList(i,:));
            end
            legend(legendList,'Location','southeast','FontName','Helvetica','FontSize',36);
            set(gca,'FontName','Helvetica','FontSize',36);
            xticks([0.25,0.5,1]);
            xticklabels({'0.25','0.50','1.00'});
            % xlabel('Normalised constriction diameter (D_F/D_0_F)','FontName','Helvetica','FontSize',36);
            ylabel('Percentage smaller','FontName','Helvetica','FontSize',36);
        end
        
        % Section 3 ===== Introduction of fines
        function [] = plotIntroduced(simulator, varargin)
            % Input parser
            parser = inputParser;
            addParameter(parser, 'run', 1);
            parse(parser, varargin{:});
            run = parser.Results.run;
            
            % Expected PSDs
            expectedPSDVol = simulator.base.(sprintf('base%i',run)).psdFineVol;
            expectedPSDNum = simulator.base.(sprintf('base%i',run)).psdFineNum;
            
            % Generated PSDs
            generatedPSDVol = Statistics.createPSDVolBucket(simulator.statistics.(sprintf('statistics%i',run)).timeIntroduced(end,4:end), simulator.base.(sprintf('base%i',run)).psdFineBucket);
            generatedPSDNum = Statistics.createPSDNumBucket(simulator.statistics.(sprintf('statistics%i',run)).timeIntroduced(end,4:end));
            
            % Plot properties
            colourList = Constants.colourList;
            
            % Plot PSD by volume
            figure;
            set(gcf, 'WindowState', 'maximized');
            semilogx(expectedPSDVol(:,2),expectedPSDVol(:,1),'Color',colourList(1,:),'LineWidth',2.5);
            hold on;
            semilogx(simulator.base.(sprintf('base%i',run)).psdFineBucket,generatedPSDVol,'x','Color',colourList(1,:),'LineWidth',2.5);
            semilogx(expectedPSDNum(:,2),expectedPSDNum(:,1), 'Color', colourList(2,:),'LineWidth',2.5);
            semilogx(simulator.base.(sprintf('base%i',run)).psdFineBucket,generatedPSDNum,'x','Color',colourList(2,:),'LineWidth',2.5);
            set(gca, 'FontSize', 14);
            legend('Expected PSD by volume', 'Generated PSD by volume', 'Expected PSD by number', 'Generated PSD by number', 'FontSize', 14);
            xlabel('Fine particle size scaled by D0', 'FontSize', 18);
            ylabel('Frequency', 'FontSize', 18);
            title('Validation: Generated PSDs', 'FontSize', 24);
        end
        
        % Section 4 ===== Flow Solver
        function [bool] = compareWithPythonRegular()
            % Generate a regular network (cubic) 10x10x10
            filter = Filter(1,1.2,'dense',[1,1,1],'Del','cubic',[],[10,10,10],[],'undirected');
            filter = filter.loadDTM('read');
            filter = filter.generateNetwork();
            
            % Flow solver parameters
            pressureEntry = 2;
            pressureExit = 1;
            
            % Unpack parameters
            [nodeNum, nodeConnectivity, nodeCoord,nodeEdgeTo, nodeInlet, nodeNodeTo, nodeOutlet, nodeSurfaceArea, nodeVolume] = Simulator.unpack(filter, 'nodeNum', 'nodeConnectivity', 'nodeCoord', 'nodeEdgeTo', 'nodeInlet', 'nodeNodeTo', 'nodeOutlet', 'nodeSurfaceArea', 'nodeVolume');
            [edgeNum, edgeArea, edgeLength, edgeList, edgeReduced] = Simulator.unpack(filter, 'edgeNum', 'edgeArea', 'edgeLength', 'edgeList', 'edgeReduced');
            [flowDirection, networkType] = Simulator.unpack(filter.networkParameters, 'flowDirection', 'networkType');
            nodeActive = zeros(filter.nodeNum,3);
            edgeFlowing = ones(filter.edgeNum,1);
            nodeRetained = 0;
            
            % Determine local node pressures with Matlab and Python
            [nodePressureMatlab, ~, ~] = Simulator.determinePressure(nodeNum, edgeNum, nodeActive, nodeConnectivity, nodeCoord,nodeEdgeTo, nodeInlet, nodeNodeTo, nodeOutlet, nodeRetained, nodeSurfaceArea, nodeVolume, edgeArea, edgeFlowing, edgeLength, edgeList, edgeReduced, pressureEntry, pressureExit, 'matlab', flowDirection, networkType);
            [nodePressurePython, ~, ~] = Simulator.determinePressure(nodeNum, edgeNum, nodeActive, nodeConnectivity, nodeCoord,nodeEdgeTo, nodeInlet, nodeNodeTo, nodeOutlet, nodeRetained, nodeSurfaceArea, nodeVolume, edgeArea, edgeFlowing, edgeLength, edgeList, edgeReduced, pressureEntry, pressureExit, 'python', flowDirection, networkType);
            
            % Determine the maximum percentage difference
            maxDiff = max(abs(nodePressureMatlab-nodePressurePython)./nodePressureMatlab);
            
            % Output
            bool = maxDiff<1e-4;
        end
        function [bool] = compareWithPythonActual()
            % Generate an actual network cu1pt2 Mod
            filter = Filter(1,1.2,'dense',[1,1,1],'Del','actual','actual',[],'zpos','undirected');
            filter = filter.loadDTM('read');
            filter = filter.generateNetwork();
            
            % Flow solver parameters
            pressureEntry = 2;
            pressureExit = 1;
            
            % Unpack parameters
            [nodeNum, nodeConnectivity, nodeCoord,nodeEdgeTo, nodeInlet, nodeNodeTo, nodeOutlet, nodeSurfaceArea, nodeVolume] = Simulator.unpack(filter, 'nodeNum', 'nodeConnectivity', 'nodeCoord', 'nodeEdgeTo', 'nodeInlet', 'nodeNodeTo', 'nodeOutlet', 'nodeSurfaceArea', 'nodeVolume');
            [edgeNum, edgeArea, edgeLength, edgeList, edgeReduced] = Simulator.unpack(filter, 'edgeNum', 'edgeArea', 'edgeLength', 'edgeList', 'edgeReduced');
            [flowDirection, networkType] = Simulator.unpack(filter.networkParameters, 'flowDirection', 'networkType');
            nodeActive = zeros(filter.nodeNum,3);
            edgeFlowing = ones(filter.edgeNum,1);
            nodeRetained = 0;
            
            % Determine local node pressures with Matlab and Python
            [nodePressureMatlab, ~, ~] = Simulator.determinePressure(nodeNum, edgeNum, nodeActive, nodeConnectivity, nodeCoord,nodeEdgeTo, nodeInlet, nodeNodeTo, nodeOutlet, nodeRetained, nodeSurfaceArea, nodeVolume, edgeArea, edgeFlowing, edgeLength, edgeList, edgeReduced, pressureEntry, pressureExit, 'matlab', flowDirection, networkType);
            [nodePressurePython, ~, ~] = Simulator.determinePressure(nodeNum, edgeNum, nodeActive, nodeConnectivity, nodeCoord,nodeEdgeTo, nodeInlet, nodeNodeTo, nodeOutlet, nodeRetained, nodeSurfaceArea, nodeVolume, edgeArea, edgeFlowing, edgeLength, edgeList, edgeReduced, pressureEntry, pressureExit, 'python', flowDirection, networkType);
            
            % Determine the maximum percentage difference
            maxDiff = max(abs(nodePressureMatlab-nodePressurePython)./nodePressureMatlab);
            fprintf('The maximum percentage difference is %f. \n', maxDiff*100);
            
            % Output
            bool = maxDiff<1e-5;
        end
        function [bool] = flowConstantConductance()
            % Initialise test
            test = zeros(4,1);
            tolerance = 1e-8;
            
            % Generate regular cubic network 40x40x60 with conductance of 1
            dynamicViscosity = Constants.dynamicViscosity;
            shapeFactor = Constants.shapeFactor;
            filter = Filter([],[],[],[],[],'cubic',[],[40,40,60],[],'undirected');
            filter = filter.generateNetworkTopologyRegular();
            filter = filter.assignConstantConductance(1);
            filter.edgeArea = filter.edgeArea*dynamicViscosity/shapeFactor;
            
            % Flow solver parameters
            pressureEntry = 60;
            pressureExit = 1;
            flowMethod = 'matlab';
            
            % Unpack parameters
            [nodeNum, nodeConnectivity, nodeCoord,nodeEdgeTo, nodeInlet, nodeNodeTo, nodeOutlet, nodeSurfaceArea, nodeVolume] = Simulator.unpack(filter, 'nodeNum', 'nodeConnectivity', 'nodeCoord', 'nodeEdgeTo', 'nodeInlet', 'nodeNodeTo', 'nodeOutlet', 'nodeSurfaceArea', 'nodeVolume');
            [edgeNum, edgeArea, edgeLength, edgeList, edgeReduced] = Simulator.unpack(filter, 'edgeNum', 'edgeArea', 'edgeLength', 'edgeList', 'edgeReduced');
            [flowDirection, networkType] = Simulator.unpack(filter.networkParameters, 'flowDirection', 'networkType');
            nodeActive = zeros(filter.nodeNum,3);
            edgeFlowing = ones(filter.edgeNum,1);
            nodeRetained = 0;
            
            % Determine local node pressures with Matlab and Python
            [nodePressure, nodeSource, edgeConductance] = Simulator.determinePressure(nodeNum, edgeNum, nodeActive, nodeConnectivity, nodeCoord,nodeEdgeTo, nodeInlet, nodeNodeTo, nodeOutlet, nodeRetained, nodeSurfaceArea, nodeVolume, edgeArea, edgeFlowing, edgeLength, edgeList, edgeReduced, pressureEntry, pressureExit, flowMethod, flowDirection, networkType);
            
            % Test 1 - Perfect linear variation
            maxDiff = max(abs(nodePressure+filter.nodeLattice(:,3)-61));
            test(1) = maxDiff<tolerance;
            
            % Test 2 - Flow of 1 in edges parallel to direction of pressure gradient
            % Determine the absolute flows in each edge
            edgeParallelFlow = zeros(edgeNum,1);
            for i = 1:edgeNum
                if filter.nodeLattice(edgeList(i,1),3)~=filter.nodeLattice(edgeList(i,2),3)
                    edgeParallelFlow(i) = abs(edgeConductance(i)*abs(nodePressure(edgeList(i,1))-nodePressure(edgeList(i,2)))-1);
                end
            end
            maxDiff = max(edgeParallelFlow);
            test(2) = maxDiff<tolerance;
            
            % Test 3 - Flow of 0 in edges perpendicular to direction of pressure gradient
            % Determine the absolute flows in each edge
            edgePerpendicularFlow = zeros(edgeNum,1);
            for i = 1:edgeNum
                if filter.nodeLattice(edgeList(i,1),3)==filter.nodeLattice(edgeList(i,2),3)
                    edgePerpendicularFlow(i) = edgeConductance(i)*abs(nodePressure(edgeList(i,1))-nodePressure(edgeList(i,2)));
                end
            end
            maxDiff = max(edgePerpendicularFlow);
            test(3) = maxDiff<tolerance;
            
            % Test 4 - Total flow that is equal to 1600
            maxDiff = abs(sum(nodeSource(nodeInlet==1))-1600);
            test(4) = maxDiff<tolerance;
            
            % Output
            bool = all(test);
        end
        function [bool] = massConservation()
            % Initialise array(s)
            test = zeros(2,1);
            tolerance = 1e-8;
            
            % Generate a network with random conductance
            filter = Filter(1,1.2,'dense',[1,1,1],'Del','cubic',[],[40,40,60],[],'undirected');
            filter = filter.loadDTM('read');
            filter = filter.generateNetwork();
            
            % Flow solver parameters
            pressureEntry = 1;
            pressureExit = 0;
            flowMethod = 'matlab';
            
            % Unpack parameters
            [nodeNum, nodeConnectivity, nodeCoord,nodeEdgeTo, nodeInlet, nodeNodeTo, nodeOutlet, nodeSurfaceArea, nodeVolume] = Simulator.unpack(filter, 'nodeNum', 'nodeConnectivity', 'nodeCoord', 'nodeEdgeTo', 'nodeInlet', 'nodeNodeTo', 'nodeOutlet', 'nodeSurfaceArea', 'nodeVolume');
            [edgeNum, edgeArea, edgeLength, edgeList, edgeReduced] = Simulator.unpack(filter, 'edgeNum', 'edgeArea', 'edgeLength', 'edgeList', 'edgeReduced');
            [flowDirection, networkType] = Simulator.unpack(filter.networkParameters, 'flowDirection', 'networkType');
            nodeActive = zeros(filter.nodeNum,3);
            edgeFlowing = ones(filter.edgeNum,1);
            nodeRetained = 0;
            
            % Determine local node pressures with Matlab and Python
            [nodePressure, nodeSource, edgeConductance] = Simulator.determinePressure(nodeNum, edgeNum, nodeActive, nodeConnectivity, nodeCoord,nodeEdgeTo, nodeInlet,nodeNodeTo,nodeOutlet,nodeRetained,nodeSurfaceArea, nodeVolume, edgeArea, edgeFlowing, edgeLength, edgeList, edgeReduced, pressureEntry, pressureExit, flowMethod, flowDirection, networkType);
            
            % Test 1 - Total inflow equals to total outflow
            totalInflow = abs(sum(nodeSource(nodeInlet==1)));
            totalOutflow = abs(sum(nodeSource(nodeOutlet==1)));
            maxDiff = abs(totalInflow-totalOutflow);
            test(1) = maxDiff<tolerance;
            
            % Test 2 - For each node, inflow equals to outflow
            for i = 1:nodeNum
                nodeFrom = i;
                n = nodeConnectivity(i);
                nodeTo = nodeNodeTo(i,1:n)';
                edgeTo = nodeEdgeTo(i,1:n)';
                [outflow] = Simulator.determineOutflow(nodeFrom, nodeTo, edgeTo, edgeConductance, nodePressure);
                totalOutflow(i) = sum(outflow)-nodeSource(i);
            end
            maxDiff = max(abs(totalOutflow));
            test(2) = maxDiff<tolerance;
            
            % Output
            bool = all(test);
        end
        function [] = plotPapersRegular()
            % Compare the permeabilities of fluid in regular simple cubic
            % packings
            
            % Input the concentration
            phiList = linspace(0.000125,0.52359,100);
            
            % Generate regular network
            filter = Filter([],[],[],[],[],'cubic',[],[10,10,10],[],'undirected');
            filter = filter.generateNetworkTopologyRegular();
            filter = filter.assignConstantConductance(1);
            
            % Flow solver parameters
            pressureEntry = 1;
            pressureExit = 0;
            flowMethod = 'matlab';
            
            % Unpack parameters
            [nodeNum, nodeConnectivity, nodeCoord,nodeEdgeTo, nodeInlet, nodeNodeTo, nodeOutlet] = Simulator.unpack(filter, 'nodeNum', 'nodeConnectivity', 'nodeCoord', 'nodeEdgeTo', 'nodeInlet', 'nodeNodeTo', 'nodeOutlet');
            [edgeNum, edgeList, edgeReduced] = Simulator.unpack(filter, 'edgeNum', 'edgeList', 'edgeReduced');
            [flowDirection, networkType, nx, ny, nz] = Simulator.unpack(filter.networkParameters, 'flowDirection', 'networkType', 'nx', 'ny', 'nz');
            nodeActive = zeros(filter.nodeNum,3);
            edgeFlowing = ones(filter.edgeNum,1);
            nodeRetained = 0;
            
            % Determine the permeability with the adopted conductance model
            % and with Kozeny-Carman equation
            k = zeros(size(phiList));
            kkc = zeros(size(phiList));
            eList = zeros(size(phiList));
            for i = 1:length(phiList)
                % Current phi
                phi = phiList(i);
                
                % Determine the diameter of spheres and distance between spheres
                r = 1;
                l = ((4/3)*pi*(r^3)/phi)^(1/3);
                d = r*2;
                assert(l>=d);
                
                % Determine the void ratio
                Vtot = l^3;
                Vs = (4/3)*pi*(r^3);
                Vv = Vtot-Vs;
                eList(i) = Vv/Vs;
                
                % Determine the node volume and surface area
                nodeVolume = ones(nodeNum,1)*(l^3-(4/3)*pi*(r^3));
                nodeSurfaceArea = ones(nodeNum,1)*(4*pi*(r^2));
                
                % Determine the edge length and area
                edgeArea = ones(edgeNum,1)*((l^2)-(pi*(r^2)));
                edgeLength = ones(edgeNum,1)*l;
                
                % Determine the flow through the network
                dynamicViscosity = Constants.dynamicViscosity;
                [~, nodeSource, ~] = Simulator.determinePressure(nodeNum, edgeNum, nodeActive, nodeConnectivity, nodeCoord, nodeEdgeTo, nodeInlet, nodeNodeTo, nodeOutlet, nodeRetained,nodeSurfaceArea, nodeVolume, edgeArea, edgeFlowing, edgeLength, edgeList, edgeReduced, pressureEntry, pressureExit, flowMethod, flowDirection, networkType);
                Q = sum(nodeSource(nodeInlet==1));
                k(i) = (Q*dynamicViscosity*l*(nz-1))/((l*nx)*(l*ny)*(pressureEntry-pressureExit));
                
                % For Kozeny-Carman estimate of permeability
                Sv = (pi*(d^2))/((1/6)*pi*(d^3));
                kkc(i) = ((1-phi)^3)/(5*(phi^2)*(Sv^2));
            end
            
            % Other results
            phizh = [0.000125,0.001,0.008,0.027,0.064,0.125,0.216,0.343,0.45,0.52359];
            Kzh = [1.096,1.212,1.525,2.008,2.81,4.292,7.442,15.4,28.1,42.1];
            Kss1 = 42.6;
            Kss2 = 40.4;
            Km = 42.5;
            % phih = [0.000125,0.001,0.008,0.027,0.064,0.125,0.216];
            % Kh = [1.096,1.212,1.525,2.009,2.83,4.5,12];
            
            % Fit a cubic polynomial to Kzh
            p = polyfit(phizh,Kzh,3);
            KzhContinuous = polyval(p,phiList);
            
            % Convert to permeabilities
            phi = phiList(end);
            kzhContinuous = (2/9)*(r^2)./(KzhContinuous.*phiList);
            % kzh = (2/9)*(r^2)./(Kzh.*phizh);
            kss1 = (2/9)*(r^2)./(Kss1*phi);
            kss2 = (2/9)*(r^2)./(Kss2*phi);
            % kh = (2/9)*(r^2)./(Kh.*phih);
            km = (2/9)*(r^2)./(Km*phi);
            
            % Determine the density where the adopted model intersects with
            % 1.2*kzh
            f = @(x)(interp1(phiList,1.2*kzhContinuous,x)-interp1(phiList,k,x));
            phiIntersect = fzero(f,0.48);
            fprintf('Intersect at %.3f. \n', phiIntersect);
            
            % Plot parameters
            a = 4;
            b = 700;
            
            % Plot
            figure;
            set(gcf, 'WindowState', 'maximized');
            plot(phiList,k,'LineWidth',Constants.lineWidth);
            hold on;
            ylim([0.005,0.045]);
            plot(phiList,kkc,'LineWidth',Constants.lineWidth,'LineStyle','--');
            plot(phiList,kzhContinuous,'LineWidth',Constants.lineWidth,'Color','k');
            plot(phiList,0.8*kzhContinuous,'LineWidth',Constants.lineWidth,'LineStyle',':','Color','k');
            scatter(phiList(end),kss1,'LineWidth',a,'Marker','d','SizeData',b);
            scatter(phiList(end),kss2,'LineWidth',a,'Marker','^','SizeData',b);
            scatter(phiList(end),km,'LineWidth',a,'Marker','x','SizeData',b);
            % plot(czh,kzh,'LineWidth',2.5,'LineStyle',':','Color',colourList(3,:));
            % scatter(czh,kzh,'LineWidth',a,'Marker','x','SizeData',b);
            % scatter(ch,kh,'LineWidth',a,'Marker','s','SizeData',b); ,'Hasimoto (1959)'
            plot(phiList,1.2*kzhContinuous,'LineWidth',Constants.lineWidth,'LineStyle',':','Color','k');
            plot(ones(2,1)*phiList(end),ylim,'Color','k','LineWidth',Constants.lineWidth,'LineStyle','--');
            hold off;
            set(gca, 'FontName','Helvetica','FontSize',Constants.fontSize);
            xlim([0.4,0.56]);
            [~, objh] = legend('Adopted conductance model','Kozeny-Carman estimate', ...
                'Zick & Homsy (1982)','30% above and below Zick & Homsy (1982)', ...
                'Martin et. al. (1951)','Snyder & Stewart (1966)','Sorensen & Stewart (1974)', ...
                'Location','northeast','FontName','Helvetica','FontSize',30);
            objhl = findobj(objh, 'type', 'patch');
            set(objhl, 'Markersize', 20, 'LineWidth', 5);
            xlabel('Packing density \phi','FontName','Helvetica','FontSize',Constants.fontSize);
            ylabel('Permeability k [m^2]','FontName','Helvetica','FontSize',Constants.fontSize);
            text(0.525,0.015,'Densest simple cubic packing','FontName','Helvetica','FontSize',18)
            text(0.525,0.0135,'\phi = 0.5236','FontName','Helvetica','FontSize',18)
        end
        function [] = plotPapersActual()
            % Which set of data
            cuCoarseList = Constants.cuCoarse;
            densityList = Constants.density;
            size = [1,1,1];
            triangulationMethod = 'Del';
            flowDirection = Constants.flowDirection;
            
            % Load the void ratios of the filters
            e = Constants.voidRatio(cuCoarseList,densityList);
            n = e./(1+e);
            phi = 1-n;
            
            % Initialise empty array(s)
            permeability = zeros(length(cuCoarseList),length(densityList),length(flowDirection));
            voidRatio = permeability;
            packingDensity = permeability;
            
            % Permeability with adopted conductance model
            for i = 1:length(cuCoarseList)
                for j = 1:length(densityList)
                    for k = 1:length(flowDirection)
                        % Generate network
                        filter = Filter(1,cuCoarseList(i),densityList{j},size,triangulationMethod,'actual','actual',[],flowDirection{k},'undirected');
                        filter = filter.loadDTM('read');
                        filter = filter.generateNetwork();
                        
                        % Flow solver parameters
                        pressureEntry = 1;
                        pressureExit = 0;
                        flowMethod = 'matlab';
                        
                        % Unpack parameters
                        [nodeNum, nodeConnectivity, nodeCoord,nodeEdgeTo, nodeInlet, nodeNodeTo, nodeOutlet, nodeSurfaceArea, nodeVolume] = Simulator.unpack(filter, 'nodeNum', 'nodeConnectivity', 'nodeCoord', 'nodeEdgeTo', 'nodeInlet', 'nodeNodeTo', 'nodeOutlet', 'nodeSurfaceArea', 'nodeVolume');
                        [edgeNum, edgeArea, edgeLength, edgeList, edgeReduced] = Simulator.unpack(filter, 'edgeNum', 'edgeArea', 'edgeLength', 'edgeList', 'edgeReduced');
                        [networkType] = Simulator.unpack(filter.networkParameters, 'networkType');
                        nodeActive = zeros(filter.nodeNum,3);
                        edgeFlowing = ones(filter.edgeNum,1);
                        nodeRetained = 0;
                        
                        % Determine L and A
                        switch flowDirection{k}(1)
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
                        
                        % Determine the flow through the network
                        dynamicViscosity = Constants.dynamicViscosity;
                        [~, nodeSource, ~] = Simulator.determinePressure(nodeNum, edgeNum, nodeActive, nodeConnectivity, nodeCoord, nodeEdgeTo, nodeInlet, nodeNodeTo, nodeOutlet, nodeRetained, nodeSurfaceArea, nodeVolume, edgeArea, edgeFlowing, edgeLength, edgeList, edgeReduced, pressureEntry, pressureExit, flowMethod, flowDirection{k}, networkType);
                        Q = sum(nodeSource(nodeInlet==1));
                        permeability(i,j,k) = (Q*dynamicViscosity*L)/(A*(pressureEntry-pressureExit));
                        voidRatio(i,j,k) = e(i,j);
                        packingDensity(i,j,k) = phi(i,j);
                    end
                end
            end
            
            % Plot for each cuCoarse
            for i = 1:length(cuCoarseList)
                % Read particle data
                filter = Filter(1,cuCoarseList(i),'dense',[1,1,1],'Del','actual','actual',[],'zpos','undirected');
                filter = filter.loadDTM('read');
                particleRadius = filter.particleRadius;
                particleDiameter = particleRadius*2;
                
                % Create cdf
                [CDF] = Statistics.createCDF(particleDiameter(:,end),1e4);
                
                % Ratio of surface area to volume
                Sv = sum(pi*(CDF(:,2)).^2)/sum((1/6)*pi*(CDF(:,2).^3));
                
                % Permeability with Kozeny-Carman
                xl = [phi(i,1)-0.02,phi(i,3)+0.02];
                phiList = linspace(phi(i,1)-0.02,phi(i,3)+0.02,1000);
                permeabilityKC = ((1-phiList).^3)./(5*(phiList.^2).*(Sv.^2));
                
                % Figure properties
                sizeData = 700;
                
                % Create legend
                legendList = {'Loose filters with the wDTM','Medium filters with the wDTM','Dense filters with the wDTM', ...
                    'Kozeny-Carman estimate', '25% above and below Kozeny-Carman estimate'};%, 'Kozeny-Carman estimate for monodisperse assembly'};
                colourList = Constants.colourList;
                
                figure;
                set(gcf, 'WindowState', 'maximized');
                
                % Permeabilities
                for j = 1:length(densityList)
                    x = packingDensity(i,j,:);
                    x = x(:);
                    y = permeability(i,j,:);
                    y = y(:);
                    scatter(x,y,'LineWidth',Constants.lineWidth,'Marker','x','SizeData',sizeData,'MarkerEdgeColor',colourList(j,:),'MarkerFaceColor',colourList(j,:));
                    hold on;
                end
                plot(phiList,permeabilityKC,'LineWidth',Constants.lineWidth,'Color','k');
                plot(phiList,0.75*permeabilityKC,'LineWidth',Constants.lineWidth,'Color','k','LineStyle',':');
                plot(phiList,1.25*permeabilityKC,'LineWidth',Constants.lineWidth,'Color','k','LineStyle',':');
                xlim(xl);
                set(gca, 'FontName','Helvetica','FontSize',Constants.fontSize);
                if i==2
                    [~, objh] = legend(legendList, 'Location','northeast','FontName','Helvetica','FontSize',28);
                    objhl = findobj(objh, 'type', 'patch');
                    set(objhl, 'Markersize', 20, 'LineWidth', 5);
                end
                xlabel('Packing density \phi','FontName','Helvetica','FontSize',Constants.fontSize);
                ylabel('Permeability k [m^2]','FontName','Helvetica','FontSize',Constants.fontSize);
            end
        end
        function [] = plotFlowVsPressure()
            % Generate an actual network
            filter = Filter(1,1.2,'dense',[1,1,1],'Del','actual','actual',[],'zpos','undirected');
            filter = filter.loadDTM('read');
            filter = filter.generateNetwork();
            
            % Flow solver parameters
            pressureEntry = 1:10;
            pressureExit = 0;
            flowMethod = 'matlab';
            
            % Unpack parameters
            [nodeNum, nodeConnectivity, nodeCoord,nodeEdgeTo, nodeInlet, nodeNodeTo, nodeOutlet, nodeSurfaceArea, nodeVolume] = Simulator.unpack(filter, 'nodeNum', 'nodeConnectivity', 'nodeCoord', 'nodeEdgeTo', 'nodeInlet', 'nodeNodeTo', 'nodeOutlet', 'nodeSurfaceArea', 'nodeVolume');
            [edgeNum, edgeArea, edgeLength, edgeList, edgeReduced] = Simulator.unpack(filter, 'edgeNum', 'edgeArea', 'edgeLength', 'edgeList', 'edgeReduced');
            [flowDirection, networkType] = Simulator.unpack(filter.networkParameters, 'flowDirection', 'networkType');
            nodeActive = zeros(filter.nodeNum,3);
            edgeFlowing = ones(filter.edgeNum,1);
            nodeRetained = 0;
            
            % Initialise array(s)
            totalFlow = zeros(length(pressureEntry),1);
            
            for i = 1:length(pressureEntry)
                % Determine local node pressures with Matlab and Python
                [~, nodeSource, ~] = Simulator.determinePressure(nodeNum, edgeNum, nodeActive, nodeConnectivity, nodeCoord,nodeEdgeTo, nodeInlet, nodeNodeTo, nodeOutlet, nodeRetained,nodeSurfaceArea, nodeVolume, edgeArea, edgeFlowing, edgeLength, edgeList, edgeReduced, pressureEntry(i), pressureExit, flowMethod, flowDirection, networkType);
                totalFlow(i) = sum(nodeSource(nodeInlet==1));
            end
            
            % Plot
            figure;
            set(gcf, 'WindowState', 'maximized');
            plot(pressureEntry-pressureExit,totalFlow,'LineWidth',2);
            hold on;
            plot(pressureEntry-pressureExit,totalFlow,'LineWidth',0.1,'LineStyle','none','Marker','o','MarkerFaceColor','k','Color','k');
            set(gca, 'FontName','Helvetica','FontSize',16);
            xlabel('Pressure difference at inlets and outlets','FontName','Helvetica','FontSize',16);
            ylabel('Total flow through the network','FontName','Helvetica','FontSize',16);
        end
        
        % Section 5 ===== Transport Logic
        function [bool] = constantConductance()
            % The function validates that given a cubic network with
            % constant conductance in each edge, particles should move
            % directly from inlets to outlets in straight lines
            
            % Initialise test
            test = zeros(3,1);
            
            % Perform simulation with d50 of 0.00001 and edgeDiameter of 1
            simulator = Simulator('maxFinesTimeStep',1);
            simulator = simulator.createBase(0.00001,'cuFine',1);
            simulator.filter = Filter([],[],[],[],[],'cubic',[],[10,10,10],[],'undirected');
            simulator.filter = simulator.filter.generateNetworkTopologyRegular();
            simulator.filter = simulator.filter.assignConstantConductance(1);
            simulator = simulator.simulateFiltration();
            
            % Determine the horizontal and vertical edges
            edgeList = simulator.filter.edgeList;
            nodeCoord = simulator.filter.nodeCoord;
            zCoord = [nodeCoord(edgeList(:,1),3),nodeCoord(edgeList(:,2),3)];
            idxVerticalEdge = zCoord(:,1)~=zCoord(:,2);
            idxHorizontalEdge = 1-idxVerticalEdge;
            
            % Test 1 - No particle flows horizontally
            edgePath = simulator.statistics.statistics1.edgePath;
            test(1) = all(edgePath(idxHorizontalEdge==1)<1e-10);
            
            % Test 2 - All particles flow vertically
            test(2) = all(edgePath(idxVerticalEdge==1)>0.9);
            
            % Test 3 - Ensure continuity
            edgeListVertical = edgeList(idxVerticalEdge==1,:);
            edgePathVertical = edgePath(idxVerticalEdge==1);
            nodeInlet = simulator.filter.nodeInlet;
            nodeOutlet = simulator.filter.nodeOutlet;
            idxInlet = find(nodeInlet==1);
            idxOutlet = find(nodeOutlet==1);
            inletNum = sum(nodeInlet);
            isNotOutlet = 1;
            test(3) = 1;
            for i = 1:inletNum
                inlet = idxInlet(i);
                idxEdgeTo = find(edgeListVertical(:,1)==inlet,1);
                idxNodeTo = edgeListVertical(idxEdgeTo,2);
                edgePathInlet = edgePathVertical(idxEdgeTo);
                while isNotOutlet
                    % Find the next path
                    idxEdgeTo = find(edgeListVertical(:,1)==idxNodeTo,1);
                    idxNodeTo = edgeListVertical(idxEdgeTo,2);
                    edgePathCurrent = edgePathVertical(idxEdgeTo);
                    
                    % Ensure continuity
                    if edgePathCurrent~=edgePathInlet
                        test(3) = 0;
                        break;
                    end
                    
                    % Is this the outlet
                    if any(idxOutlet==idxNodeTo)
                        isNotOutlet = 0;
                    end
                end
                if test(3)==0
                    break;
                end
            end
            
            % Output
            bool = all(test);
        end
        function [bool] = cloggingMechanism(simulator, varargin)
            % Input parser
            parser = inputParser;
            addParameter(parser, 'run', 1);
            parse(parser, varargin{:});
            run = parser.Results.run;
            
            % Initialise array(s)
            test = zeros(3,1);
            
            % Unpack
            statistics = simulator.statistics.(sprintf('statistics%i',run));
            
            % Test 1 - Mechanism 1: Particle must be larger
            test(1) = all(statistics.edgeFlowing(statistics.edgeFlowing(:,2)==1,4)>simulator.filter.edgeDiameter(statistics.edgeFlowing(:,2)==1));
            
            % Test 2 - Mechanism 3: Particle must be smaller
            test(2) = all(statistics.edgeFlowing(statistics.edgeFlowing(:,2)==3,4)<simulator.filter.edgeDiameter(statistics.edgeFlowing(:,2)==3));
            
            % Test 3 - Mechanism 3: Opposite node must be fully filled
            idxMech3 = find(statistics.edgeFlowing(:,2)==3);
            nodeTo = zeros(size(idxMech3));
            for i = 1:length(idxMech3)
                j = idxMech3(i);
                nodeFrom = statistics.edgeFlowing(j,3);
                nodeTo(i) = sum(simulator.filter.edgeList(j,:)-nodeFrom)+nodeFrom;
            end
            test(3) = all((pi/6)*(statistics.edgeFlowing(idxMech3,4).^3)>simulator.filter.nodeVolume(nodeTo)*simulator.simulatorParameters.densityLimit-statistics.nodeRetained(nodeTo,3));
            
            % Output
            bool = all(test);
        end
        function [bool] = statistics(simulator, varargin)
            % Input parser
            parser = inputParser;
            addParameter(parser, 'run', 1);
            parse(parser, varargin{:});
            run = parser.Results.run;
            
            % Initialise array(s)
            test = zeros(3,1);
            tolerance = 1e-5;
            
            % Unpack
            statistics = simulator.statistics.(sprintf('statistics%i',run));
            
            % Test 1 - nodeRetained against timeRetained
            nodeRetainedNum = sum(statistics.nodeRetained(:,1));
            nodeRetainedVol = sum(statistics.nodeRetained(:,3));
            timeRetainedNum = statistics.timeRetained(end,1);
            timeRetainedVol = statistics.timeRetained(end,2);
            test(1) = nodeRetainedNum==timeRetainedNum;
            maxDiff = abs(nodeRetainedVol-timeRetainedVol);
            test(1) = and(test(1),maxDiff<tolerance);
            
            % Test 2 - timeRetained against layerRetained
            layerRetainedPSD = sum(statistics.layerRetained);
            timeRetainedPSD = statistics.timeRetained(end,4:end);
            test(2) = all(layerRetainedPSD==timeRetainedPSD);
            
            % Test 3 - timeIntroduced against (timeEroded+timeRetained)
            timeIntroducedNum = statistics.timeIntroduced(end,1);
            timeIntroducedVol = statistics.timeIntroduced(end,2);
            timeIntroducedPSD = statistics.timeIntroduced(end,4:end);
            timeErodedNum = statistics.timeEroded(end,1);
            timeErodedVol = statistics.timeEroded(end,2);
            timeErodedPSD = statistics.timeEroded(end,4:end);
            test(3) = all(timeIntroducedPSD==(timeErodedPSD+timeRetainedPSD));
            test(3) = and(test(3),timeIntroducedNum==(timeErodedNum+timeRetainedNum));
            maxDiff = abs(timeIntroducedVol-timeErodedVol-timeRetainedVol);
            test(3) = and(test(3),maxDiff<tolerance);
            
            % Output
            bool = all(test);
        end
        
        % Data from Shire (2014)
        function [dist_CSD_num, dist_PSD_coarse, dist_PSD_fine, dist_PSD_fine_num, dist_PSD_fine_bucket, D50_coarse, result_Terzaghi_coeff] = generateInputShire(fine_cu, d50, coarse_cu, coarse_density)
            % The function generates required network geometries input using functions
            % originally created by Tom Shire
            
            percent_finer = 5;
            coarse_min = 1;             % Keep as 1, vary fine size
            
            [fine_min, fine_max] = determineFinesParameters(fine_cu, d50);
            [dist_CSD_num, dist_PSD_coarse, dist_PSD_fine, dist_PSD_fine_num, D50_coarse, ~, result_Terzaghi_coeff, dist_PSD_fine_bucket] = generateInput(fine_min, fine_max, percent_finer, coarse_min, coarse_cu, coarse_density);
        end
    end
end
% ===== PRIVATE FUNCTIONS =================================================
% Script from Shire (2014), to load his data, specifically CSD
function [fine_min, fine_max] = determineFinesParameters(fine_cu, d50)
% The function determines fine_min and fine_max based on Cu,base and d50/D0
m = 50/log(fine_cu);
if m==Inf
    fine_min = d50;
    fine_max = d50;
else
    c = 50-(log(d50)*m);
    fine_min = exp((0-c)/m);
    fine_max = exp((100-c)/m);
end
assert(fine_min<=d50);
assert(fine_max>=d50);
end
function [dist_CSD_num, dist_PSD_coarse, dist_PSD_fine, dist_PSD_fine_num, D50_coarse, fines_per_void, result_Terzaghi_coeff, dist_PSD_fine_bucket] = generateInput(fine_min, fine_max, percent_finer, coarse_min, coarse_cu, coarse_density)
% The function generates the CSD and PSD by number, as well as D50F
% Original author: Tom Shire

% coarse_min = 1; % best to leave this as 1 and just change fine radii
% coarse_cu = 'Cu4pt5'; % options: cu 1.5, 2.0, 3.0, 4.5, 6.0
% coarse_density = 'loose'; % options: 'loose', 'medium', 'dense', 'very_loose'

% (1) ==== Extract data for CSD_num and PSD_coarse
[dist_CSD_num, dist_PSD_coarse, part_per_void] = PSD_CSD_Coarse_Maker(coarse_cu, coarse_density, coarse_min);

% (2) ==== PSD_fine by VOLUME
% Calculate volumes
vol0 = (fine_min.^3) * pi * (4/3);
volfinemax= (fine_max.^3) * pi * (4/3);

% Calculate number of particles
vol_part_fine = logspace(log10(vol0),log10(volfinemax),10001);
dist = vol_part_fine ./ (pi * (4/3));
dist = dist.^(1/3);
for k = 1 : 10000
    num1 = k;
    num2 = k + 1;
    dist_PSD_fine(num1:num2) = logspace(log10(dist(num1)),log10(dist(num2)),(num2-num1+1));
end
clear dist vol* num*
dist_PSD_fine = dist_PSD_fine';
dist_PSD_fine(:,2) = 0:0.01:100;

% (3) ==== PSD_fine_eroded 'bucket'
erod_ind = 1:100:length(dist_PSD_fine);
dist_PSD_fine_bucket = dist_PSD_fine(erod_ind,1);

% (4) ==== Determine fines_per_void
% Create 1 percent finer
PSD_fine_in_full = dist_PSD_fine;
PSD_fine_in_full(:,2) =  PSD_fine_in_full(:,2) ./ (100/percent_finer);

% Create full PSD by vol
PSD_full = [PSD_fine_in_full; dist_PSD_coarse];

% Get PSD by number - then number of fines particles for each coarse one
% And fine particles for each node in the system!

Vol_part = (PSD_full(:,1)).^3 * pi * (4/3) ;
MeanVol_and_PC_passing = zeros(size(Vol_part,1)-1 , 3);
for m = 1:size(MeanVol_and_PC_passing,1)
    MeanVol_and_PC_passing(m,1) = 0.5 * (Vol_part(m) + Vol_part(m+1));
    MeanVol_and_PC_passing(m,2) = PSD_full(m+1,2) - PSD_full(m,2);
end
max_vol = MeanVol_and_PC_passing(length(MeanVol_and_PC_passing),1);
max_vol_PP = MeanVol_and_PC_passing(length(MeanVol_and_PC_passing),2);
MeanVol_and_PC_passing(:,3) = 100 * (max_vol./MeanVol_and_PC_passing(:,1)) .* (MeanVol_and_PC_passing(:,2)./max_vol_PP);
fines_total = sum(MeanVol_and_PC_passing(1:length(dist_PSD_fine)-1 , 3));
coarse_total = sum(MeanVol_and_PC_passing( length(dist_PSD_fine) : length(MeanVol_and_PC_passing) , 3));
fines_per_coarse = fines_total/coarse_total;

% Take 50% of fines as mobile
fines_per_void =  fines_per_coarse * part_per_void;

% TOTFINE = fines_per_void * 96000

% (5) ==== Create PSD of fines by number
% convert PSD of fines by vol to number
PSD_num_temp = zeros(length(dist_PSD_fine)-1,2);
for n = 1:(length(dist_PSD_fine)-1)
    PSD_num_temp(n,1) = mean(dist_PSD_fine(n:n+1,1));
    PSD_num_temp(n,2) = dist_PSD_fine(n+1,2)-dist_PSD_fine(n,2);
end
Vol_temp = PSD_num_temp(:,1).^3 * pi * (4/3);
Num_temp = (PSD_num_temp(:,2)./Vol_temp) / (PSD_num_temp(1,2)./Vol_temp(1));

% Num_temp = Vol_temp ./ min(Vol_temp);
PP_num = 100 * cumsum(Num_temp) ./ sum(Num_temp);
PP_num = [0; PP_num];
PP_num(end) = 100;
dist_PSD_fine_num = zeros(10001,2);
dist_PSD_fine_num(:,2) = 0:0.01:100;
for k = 1:10001
    ind = find(PP_num >= (k-1)/100, 1, 'first');
    dist_PSD_fine_num(k,1) = dist_PSD_fine(ind,1);
end

% (6) ==== Terzaghi calculation - Get d85 and D15 and calculate coefficient
d85_fine_ind = find(dist_PSD_fine(:,2) <= 85, 1, 'last');
d85_fine = dist_PSD_fine(d85_fine_ind,1);
D15_coarse_ind = find(dist_PSD_coarse(:,2) >= 15, 1, 'first');
D15_coarse = dist_PSD_coarse(D15_coarse_ind,1);
result_Terzaghi_coeff = D15_coarse / d85_fine;

% (7) ==== Determine D50
D50_coarse_ind = find(dist_PSD_coarse(:,2) >= 50, 1, 'first');
D50_coarse = dist_PSD_coarse(D50_coarse_ind,1);
end
function [CSD_num, PSD_coarse, part_per_void] = PSD_CSD_Coarse_Maker(coarse_cu, coarse_density, coarse_min)
% Loads output from delaunay voids analysis and creates input, used in
% generateInput function
% Original author: Tom Shire

coarse_density = coarse_density(1:5);

switch coarse_cu
    case 'Cu1pt2'
        switch coarse_density
            case 'dense'
                PSD_coarse = load(['..', filesep, 'data', filesep, 'dataShire', filesep, 'EffLVL2particles-Cu1pt2_Dense_ReboulJ_LEVEL2_ed']);
                CSD = load(['..', filesep, 'data', filesep, 'dataShire', filesep, 'EffLVL2constrictions-Cu1pt2_Dense_ReboulJ_LEVEL2_ed']);
                Voids = load(['..', filesep, 'data', filesep, 'dataShire', filesep, 'Cu1pt2_Dense_ReboulJ_LEVEL2_edgesOK']);
        end
    case 'Cu4pt5'
        switch coarse_density
            case 'loose'
                PSD_coarse = load(['..', filesep, 'data', filesep, 'dataShire', filesep, 'EffLVL2particles-Cu4pt5_Loose_ReboulJ_LEVEL2_ed']);
                CSD = load(['..', filesep, 'data', filesep, 'dataShire', filesep, 'EffLVL2constrictions-Cu4pt5_Loose_ReboulJ_LEVEL2_ed']);
                Voids = load(['..', filesep, 'data', filesep, 'dataShire', filesep, 'Cu4pt5_Loose_ReboulJ_LEVEL2_edgesOK']);
            case 'mediu'
                PSD_coarse = load(['..', filesep, 'data', filesep, 'dataShire', filesep, 'EffLVL2particles-Cu4pt5_Medium_ReboulJ_LEVEL2_ed']);
                CSD = load(['..', filesep, 'data', filesep, 'dataShire', filesep, 'EffLVL2constrictions-Cu4pt5_Medium_ReboulJ_LEVEL2_ed']);
                Voids = load(['..', filesep, 'data', filesep, 'dataShire', filesep, 'Cu4pt5_Medium_ReboulJ_LEVEL2_edgesOK']);
            case 'dense'
                PSD_coarse = load(['..', filesep, 'data', filesep, 'dataShire', filesep, 'EffLVL2particles-Cu4pt5_Dense_ReboulJ_LEVEL2_ed']);
                CSD = load(['..', filesep, 'data', filesep, 'dataShire', filesep, 'EffLVL2constrictions-Cu4pt5_Dense_ReboulJ_LEVEL2_ed']);
                Voids = load(['..', filesep, 'data', filesep, 'dataShire', filesep, 'Cu4pt5_Dense_ReboulJ_LEVEL2_edgesOK']);
            case 'very_'
                % PSD_coarse = load(['..', filesep, 'data', filesep, 'dataShire', filesep, 'EffLVL2particles-Cu4pt5_Loose_ReboulJ_LEVEL2_ed']);
                % CSD = load(['..', filesep, 'data', filesep, 'dataShire', filesep, 'EffLVL2constrictions-Cu4pt5_Loose_ReboulJ_LEVEL2_ed']);
        end
    case 'Cu1pt5'
        switch coarse_density
            case 'loose'
                PSD_coarse = load(['..', filesep, 'data', filesep, 'dataShire', filesep, 'EffLVL2particles-Cu1pt5_Loose_ReboulJ_LEVEL2_ed']);
                CSD = load(['..', filesep, 'data', filesep, 'dataShire', filesep, 'EffLVL2constrictions-Cu1pt5_Loose_ReboulJ_LEVEL2_ed']);
                Voids = load(['..', filesep, 'data', filesep, 'dataShire', filesep, 'Cu1pt5_Loose_ReboulJ_LEVEL2_edgesOK']);
            case 'mediu'
                PSD_coarse = load(['..', filesep, 'data', filesep, 'dataShire', filesep, 'EffLVL2particles-Cu1pt5_Medium_ReboulJ_LEVEL2_ed']);
                CSD = load(['..', filesep, 'data', filesep, 'dataShire', filesep, 'EffLVL2constrictions-Cu1pt5_Medium_ReboulJ_LEVEL2_ed']);
                Voids = load(['..', filesep, 'data', filesep, 'dataShire', filesep, 'Cu1pt5_Medium_ReboulJ_LEVEL2_edgesOK']);
            case 'dense'
                PSD_coarse = load(['..', filesep, 'data', filesep, 'dataShire', filesep, 'EffLVL2particles-Cu1pt5_Dense_ReboulJ_LEVEL2_ed']);
                CSD = load(['..', filesep, 'data', filesep, 'dataShire', filesep, 'EffLVL2constrictions-Cu1pt5_Dense_ReboulJ_LEVEL2_ed']);
                Voids = load(['..', filesep, 'data', filesep, 'dataShire', filesep, 'Cu1pt5_Dense_ReboulJ_LEVEL2_edgesOK']);
            case 'very_'
                % PSD_coarse = load(['..', filesep, 'data', filesep, 'dataShire', filesep, 'EffLVL2particles-Cu4pt5_Loose_ReboulJ_LEVEL2_ed']);
                % CSD = load(['..', filesep, 'data', filesep, 'dataShire', filesep, 'EffLVL2constrictions-Cu4pt5_Loose_ReboulJ_LEVEL2_ed']);
        end
    case 'Cu2pt0'
        switch coarse_density
            case 'loose'
                PSD_coarse = load(['..', filesep, 'data', filesep, 'dataShire', filesep, 'EffLVL2particles-Cu2_Loose_ReboulJ_LEVEL2_ed']);
                CSD = load(['..', filesep, 'data', filesep, 'dataShire', filesep, 'EffLVL2constrictions-Cu2_Loose_ReboulJ_LEVEL2_ed']);
                Voids = load(['..', filesep, 'data', filesep, 'dataShire', filesep, 'Cu2_Loose_ReboulJ_LEVEL2_edgesOK']);
            case 'mediu'
                PSD_coarse = load(['..', filesep, 'data', filesep, 'dataShire', filesep, 'EffLVL2particles-Cu2_Medium_ReboulJ_LEVEL2_ed']);
                CSD = load(['..', filesep, 'data', filesep, 'dataShire', filesep, 'EffLVL2constrictions-Cu2_Medium_ReboulJ_LEVEL2_ed']);
                Voids = load(['..', filesep, 'data', filesep, 'dataShire', filesep, 'Cu2_Medium_ReboulJ_LEVEL2_edgesOK']);
            case 'dense'
                PSD_coarse = load(['..', filesep, 'data', filesep, 'dataShire', filesep, 'EffLVL2particles-Cu2_Dense_ReboulJ_LEVEL2_ed']);
                CSD = load(['..', filesep, 'data', filesep, 'dataShire', filesep, 'EffLVL2constrictions-Cu2_Dense_ReboulJ_LEVEL2_ed']);
                Voids = load(['..', filesep, 'data', filesep, 'dataShire', filesep, 'Cu2_Dense_ReboulJ_LEVEL2_edgesOK']);
            case 'very_'
                % PSD_coarse = load(['..', filesep, 'data', filesep, 'dataShire', filesep, 'EffLVL2particles-Cu4pt5_Loose_ReboulJ_LEVEL2_ed']);
                % CSD = load(['..', filesep, 'data', filesep, 'dataShire', filesep, 'EffLVL2constrictions-Cu4pt5_Loose_ReboulJ_LEVEL2_ed']);
        end
    case 'Cu3pt0'
        switch coarse_density
            case 'loose'
                PSD_coarse = load(['..', filesep, 'data', filesep, 'dataShire', filesep, 'EffLVL2particles-Cu3_Loose_ReboulJ_LEVEL2_ed']);
                CSD = load(['..', filesep, 'data', filesep, 'dataShire', filesep, 'EffLVL2constrictions-Cu3_Loose_ReboulJ_LEVEL2_ed']);
                Voids = load(['..', filesep, 'data', filesep, 'dataShire', filesep, 'Cu3_Loose_ReboulJ_LEVEL2_edgesOK']);
            case 'mediu'
                PSD_coarse = load(['..', filesep, 'data', filesep, 'dataShire', filesep, 'EffLVL2particles-Cu3_Medium_ReboulJ_LEVEL2_ed']);
                CSD = load(['..', filesep, 'data', filesep, 'dataShire', filesep, 'EffLVL2constrictions-Cu3_Medium_ReboulJ_LEVEL2_ed']);
                Voids = load(['..', filesep, 'data', filesep, 'dataShire', filesep, 'Cu3_Medium_ReboulJ_LEVEL2_edgesOK']);
            case 'dense'
                PSD_coarse = load(['..', filesep, 'data', filesep, 'dataShire', filesep, 'EffLVL2particles-Cu3_Dense_ReboulJ_LEVEL2_ed']);
                CSD = load(['..', filesep, 'data', filesep, 'dataShire', filesep, 'EffLVL2constrictions-Cu3_Dense_ReboulJ_LEVEL2_ed']);
                Voids = load(['..', filesep, 'data', filesep, 'dataShire', filesep, 'Cu3_Dense_ReboulJ_LEVEL2_edgesOK']);
            case 'very_'
                % PSD_coarse = load(['..', filesep, 'data', filesep, 'dataShire', filesep, 'EffLVL2particles-Cu4pt5_Loose_ReboulJ_LEVEL2_ed']);
                % CSD = load(['..', filesep, 'data', filesep, 'dataShire', filesep, 'EffLVL2constrictions-Cu4pt5_Loose_ReboulJ_LEVEL2_ed']);
        end
    case 'Cu6pt0'
        switch coarse_density
            case 'loose'
                PSD_coarse = load(['..', filesep, 'data', filesep, 'dataShire', filesep, 'EffLVL2particles-Cu6_ALL_Loose_ReboulJ_LEVEL2_ed']);
                CSD = load(['..', filesep, 'data', filesep, 'dataShire', filesep, 'EffLVL2constrictions-Cu6_ALL_Loose_ReboulJ_LEVEL2_ed']);
                Voids = load(['..', filesep, 'data', filesep, 'dataShire', filesep, 'Cu6_Loose_ReboulJ_LEVEL2_edgesOK']);
            case 'mediu'
                PSD_coarse = load(['..', filesep, 'data', filesep, 'dataShire', filesep, 'EffLVL2particles-Cu6_Medium_ReboulJ_LEVEL2_ed']);
                CSD = load(['..', filesep, 'data', filesep, 'dataShire', filesep, 'EffLVL2constrictions-Cu6_Medium_ReboulJ_LEVEL2_ed']);
                Voids = load(['..', filesep, 'data', filesep, 'dataShire', filesep, 'Cu6_Medium_ReboulJ_LEVEL2_edgesOK']);
            case 'dense'
                PSD_coarse = load(['..', filesep, 'data', filesep, 'dataShire', filesep, 'EffLVL2particles-Cu6_Dense_ReboulJ_LEVEL2_ed']);
                CSD = load(['..', filesep, 'data', filesep, 'dataShire', filesep, 'EffLVL2constrictions-Cu6_Dense_ReboulJ_LEVEL2_ed']);
                Voids = load(['..', filesep, 'data', filesep, 'dataShire', filesep, 'Cu6_Dense_ReboulJ_LEVEL2_edgesOK']);
            case 'very_'
                % PSD_coarse = load(['..', filesep, 'data', filesep, 'dataShire', filesep, 'EffLVL2particles-Cu4pt5_Loose_ReboulJ_LEVEL2_ed']);
                % CSD = load(['..', filesep, 'data', filesep, 'dataShire', filesep, 'EffLVL2constrictions-Cu4pt5_Loose_ReboulJ_LEVEL2_ed']);
        end
end

Voids = Voids.merge;
Voids = max(Voids);
PSD_coarse = PSD_coarse.Q;
PSD_coarse = PSD_coarse(:,5);
PSD_coarse = sortrows(PSD_coarse);
PSD_scale = coarse_min./PSD_coarse(1);
PSD_coarse = PSD_coarse * PSD_scale;
vol_coarse = PSD_coarse.^3 * pi * (4/3);
PSD_coarse(:,2) = 100 * cumsum(vol_coarse)./sum(vol_coarse);

% Create full CSD
CSD = CSD.SOL2;
CSD_full_num(:,2)  = CSD(:,4);
CSD_full_num(:,2)  = CSD_full_num(:,2) * PSD_scale;
CSD_full_num = sortrows(CSD_full_num,2);
clear CSD;
CSD_full_num(:,1) = (1:length(CSD_full_num)) ./ length(CSD_full_num);
CSD_num = CSD_full_to_1000(CSD_full_num);

% Particles per bond (used to calcaulte how many fines are needed
part_per_void = size(PSD_coarse,1) / Voids;

end
function [CSD_num] = CSD_full_to_1000(CSD_num_full)
% Converts a full CSD to one with 1000 increments of 0.1% each
% Original author: Tom Shire

CSD_num = zeros(1000,2);
CSD_num(:,1) = 0.001:0.001:1;
for k= 1:size(CSD_num,1)
    ind = find(CSD_num_full(:,1) >= (k/1000), 1, 'first');
    CSD_num(k,2) = CSD_num_full(ind,2);
end
end