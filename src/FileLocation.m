classdef FileLocation
    % The class unifies the file locations and names
    
    properties (Constant)
        % Main data folder
        dataFolder = ['..', filesep, 'data', filesep];
        
        % Sub-folders in data folder
        dtmRawFolder = ['dtmRaw', filesep];
        dtmCleanedFolder = ['dtmCleaned', filesep];
        networkFolder = ['network', filesep];
        analysisFolder = ['analysis', filesep];
        flowSolverFolder = ['flowSolver', filesep];
        simulationFolder = ['simulation', filesep];
        
        % Sub-folder(s) in analysis folder
        percolationThresholdFolder = ['percolationThreshold', filesep];
        
        % Sub-folder(s) in simulation folder
        monodisperseFolder = ['monodisperse', filesep];
        polydisperseFolder = ['polydisperse', filesep];
        validationFolder = ['validation', filesep];
        
        % Sub-folder(s) in monodisperse and polydisperse folder
        normalFolder = ['normal', filesep];
        lessFolder = ['less', filesep]; % For simulations with t = 5
        partsFolder = ['parts', filesep]; % Simulations in 5 parts
        storePathFolder = ['storePath', filesep];
        
        % Files in flowSolver folder
        flowCSV = 'data_outputPNM.csv';
        flowTXT = 'output_openPNM.txt';
    end
    
    methods (Static)
        % Filter
        function [s] = filterDTMRaw(filter, component)
            % Unpack parameters
            cuCoarse = filter.filterParameters.cuCoarse;
            density = filter.filterParameters.density;
            size = filter.filterParameters.size;
            triangulationMethod = filter.filterParameters.triangulationMethod;
            
            cuFolder = ['cu', num2str(floor(cuCoarse)), 'pt', num2str(round(cuCoarse*10-floor(cuCoarse)*10)), filesep];
            densityFolder = [density, filesep];
            sizeFolder = [sprintf('%ix%ix%i',size(1),size(2),size(3)), filesep];
            switch component
                case 'node'
                    fileName = ['node', triangulationMethod, '_rev.mat'];
                case 'edge'
                    fileName = ['edge', triangulationMethod, '_rev.mat'];
                case 'particle'
                    fileName = 'particle.mat';
            end
            s = [FileLocation.dataFolder, FileLocation.dtmRawFolder, cuFolder, densityFolder, sizeFolder, fileName];
        end
        function [s] = filterDTMRawVariable(filter, component)
            % Unpack parameters
            triangulationMethod = filter.filterParameters.triangulationMethod;
            
            switch component
                case 'node'
                    s = ['node', triangulationMethod];
                case 'edge'
                    s = ['edge', triangulationMethod];
                case 'particle'
                    s = 'Q';
            end
        end
        function [s] = filterDTMCleaned(filter)
            % Unpack parameters
            cuCoarse = filter.filterParameters.cuCoarse;
            density = filter.filterParameters.density;
            size = filter.filterParameters.size;
            triangulationMethod = filter.filterParameters.triangulationMethod;
            
            cuFolder = ['cu', num2str(floor(cuCoarse)), 'pt', num2str(round(cuCoarse*10-floor(cuCoarse)*10)), filesep];
            densityFolder = [density, filesep];
            sizeFolder = [sprintf('%ix%ix%i',size(1),size(2),size(3)), filesep];
            dtmCleaned = ['filter', triangulationMethod, '.mat'];
            s = [FileLocation.dataFolder, FileLocation.dtmCleanedFolder, cuFolder, densityFolder, sizeFolder, dtmCleaned];
        end
        function [s] = filterCut(filter)
            % Unpack parameters
            cuCoarse = filter.filterParameters.cuCoarse;
            density = filter.filterParameters.density;
            size = filter.filterParameters.size;
            triangulationMethod = filter.filterParameters.triangulationMethod;
            sizeCutFrom = Constants.cutFrom(size);
            
            cuFolder = ['cu', num2str(floor(cuCoarse)), 'pt', num2str(round(cuCoarse*10-floor(cuCoarse)*10)), filesep];
            densityFolder = [density, filesep];
            sizeFolder = [sprintf('%ix%ix%i',sizeCutFrom(1),sizeCutFrom(2),sizeCutFrom(3)), filesep];
            dtmCleaned = ['filter', triangulationMethod, sprintf('%i%i%i',size(1),size(2),size(3)), '.mat'];
            s = [FileLocation.dataFolder, FileLocation.dtmCleanedFolder, cuFolder, densityFolder, sizeFolder, dtmCleaned];
        end
        function [s] = filterBoundarySize(filter)
            % Unpack parameters
            cuCoarse = filter.filterParameters.cuCoarse;
            density = filter.filterParameters.density;
            size = filter.filterParameters.size;
            
            cuFolder = ['cu', num2str(floor(cuCoarse)), 'pt', num2str(round(cuCoarse*10-floor(cuCoarse)*10)), filesep];
            densityFolder = [density, filesep];
            sizeFolder = [sprintf('%ix%ix%i',size(1),size(2),size(3)), filesep];
            fileName = 'log.file';
            s = [FileLocation.dataFolder, FileLocation.dtmRawFolder, cuFolder, densityFolder, sizeFolder, fileName];
        end
        function [s] = filterNetwork(filter)
            % Unpack parameters
            cuCoarse = filter.filterParameters.cuCoarse;
            density = filter.filterParameters.density;
            size = filter.filterParameters.size;
            triangulationMethod = filter.filterParameters.triangulationMethod;
            networkType = filter.networkParameters.networkType;
            flowDirection = filter.networkParameters.flowDirection;
            directedOrUndirected = filter.networkParameters.directedOrUndirected;
            nx = filter.networkParameters.nx;
            ny = filter.networkParameters.ny;
            nz = filter.networkParameters.nz;
            
            networkTypeFolder = [networkType, filesep];
            if strcmp(networkType, 'actual')
                % For actual networks
                cuFolder = ['cu', num2str(floor(cuCoarse)), 'pt', num2str(round(cuCoarse*10-floor(cuCoarse)*10)), filesep];
                densityFolder = [density, filesep];
                network = ['networkactualcu', num2str(floor(cuCoarse)), 'pt', num2str(round(cuCoarse*10-floor(cuCoarse)*10)), density, sprintf('%i%i%i',size(1),size(2),size(3)), triangulationMethod, flowDirection, directedOrUndirected, '.mat'];
                s = [FileLocation.dataFolder, FileLocation.networkFolder, networkTypeFolder, cuFolder, densityFolder, network];
            else
                % For regular networks
                network = ['network', networkType, num2str(nx), num2str(ny), num2str(nz), directedOrUndirected, '.mat'];
                s = [FileLocation.dataFolder, FileLocation.networkFolder, networkTypeFolder, network];
            end
        end
        function [s] = filterShortestPath(filter)
            % Unpack parameters
            cuCoarse = filter.filterParameters.cuCoarse;
            density = filter.filterParameters.density;
            size = filter.filterParameters.size;
            triangulationMethod = filter.filterParameters.triangulationMethod;
            networkType = filter.networkParameters.networkType;
            flowDirection = filter.networkParameters.flowDirection;
            directedOrUndirected = filter.networkParameters.directedOrUndirected;
            
            networkTypeFolder = [networkType, filesep];
            if strcmp(networkType, 'actual')
                % For actual networks
                cuFolder = ['cu', num2str(floor(cuCoarse)), 'pt', num2str(round(cuCoarse*10-floor(cuCoarse)*10)), filesep];
                densityFolder = [density, filesep];
                fileName = ['shortestpathactualcu', num2str(floor(cuCoarse)), 'pt', num2str(round(cuCoarse*10-floor(cuCoarse)*10)), density, sprintf('%i%i%i',size(1),size(2),size(3)), triangulationMethod, flowDirection, directedOrUndirected, '.mat'];
                s = [FileLocation.dataFolder, FileLocation.networkFolder, networkTypeFolder, cuFolder, densityFolder, fileName];
            end
        end
        function [s] = filterAnalyticalCSD(filter)
            % Unpack parameters
            cuCoarse = filter.filterParameters.cuCoarse;
            
            cuFolder = ['cu', num2str(floor(cuCoarse)), 'pt', num2str(round(cuCoarse*10-floor(cuCoarse)*10)), filesep];
            fileName = 'filterAnalyticalCSD.mat';
            s = [FileLocation.dataFolder, FileLocation.dtmCleanedFolder, cuFolder, fileName];
        end
        
        % Flow solver (Only when using Python solver)
        function [s] = flowSolverTXT()
            s = [FileLocation.dataFolder, FileLocation.flowSolverFolder, FileLocation.flowTXT];
        end
        function [s] = flowSolverCSV()
            s = [FileLocation.dataFolder, FileLocation.flowSolverFolder, FileLocation.flowCSV];
        end
        
        % Percolation threshold (Section 3.2.2)
        function [s] = percolationThreshold(networkType, assignmentMethod, cuCoarse, triangulationMethod, nxyz, isPeriodic)
            switch networkType
                case 'actual'
                    switch assignmentMethod
                        case 'actual'
                            fileName = sprintf('results_actual%i%s.mat',cuCoarse,triangulationMethod);
                        case 'random'
                            fileName = sprintf('results_%i%s.mat',cuCoarse,triangulationMethod);
                    end
                case 'cubic'
                    if isPeriodic
                        fileName = sprintf('results_%i%i%i.mat',nxyz(1),nxyz(2),nxyz(3));
                    else
                        fileName = sprintf('results_np%i%i%i.mat',nxyz(1),nxyz(2),nxyz(3));
                    end
            end
            s = [FileLocation.dataFolder, FileLocation.analysisFolder, FileLocation.percolationThresholdFolder, fileName];
        end
        
        % For monodisperse simulations
        function [s] = simulatorMonodisperse(transportLogicType, networkType, assignmentMethod, cuCoarse, density, size, triangulationMethod, varargin)
            if isempty(varargin)
                s = [FileLocation.dataFolder, FileLocation.simulationFolder, FileLocation.monodisperseFolder, FileLocation.normalFolder, ...
                    transportLogicType, networkType, assignmentMethod, num2str(round(cuCoarse*10)), density, ...
                    sprintf('%i%i%i',size(1),size(2),size(3)), triangulationMethod, '.mat'];
            else
                s = [FileLocation.dataFolder, FileLocation.simulationFolder, FileLocation.monodisperseFolder, FileLocation.normalFolder, ...
                    transportLogicType, networkType, assignmentMethod, num2str(round(cuCoarse*10)), density, ...
                    sprintf('%i%i%i',size(1),size(2),size(3)), triangulationMethod, varargin{1}, '.mat'];
            end
        end
        function [s] = simulatorMonodisperseLess(transportLogicType, networkType, assignmentMethod, cuCoarse, density, size, triangulationMethod, varargin)
            if isempty(varargin)
                s = [FileLocation.dataFolder, FileLocation.simulationFolder, FileLocation.monodisperseFolder, FileLocation.lessFolder, ...
                    transportLogicType, networkType, assignmentMethod, num2str(round(cuCoarse*10)), density, ...
                    sprintf('%i%i%i',size(1),size(2),size(3)), triangulationMethod, '.mat'];
            else
                s = [FileLocation.dataFolder, FileLocation.simulationFolder, FileLocation.monodisperseFolder, FileLocation.lessFolder, ...
                    transportLogicType, networkType, assignmentMethod, num2str(round(cuCoarse*10)), density, ...
                    sprintf('%i%i%i',size(1),size(2),size(3)), triangulationMethod, varargin{1}, '.mat'];
            end
        end
        function [s] = simulatorMonodisperseParts(transportLogicType, networkType, assignmentMethod, cuCoarse, density, size, triangulationMethod, d, divide)
            if divide
                endWith = 'divided';
            else
                endWith = '';
            end
            
            s = [FileLocation.dataFolder, FileLocation.simulationFolder, FileLocation.monodisperseFolder, FileLocation.partsFolder, ...
                transportLogicType, networkType, assignmentMethod, num2str(round(cuCoarse*10)), density, ...
                sprintf('%i%i%i',size(1),size(2),size(3)), triangulationMethod, ...
                sprintf('%i%s', round(d*1000), endWith), '.mat'];
        end
        function [s] = simulatorMonodisperseStorePaths(transportLogicType, networkType, assignmentMethod, cuCoarse, density, size, triangulationMethod, d)
            s = [FileLocation.dataFolder, FileLocation.simulationFolder, FileLocation.monodisperseFolder, FileLocation.storePathFolder, ...
                transportLogicType, networkType, assignmentMethod, num2str(round(cuCoarse*10)), density, ...
                sprintf('%i%i%i',size(1),size(2),size(3)), triangulationMethod, ...
                sprintf('%i', round(d*1000)), '.mat'];
        end
        function [s] = simulatorValidationPercentageClogged()
            [s] = [FileLocation.dataFolder, FileLocation.simulationFolder, FileLocation.validationFolder, 'percentageClogged.mat'];
        end
        
        % For polydisperse simulations
        function [s] = simulatorPolydisperse(tranportLogicType, networkType, assignmentMethod, cuCoarse, density, size, triangulationMethod, cuFine, varargin)
            if isempty(varargin)
                s = [FileLocation.dataFolder, FileLocation.simulationFolder, FileLocation.polydisperseFolder, FileLocation.normalFolder, ...
                    transportLogicType, networkType, assignmentMethod, num2str(round(cuCoarse*10)), density, ...
                    sprintf('%i%i%i',size(1),size(2),size(3)), triangulationMethod, sprintf('%i',round(cuFine*10)), '.mat'];
            else
                s = [FileLocation.dataFolder, FileLocation.simulationFolder, FileLocation.polydisperseFolder, FileLocation.normalFolder, ...
                    transportLogicType, networkType, assignmentMethod, num2str(round(cuCoarse*10)), density, ...
                    sprintf('%i%i%i',size(1),size(2),size(3)), triangulationMethod, sprintf('%i',round(cuFine*10)), varargin{1}, '.mat'];
            end
        end
        function [s] = simulatorPolydisperseRepeat(transportLogicType, networkType, assignmentMethod, cuCoarse, density, size, triangulationMethod, cuFine, r, varargin)
            if isempty(varargin)
                s = [FileLocation.dataFolder, FileLocation.simulationFolder, FileLocation.polydisperseFolder, FileLocation.partsFolder, ...
                    transportLogicType, networkType, assignmentMethod, num2str(round(cuCoarse*10)), density, ...
                    sprintf('%i%i%i',size(1),size(2),size(3)), triangulationMethod, sprintf('%i',round(cuFine*10)), sprintf('%i',round(r*10)), '.mat'];
            else
                s = [FileLocation.dataFolder, FileLocation.simulationFolder, FileLocation.polydisperseFolder, FileLocation.partsFolder, ...
                    transportLogicType, networkType, assignmentMethod, num2str(round(cuCoarse*10)), density, ...
                    sprintf('%i%i%i',size(1),size(2),size(3)), triangulationMethod, sprintf('%i',round(cuFine*10)), sprintf('%i',round(r*10)), varargin{1}, '.mat'];
            end
        end
        
        % For analysis results
        function [s] = analysisPolydisperseMechanisms()
            fileName = 'analysisPolydisperseMechanisms.mat';
            s = [FileLocation.dataFolder, FileLocation.analysisFolder, FileLocation.polydisperseFolder, fileName];
        end
    end
end