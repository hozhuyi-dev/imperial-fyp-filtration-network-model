% Tested methods as of 13th June 2019

clear
clc
close all

% Add src to path
path(path, ['..', filesep, 'src'])

%% Methods available for Filter objects

% Creating a Filter object
simulator = Simulator;
simulator = simulator.createFilter(1.2,'density','dense','size',[1,1,1],'triangulationMethod','Del');

% Percolation threshold =====
% Determining percolation threshold using undirected network
p = simulator.filter.percolationThresholdAssigned('minP',0.5,'display','on');

% Using directed network (d is the length of the percolating path)
[p, d] = simulator.percolationThresholdAssignedDirected('minP',0.58,'display','on');

% Determine Dxc* and DxF given x =====
D = simulator.filter.percentileEdgeDiameter(50);
D = simulator.filter.percentileParticleDiameter(50);
D = simulator.filter.percentileParticleDiameterVol(50);

% Shortest path =====
% weight can be 'length' or 'number'
[P, d] = simulator.filter.shortestPath('weight', 'length');
simulator.filter.plotShortestPath('weight','length');

% Analytical CSD =====
% by Silveira et. al. 1975 Locke et. al. 2001 Indraratna et. al. 2007
relativeDensity = 0.5;
particleNum = 3;
[csd] = simulator.filter.analyticalCSDCombined(particleNum, relativeDensity);
figure
semilogx(csd(:,2),csd(:,1))

% Path length given node indices =====
[P, d1] = simulator.filter.shortestPath('weight', 'length');
[d2] = simulator.filter.pathLength(P); % d1=d2

% Clusters of nodes with all connected edges smaller than specified base particle diameter d =====
d = 0.4;
[s] = simulator.filter.nodeSmaller(d); % s gives the sizes of clusters in descending order

% Two-point spatial correlation =====
% If using this, allow the script to update nodeList in private Filter 
% method 'generateNetworkTopologyActual', create the filter again

% Updating nodeList is a one-time calculation, as it will be saved after
% the first run

% Given any nodeProperty
nodeProperty = simulator.filter.nodeSurfaceArea;
distance = 1; % Number of edges apart
[r] = simulator.filter.spatialCorrelation(nodeProperty, distance);

% Specifically for nodeVolume
% nodeList is generated for up to 5 edges apart only
[c] = simulator.filter.spatialCorrelationNodeVolume('distanceUpTo',5);

% Plotting edges =====
% (1) idxEdge: Boolean for which edges to plot
% (2) idxColour: Which colour to plot in (up to 10 different colours)

% Example, plotting first 100 edges
idxEdge = zeros(simulator.filter.edgeNum,1);
idxEdge(1:100) = 1;
% 1:50 has the same colour, 51:100 has different colour, and so on
idxColour = zeros(simulator.filter.edgeNum,1);
idxColour(1:50) = 1;
idxColour(51:100) = 2;
% Plot (plotView represents the direction to plot)
simulator.filter.plotEdges(idxEdge, idxColour, 'plotView', 1:3);

%% Methods available for Base objects

% Creating a Base object
simulator = Simulator;
simulator = simulator.createBase(0.3,'cuFine',1,'run',1);

%% Methods available for Simulator objects

% Example run =====
% Creating a Simulator object (see the src Simulator file for lists of
% options)
simulator = Simulator('transportTypeLogic','flowbiased');

% Creating filter and base
simulator = simulator.createFilter(1.2);
simulator = simulator.createBase(0.2,'run',1);
simulator = simulator.createBase(0.3,'run',2);

% Running simulation with base 1 with clean filter
simulator = simulator.simulateFiltration('run',1,'followingRun',1);

% Running simulation with base 2 with clean filter
simulator = simulator.simulateFiltration('run',2,'followingRun',2);
% OR
% Running simulation with base 2 after simulation run 1
simulator = simulator.simulateFiltration('run',2,'followingRun',1);

% Save results
fileLocation = '../test.m'; % In network-transport folder
simulator.saveResults(fileLocation);

% Load from a saved data =====
simulator = Simulator.loadResults(fileLocation);

% Number of particles introduced =====
% Extracted for every run given in an array
n = simulator.introducedNum;

%% Methods available for Statistics object

% Statistics objects are created automatically after simulations

%% Methods available for Validate class

% Constriction size distribution, compared with Shire (2014)
% Shire (2014) data is in the folder '../data/dataShire/'
Validate.plotCompareShire;

% Plot the comparison of permeabilities in regular lattices with papers
Validate.plotPapersRegular;

% Plot the comparison of permeabilities in random sphere packings with papers
Validate.plotPapersActual;