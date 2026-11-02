% Script to plot graphs to show spatial correlation

clear;
clc;
close all;

% Add src/ to path
path(path, ['..', filesep, 'src']);

% Which set of data
cuCoarse = 3;
triangulationMethod = 'Del';
density = 'dense';
size = [1,1,1];

% Load data
filter = Filter(1,cuCoarse,density,size,triangulationMethod,'actual','actual',[],'zpos','undirected');
filter = filter.loadDTM('read');
filter = filter.generateNetwork();

%% Largest diameter against next largest

% Determine the largest radius and next largest radius for each node
idxTwoEdges = find(filter.nodeConnectivity>1);
nodeLargestRadius = zeros(length(idxTwoEdges),1);
nodeNextLargestRadius = zeros(length(idxTwoEdges),1);
for i = 1:length(idxTwoEdges)
    node = idxTwoEdges(i);
    edgeTo = filter.nodeEdgeTo(node,1:filter.nodeConnectivity(node));
    edgeToRadius = filter.edgeDiameter(edgeTo)/2;
    edgeToRadius = sort(edgeToRadius);
    nodeLargestRadius(i) = edgeToRadius(end);
    nodeNextLargestRadius(i) = edgeToRadius(end-1);
end

% Determine the Pearson correlation coefficient
R = corrcoef(nodeLargestRadius,nodeNextLargestRadius);
P = R(1,2);

% Plot
fig = figure;
set(gcf, 'WindowState', 'maximized');
scatter(nodeNextLargestRadius,nodeLargestRadius,'Marker','.','SizeData',200);
hold on;
p = min([max(xlim),max(ylim)]);
plot([0,p],[0,p],'LineWidth',Constants.lineWidth,'LineStyle','--','Color','k');
[~, objh] = legend({sprintf('Correlation coefficient = %.4f',P),'x = y'}, ...
    'Location','northwest','FontName','Helvetica','FontSize',Constants.fontSize);
objhl = findobj(objh, 'type', 'patch');
set(objhl, 'Markersize', 20, 'LineWidth', 5);
% text(0.72,0.7,'x = y','FontName','Helvetica','FontSize',Constants.fontSize);
% xl = xlim;
% yl = ylim;
% xl98 = xl(1)+0.98*(xl(2)-xl(1));
% yl05 = yl(1)+0.05*(yl(2)-yl(1));
% text(xl98,yl05,sprintf('Correlation coefficient = %.4f',P),'FontName','Helvetica','FontSize',Constants.fontSize,'HorizontalAlignment','right');
set(gca,'FontName','Helvetica','FontSize',Constants.fontSize);
xlabel('Next largest diameter','FontName','Helvetica','FontSize',Constants.fontSize);
ylabel('Largest diameter','FontName','Helvetica','FontSize',Constants.fontSize);

% Save
print(fig,[FileLocation.dataFolder, 'figures', filesep, 'SpatialCorrelation(a)'],'-dpng','-r300')

%% Largest diameter against volume

% Determine largest radius for each node
nodeLargestRadius = zeros(filter.nodeNum,1);
for i = 1:filter.nodeNum
    edgeTo = filter.nodeEdgeTo(i,1:filter.nodeConnectivity(i));
    edgeToRadius = filter.edgeDiameter(edgeTo)/2;
    nodeLargestRadius(i) = max(edgeToRadius);
end

% Determine the Pearson correlation coefficient
R = corrcoef(nodeLargestRadius,filter.nodeVolume);
P = R(1,2);

% Plot
fig = figure;
set(gcf, 'WindowState', 'maximized');
scatter(filter.nodeVolume,nodeLargestRadius,'Marker','.','SizeData',200);
[~, objh] = legend({sprintf('Correlation coefficient = %.4f',P)}, ...
    'Location','northwest','FontName','Helvetica','FontSize',Constants.fontSize);
objhl = findobj(objh, 'type', 'patch');
set(objhl, 'Markersize', 20, 'LineWidth', 5);
% xl = xlim;
% yl = ylim;
% xl98 = xl(1)+0.98*(xl(2)-xl(1));
% yl05 = yl(1)+0.05*(yl(2)-yl(1));
% text(xl98,yl05,sprintf('Correlation coefficient = %.4f',P),'FontName','Helvetica','FontSize',Constants.fontSize,'HorizontalAlignment','right');
set(gca, 'FontName','Helvetica','FontSize',Constants.fontSize);
xlabel('Node volume','FontName','Helvetica','FontSize',Constants.fontSize);
ylabel('Largest diameter','FontName','Helvetica','FontSize',Constants.fontSize);

% Save
print(fig,[FileLocation.dataFolder, 'figures', filesep, 'SpatialCorrelation(b)'],'-dpng','-r300')

%% Node volume against surface area

% Maximum surface area
maxNodeSurfaceArea = max(filter.nodeSurfaceArea);
maxRadius = sqrt(maxNodeSurfaceArea/(4*pi));

% Maximum volume
maxNodeVolume = max(filter.nodeVolume);

% Determine the line of volume to surface area for sphere
radius = linspace(0,maxRadius,1e3);
nodeSurfaceArea = 4*pi*radius.^2;
nodeVolume = (4/3)*pi*radius.^3;

% Determine the Pearson correlation coefficient
R = corrcoef(filter.nodeSurfaceArea,filter.nodeVolume);
P = R(1,2);

% Plot
fig = figure;
set(gcf, 'WindowState', 'maximized');
scatter(filter.nodeSurfaceArea,filter.nodeVolume,'Marker','.','SizeData',200);
hold on;
plot(nodeSurfaceArea,nodeVolume,'LineWidth',Constants.lineWidth,'LineStyle','--','Color','k');
[~, objh] = legend({sprintf('Correlation coefficient = %.4f',P),'Ratio of volume to surface area for spheres'}, ...
    'Location','northwest','FontName','Helvetica','FontSize',Constants.fontSize);
objhl = findobj(objh, 'type', 'patch');
set(objhl, 'Markersize', 20, 'LineWidth', 5);
% text(7,1.9,'Ratio of volume to surface area for spheres','FontName','Helvetica','FontSize',Constants.fontSize,'HorizontalAlignment','right');
xlim([0,maxNodeSurfaceArea]);
ylim([0,maxNodeVolume]);
% xl = xlim;
% yl = ylim;
% xl98 = xl(1)+0.98*(xl(2)-xl(1));
% yl05 = yl(1)+0.05*(yl(2)-yl(1));
% text(xl98,yl05,sprintf('Correlation coefficient = %.4f',P),'FontName','Helvetica','FontSize',Constants.fontSize,'HorizontalAlignment','right');
set(gca, 'FontName','Helvetica','FontSize',Constants.fontSize);
xlabel('Node surface area','FontName','Helvetica','FontSize',Constants.fontSize);
ylabel('Node volume','FontName','Helvetica','FontSize',Constants.fontSize);

% Save
print(fig,[FileLocation.dataFolder, 'figures', filesep, 'SpatialCorrelation(c)2'],'-dpng','-r300')

%% Node volume against surface area

% Determine idx of internal nodes
idxInternal = find(sum(filter.nodeBoundary,2)==0);

% Determine the surface area + edge areas
nodeInternalSurfaceArea = zeros(length(idxInternal),1);
for i = 1:length(idxInternal)
    node = idxInternal(i);
    nodeInternalSurfaceArea(i) = filter.nodeSurfaceArea(node);
    for j = 1:filter.nodeConnectivity(node)
        nodeInternalSurfaceArea(i) = nodeInternalSurfaceArea(i)+filter.edgeArea(filter.nodeEdgeTo(node,j));
    end
end

% Maximum surface area
maxNodeSurfaceArea = max(filter.nodeSurfaceArea);
maxRadius = sqrt(maxNodeSurfaceArea/(4*pi));

% Maximum volume
maxNodeVolume = max(filter.nodeVolume);

% Determine the line of volume to surface area for sphere
radius = linspace(0,maxRadius,1e3);
nodeSurfaceArea = 4*pi*radius.^2;
nodeVolume = (4/3)*pi*radius.^3;

% Determine the Pearson correlation coefficient
R = corrcoef(nodeInternalSurfaceArea,filter.nodeVolume(idxInternal));
P = R(1,2);

% Plot
fig = figure;
set(gcf, 'WindowState', 'maximized');
scatter(nodeInternalSurfaceArea,filter.nodeVolume(idxInternal),'Marker','.','SizeData',200);
hold on;
plot(nodeSurfaceArea,nodeVolume,'LineWidth',Constants.lineWidth,'LineStyle','--','Color','k');
[~, objh] = legend({sprintf('Correlation coefficient = %.4f',P),'Ratio of volume to surface area for spheres'}, ...
    'Location','northwest','FontName','Helvetica','FontSize',Constants.fontSize);
objhl = findobj(objh, 'type', 'patch');
set(objhl, 'Markersize', 20, 'LineWidth', 5);
% text(7,1.9,'Ratio of volume to surface area for spheres','FontName','Helvetica','FontSize',Constants.fontSize,'HorizontalAlignment','right');
xlim([0,maxNodeSurfaceArea]);
ylim([0,maxNodeVolume]);
% xl = xlim;
% yl = ylim;
% xl98 = xl(1)+0.98*(xl(2)-xl(1));
% yl05 = yl(1)+0.05*(yl(2)-yl(1));
% text(xl98,yl05,sprintf('Correlation coefficient = %.4f',P),'FontName','Helvetica','FontSize',Constants.fontSize,'HorizontalAlignment','right');
set(gca, 'FontName','Helvetica','FontSize',Constants.fontSize);
xlabel('Sum of node surface area and edge total area','FontName','Helvetica','FontSize',Constants.fontSize);
ylabel('Node volume','FontName','Helvetica','FontSize',Constants.fontSize);

% Save
print(fig,[FileLocation.dataFolder, 'figures', filesep, 'SpatialCorrelation(c)'],'-dpng','-r300')

%% Edge area against diameter

% Maximum area
maxEdgeArea = max(filter.edgeArea);
maxRadius = sqrt(maxEdgeArea/(pi));

% Determine the line of diameter to area for circle
radius = linspace(0,maxRadius,1e3);
edgeArea = pi*radius.^2;
diameter = 2*radius;

% Determine the Pearson correlation coefficient
R = corrcoef(filter.edgeArea,filter.edgeDiameter);
P = R(1,2);

% Plot
fig = figure;
set(gcf, 'WindowState', 'maximized');
scatter(filter.edgeArea,filter.edgeDiameter,'Marker','.','SizeData',200);
hold on;
plot(edgeArea,diameter,'LineWidth',Constants.lineWidth,'LineStyle','--','Color','k');
xlim([0,maxEdgeArea]);
[~, objh] = legend({sprintf('Correlation coefficient = %.4f',P),'Area for circles'}, ...
    'Location','northwest','FontName','Helvetica','FontSize',Constants.fontSize);
objhl = findobj(objh, 'type', 'patch');
set(objhl, 'Markersize', 20, 'LineWidth', 5);
set(gca, 'FontName','Helvetica','FontSize',Constants.fontSize);
xlabel('Edge total area','FontName','Helvetica','FontSize',Constants.fontSize);
ylabel('Edge diameter','FontName','Helvetica','FontSize',Constants.fontSize);

% Save
print(fig,[FileLocation.dataFolder, 'figures', filesep, 'SpatialCorrelation(d)'],'-dpng','-r300')

%% Edge radius against length

% Determine the Pearson correlation coefficient
R = corrcoef(filter.edgeDiameter,filter.edgeLength);
P = R(1,2);

% Plot
fig = figure;
set(gcf, 'WindowState', 'maximized');
scatter(filter.edgeLength,filter.edgeDiameter,'Marker','.','SizeData',200);
[~, objh] = legend({sprintf('Correlation coefficient = %.4f',P)}, ...
    'Location','northwest','FontName','Helvetica','FontSize',Constants.fontSize);
objhl = findobj(objh, 'type', 'patch');
set(objhl, 'Markersize', 20, 'LineWidth', 5);
set(gca, 'FontName','Helvetica','FontSize',Constants.fontSize);
xlabel('Edge length','FontName','Helvetica','FontSize',Constants.fontSize);
ylabel('Edge diameter','FontName','Helvetica','FontSize',Constants.fontSize);

% Save
print(fig,[FileLocation.dataFolder, 'figures', filesep, 'SpatialCorrelation(e)'],'-dpng','-r300')

%% Edge reduced against radius

% Determine the Pearson correlation coefficient
R = corrcoef(filter.edgeReduced,filter.edgeLength);
P = R(1,2);

% Plot
fig = figure;
set(gcf, 'WindowState', 'maximized');
scatter(filter.edgeReduced,filter.edgeDiameter,'Marker','.','SizeData',200);
[~, objh] = legend({sprintf('Correlation coefficient = %.4f',P)}, ...
    'Location','northwest','FontName','Helvetica','FontSize',Constants.fontSize);
objhl = findobj(objh, 'type', 'patch');
set(objhl, 'Markersize', 20, 'LineWidth', 5);
set(gca, 'FontName','Helvetica','FontSize',Constants.fontSize);
xlabel('Edge reduced area','FontName','Helvetica','FontSize',Constants.fontSize);
ylabel('Edge diameter','FontName','Helvetica','FontSize',Constants.fontSize);

% Save
print(fig,[FileLocation.dataFolder, 'figures', filesep, 'SpatialCorrelation(f)'],'-dpng','-r300')

%% Connectivity

% nodeConnectivity1 = nodeConnectivity+max(nodeBoundary,[],2);
% figure; 
% set(gcf, 'WindowState', 'maximized');
% histogram(nodeConnectivity1,min(nodeConnectivity1)-0.5:1:max(nodeConnectivity1)+0.5);
% set(gca, 'FontName','Helvetica','FontSize',fontSize);
% xlabel('Node connectivity','FontName','Helvetica','FontSize',fontSize);
% ylabel('Frequency','FontName','Helvetica','FontSize',fontSize);