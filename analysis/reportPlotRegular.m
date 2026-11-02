% Generate plots to show connectivity

clear;
clc;
close all;

edgeSize = 10;
nodeSize = 2e3;
az = 20;
el = 10;
axisSize = 22;
zSize = 30;

%% Square

nodeCoord = [0, 0, 0; 1, 0, 0; -1, 0, 0; 0, 0, 1; 0, 0, -1];
edgeList = ones(size(nodeCoord, 1)-1, 2);
edgeList(:,2) = 2:size(nodeCoord, 1);
nodeNum = size(nodeCoord, 1);
edgeNum = size(edgeList, 1);

figure
set(gcf, 'WindowState', 'maximized')
% The edges are represented as lines
for i = 1:edgeNum
    plot3(nodeCoord(edgeList(i,:), 1), ...
        nodeCoord(edgeList(i,:), 2), ...
        nodeCoord(edgeList(i,:), 3), ...
        'lineWidth', edgeSize, 'Color', zeros(3,1))
    hold on
end
% The nodes are represented with markers.
scatter3(nodeCoord(:,1), nodeCoord(:,2), nodeCoord(:,3), nodeSize*ones(nodeNum, 1), 'k', 'filled')
hold off
set(gca, 'DataAspectRatio', [1,1,1])
grid on
xticks(-1:0.5:1)
yticks(-1:0.5:1)
zticks(-1:0.5:1)
xlim([-1, 1])
ylim([-1, 1])
zlim([-1, 1])
view([az,el])
set(gca, 'FontSize', axisSize);
zlabel('Flow direction', 'FontSize', zSize);

%% Slanted Square

x = cosd(45);
nodeCoord = [0, 0, 0; x, 0, x; x, 0, -x; -x, 0, x; -x, 0, -x];
edgeList = ones(size(nodeCoord, 1)-1, 2);
edgeList(:,2) = 2:size(nodeCoord, 1);
nodeNum = size(nodeCoord, 1);
edgeNum = size(edgeList, 1);

figure
set(gcf, 'WindowState', 'maximized')
% The edges are represented as lines
for i = 1:edgeNum
    plot3(nodeCoord(edgeList(i,:), 1), ...
        nodeCoord(edgeList(i,:), 2), ...
        nodeCoord(edgeList(i,:), 3), ...
        'lineWidth', edgeSize, 'Color', zeros(3,1))
    hold on
end
% The nodes are represented with markers.
scatter3(nodeCoord(:,1), nodeCoord(:,2), nodeCoord(:,3), nodeSize*ones(nodeNum, 1), 'k', 'filled')
hold off
set(gca, 'DataAspectRatio', [1,1,1])
grid on
xticks(-1:0.5:1)
yticks(-1:0.5:1)
zticks(-1:0.5:1)
xlim([-1, 1])
ylim([-1, 1])
zlim([-1, 1])
view([az,el])
set(gca, 'FontSize', axisSize);
zlabel('Flow direction', 'FontSize', zSize);

%% Cubic

x = 1;
nodeCoord = [0, 0, 0; 0, 0, x; 0, x, 0; x, 0, 0; 0, 0, -x; 0, -x, 0; -x, 0, 0];
edgeList = ones(size(nodeCoord, 1)-1, 2);
edgeList(:,2) = 2:size(nodeCoord, 1);
nodeNum = size(nodeCoord, 1);
edgeNum = size(edgeList, 1);

figure
set(gcf, 'WindowState', 'maximized')
% The edges are represented as lines
for i = 1:edgeNum
    plot3(nodeCoord(edgeList(i,:), 1), ...
        nodeCoord(edgeList(i,:), 2), ...
        nodeCoord(edgeList(i,:), 3), ...
        'lineWidth', edgeSize, 'Color', zeros(3,1))
    hold on
end
% The nodes are represented with markers.
scatter3(nodeCoord(:,1), nodeCoord(:,2), nodeCoord(:,3), nodeSize*ones(nodeNum, 1), 'k', 'filled')
hold off
set(gca, 'DataAspectRatio', [1,1,1])
grid on
xticks(-1:0.5:1)
yticks(-1:0.5:1)
zticks(-1:0.5:1)
xlim([-1, 1])
ylim([-1, 1])
zlim([-1, 1])
view([az,el])
set(gca, 'FontSize', axisSize);
zlabel('Flow direction', 'FontSize', zSize);

%% Slanted Cubic

a = atand(1/sqrt(2));
nodeCoord = [0, 0, 0; 0, cosd(a), sind(a); cosd(a)*cosd(30), -cosd(a)*sind(30), sind(a); -cosd(a)*cosd(30), -cosd(a)*sind(30), sind(a); cosd(a)*cosd(30), cosd(a)*sind(30), -sind(a); -cosd(a)*cosd(30), cosd(a)*sind(30), -sind(a); 0, -cosd(a), -sind(a)];
edgeList = ones(size(nodeCoord, 1)-1, 2);
edgeList(:,2) = 2:size(nodeCoord, 1);
nodeNum = size(nodeCoord, 1);
edgeNum = size(edgeList, 1);

figure
set(gcf, 'WindowState', 'maximized')
% The edges are represented as lines
for i = 1:edgeNum
    plot3(nodeCoord(edgeList(i,:), 1), ...
        nodeCoord(edgeList(i,:), 2), ...
        nodeCoord(edgeList(i,:), 3), ...
        'lineWidth', edgeSize, 'Color', zeros(3,1))
    hold on
end
% The nodes are represented with markers.
scatter3(nodeCoord(:,1), nodeCoord(:,2), nodeCoord(:,3), nodeSize*ones(nodeNum, 1), 'k', 'filled')
hold off
set(gca, 'DataAspectRatio', [1,1,1])
grid on
xticks(-1:0.5:1)
yticks(-1:0.5:1)
zticks(-1:0.5:1)
xlim([-1, 1])
ylim([-1, 1])
zlim([-1, 1])
view([az,el])
set(gca, 'FontSize', axisSize);
zlabel('Flow direction', 'FontSize', zSize);

%% BCC

b = atand(1/sqrt(2));
a = sind(b);
a2 = cosd(b)*cosd(45);
nodeCoord = [0, 0, 0; a2, a2, a; a2, -a2, a; -a2, a2, a; -a2, -a2, a; a2, a2, -a; a2, -a2, -a; -a2, a2, -a; -a2, -a2, -a];
edgeList = ones(size(nodeCoord, 1)-1, 2);
edgeList(:,2) = 2:size(nodeCoord, 1);
nodeNum = size(nodeCoord, 1);
edgeNum = size(edgeList, 1);

figure
set(gcf, 'WindowState', 'maximized')
% The edges are represented as lines
for i = 1:edgeNum
    plot3(nodeCoord(edgeList(i,:), 1), ...
        nodeCoord(edgeList(i,:), 2), ...
        nodeCoord(edgeList(i,:), 3), ...
        'lineWidth', edgeSize, 'Color', zeros(3,1))
    hold on
end
% The nodes are represented with markers.
scatter3(nodeCoord(:,1), nodeCoord(:,2), nodeCoord(:,3), nodeSize*ones(nodeNum, 1), 'k', 'filled')
hold off
set(gca, 'DataAspectRatio', [1,1,1])
grid on
xticks(-1:0.5:1)
yticks(-1:0.5:1)
zticks(-1:0.5:1)
xlim([-1, 1])
ylim([-1, 1])
zlim([-1, 1])
view([az,el])
set(gca, 'FontSize', axisSize);
zlabel('Flow direction', 'FontSize', zSize);

%% FCC

b = 45;
a = sind(b);
nodeCoord = [0, 0, 0; a, 0, a; 0, a, a; -a, 0, a; 0, -a, a; a, 0, -a; 0, a, -a; -a, 0, -a; 0, -a, -a; a, a, 0; a, -a, 0; -a, a, 0; -a, -a, 0];
edgeList = ones(size(nodeCoord, 1)-1, 2);
edgeList(:,2) = 2:size(nodeCoord, 1);
nodeNum = size(nodeCoord, 1);
edgeNum = size(edgeList, 1);

figure
set(gcf, 'WindowState', 'maximized')
% The edges are represented as lines
for i = 1:edgeNum
    plot3(nodeCoord(edgeList(i,:), 1), ...
        nodeCoord(edgeList(i,:), 2), ...
        nodeCoord(edgeList(i,:), 3), ...
        'lineWidth', edgeSize, 'Color', zeros(3,1))
    hold on
end
% The nodes are represented with markers.
scatter3(nodeCoord(:,1), nodeCoord(:,2), nodeCoord(:,3), nodeSize*ones(nodeNum, 1), 'k', 'filled')
hold off
set(gca, 'DataAspectRatio', [1,1,1])
grid on
xticks(-1:0.5:1)
yticks(-1:0.5:1)
zticks(-1:0.5:1)
xlim([-1, 1])
ylim([-1, 1])
zlim([-1, 1])
view([az,el])
set(gca, 'FontSize', axisSize);
zlabel('Flow direction', 'FontSize', zSize);