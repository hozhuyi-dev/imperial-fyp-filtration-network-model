classdef Plot
    % The class aims to create plots that are consistent appropriate for
    % report
    
    properties (Constant)
        % Font properties
        fontName = 'Helvetica';
        fontSizeHalf = 48;
        fontSizeFull = 12;
        
        % Line properties
        lineWidthHalf = 5;
        lineWidthFull = 2;
        
        % Marker properties
        
        % Colours
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
    end
    
    methods (Static)
        % Given an array
        function [fig] = arrayColumn(A, varargin)
            % Input parser
            parser = inputParser;
            addParameter(parser, 'legendList', {})
            addParameter(parser, 'legendLocation', 'northeast');
            addParameter(parser, 'xLabelString', '')
            addParameter(parser, 'yLabelString', '')
            addParameter(parser, 'xl', [])
            addParameter(parser, 'plotOption', 'full')
            parse(parser, varargin{:})
            legendList = parser.Results.legendList;
            legendLocation = parser.Results.legendLocation;
            xLabelString = parser.Results.xLabelString;
            yLabelString = parser.Results.yLabelString;
            plotOption = parser.Results.plotOption;
            
            % Font size
            switch plotOption
                case 'full'
                    fontSize = Plot.fontSizeFull;
                    lineWidth = Plot.lineWidthFull;
                case 'half'
                    fontSize = Plot.fontSizeHalf;
                    lineWidth = Plot.lineWidthHalf;
            end
            
            % Number of rows
            rowNum = size(A,1);
            
            % Plot
            fig = figure;
            set(gcf,'WindowState','maximized');
            for i = 1:size(A,2)
                % Plot
                plot(1:rowNum,A(:,i),'LineWidth',lineWidth);
                hold on;
            end
            legend(legendList,'Location',legendLocation,'FontName',Plot.fontName,'FontSize',fontSize);
            if ~isempty(xl)
                xlim(xl);
            end
            set(gca, 'FontName','Helvetica','FontSize',16);
            xlabel(xLabelString,'FontName',Plot.fontName,'FontSize',fontSize);
            ylabel(yLabelString,'FontName',Plot.fontName,'FontSize',fontSize);
        end
        
        % Given a struct
        
        % Save figures
        
    end
end