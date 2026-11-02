classdef Base
    % The class is used to generate any property relating to the base
    % particles
    
    properties
        % Parameters
        baseType
        cuFine
        d50
        minFine
        maxFine
        
        % PSDs
        psdFineVol
        psdFineNum
        psdFineBucket
    end
    
    methods
        % ===== Initialisation ============================================
        % Creation of an object of this class
        function [obj] = Base()
        end
        
        % ===== Creation of bases =========================================
        % Linearly distributed base by volume
        function [obj] = linearVolBase(obj, d50, cuFine)
            % Input base parameters
            obj.d50 = d50;
            obj.cuFine = cuFine;
            obj.baseType = 'linearVol';
            
            % Determine max and min particle size
            [obj] = obj.determineFinesParameters();
            
            % Determine the PSDs
            [obj] = obj.generatePSDs();
        end
        
        % Linearly distributed base by number
        function [obj] = linearNumBase(obj, minFine, maxFine)
            % Input base parameters
            obj.minFine = minFine;
            obj.maxFine = maxFine;
            obj.d50 = [];
            obj.cuFine = [];
            obj.baseType = 'linearNum';
            
            % Determine the PSDs
            [obj] = obj.generatePSDs();
        end
        
        % ===== Saving mechanism ==========================================
        % Re-generate PSDs
        function [obj] = generatePSDs(obj)
            if strcmp(obj.baseType, 'linearVol')
                % Linearly graded base (by vol/mass)
                obj = obj.determinePSDVol();
                
                % PSD by number
                obj.psdFineNum = Statistics.convertPSDVolToNum(obj.psdFineVol);
            else
                % Create PSD by num
                obj.psdFineNum = NaN(10001,2);
                obj.psdFineNum(:,1) = linspace(0,100,10001);
                obj.psdFineNum(:,2) = linspace(obj.minFine,obj.maxFine,10001);
                
                % PSD by vol is not required in simulation
                obj.psdFineVol = [];
            end
            
            % PSD bucket
            obj = obj.determinePSDBucket();
        end
        
        % Remove PSDs
        function [obj] = removePSDs(obj)
            obj.psdFineBucket = [];
            obj.psdFineNum = [];
            obj.psdFineVol = [];
        end
    end
    
    methods (Access = private)
        % Determine the maximum and minimum base diameters given d50B
        function [obj] = determineFinesParameters(obj)
            % The function determines fine_min and fine_max based on Cu,base and d50/D0
            m = 50/log(obj.cuFine);
            if m==Inf
                obj.minFine = obj.d50;
                obj.maxFine = obj.d50;
            else
                c = 50-(log(obj.d50)*m);
                obj.minFine = exp((0-c)/m);
                obj.maxFine = exp((100-c)/m);
            end
            assert(obj.minFine<=obj.d50);
            assert(obj.maxFine>=obj.d50);
        end
        
        % For linearly graded base by vol
        function [obj] = determinePSDVol(obj)
            % Determine psdFineVol
            
            % Clear
            obj.psdFineVol = [];
            
            % Calculate volumes
            vol0 = (obj.minFine.^3) * pi * (4/3);
            volfinemax= (obj.maxFine.^3) * pi * (4/3);
            
            % Calculate number of particles
            vol_part_fine = logspace(log10(vol0),log10(volfinemax),10001);
            dist = vol_part_fine ./ (pi * (4/3));
            dist = dist.^(1/3);
            for k = 1 : 10000
                num1 = k;
                num2 = k + 1;
                obj.psdFineVol(num1:num2) = logspace(log10(dist(num1)),log10(dist(num2)),(num2-num1+1));
            end
            obj.psdFineVol = obj.psdFineVol';
            obj.psdFineVol(:,2) = obj.psdFineVol;
            obj.psdFineVol(:,1) = 0:0.01:100;
        end
        
        % Create psdFineBucket, which is used to save the PSDs
        function [obj] = determinePSDBucket(obj)
            % Clear
            obj.psdFineBucket = [];
            
            % Determine psdFineBucket
            erod_ind = 1:100:length(obj.psdFineVol);
            obj.psdFineBucket = obj.psdFineVol(erod_ind,2);
        end
    end
    
    methods (Static)
        % Determine d50B, given (1) D15F, (2) the ratio of D15F to d85B,
        % and (3) cuFine
        function [d50] = givenRatio(r, D15, cuFine)
            % Gradient of line
            m = 50/log(cuFine);
            
            % Determine d85
            d85 = D15/r;
            
            % Determine d50
            if m==Inf
                d50 = d85;
            else
                c = 85-(log(d85)*m);
                d50 = exp((50-c)/m);
            end
            assert(d50<=d85)
        end
    end
end