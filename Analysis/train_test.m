classdef train_test
    %train_test: This class serves to organize a matrix containing relevant
    %features which is labeled into portions conducive for use in several
    %analyses.
    
    properties
        mat
        label
        ratio
        split = 0.8
        ignore_nan = true
    end
    
    methods
        
        % Relevant functions:
        %   Load data
        %   Organize data by both ratio and split
        %   Return either some or all of the different categories
        %   Support different organizational splits
        %   Contain accuracy tests
        
        % CONSTRUCTOR FOR CLASS
        function obj = train_test(mat, label, ratio, varargin)
            %NOTE: Any empty input argument is left as a []
            
            obj.mat = mat;
            obj.label = label;
            obj.ratio = ratio;
            
            
            p = inputParser;
            addParameter(p, 'split', 0.8, @isscalar);
            addParameter(p, 'ignore_nan', true, @islogical);
            parse(p,varargin{:});
            
            if p.Results.ignore_nan
                obj = obj.remove_nan();
            end
            
        end
        
        function obj = remove_nan(obj)
            % Function which removes all rows containing NaNs from the
            % feature matrix and their corresponding rows in the label
            % vector
            obj.label(any(isnan(obj.label),2),:)=[];
            obj.mat(any(isnan(obj.label),2),:)=[];
        
                
            obj.label(any(isnan(obj.mat),2),:) = [];
            obj.mat(any(isnan(obj.mat),2),:) = [];
        end
        
        function [catgs, locs, quant] = get_categories(obj)
            %returns categories, their locations, and their quantities.
            %NaNs do not work with this, so there can be no NaN labels.
            [catgs, locs, kk] = unique(obj.label);
            quant = histc(kk, 1:numel(locs));
            
            if length(obj.ratio) ~= length(quant)
                disp(strcat('ERROR: Provided ratios for categories: ', string(length(obj.ratio)),...
                ' does not match number of categories present: ',string(length(quant))));
                return
            end
        end
        
        function num_cat = new_ratio(obj)
            % Takes the given ratio of each label and returns what the
            % sampling should be
            [catgs, locs, quant] = obj.get_categories;
            req_change = 2;
            num_cat = NaN;
            
            while req_change >= 1
                [val idx] = min(quant);
                
                % find what all other values should be from ratio
                for i = 1:length(obj.ratio)
                    num_cat(i) = obj.ratio(i)/obj.ratio(idx) * val;
                end
                
                [req_change, idz] = max(num_cat - quant);
                quant(idx) = quant(idx) - (req_change * obj.ratio(idx)/obj.ratio(idz));
            end
            
            % floor new_f so you don't deal with float issues, there is no way
            % rounding down will require you to access more data
            num_cat = floor(num_cat);
        end
        
        
        
        function obj = rand_samp(obj, type)
            % Generates the random selection of rows/labels for testing and
            % training sets
            
            %type = KFold, split, 
            
            
        end
        
        function [train_mat, train_label] = get_train(obj)
            
        end
        
        function [test_mat, test_label] = get_test(obj)
            
        end
        
        function accuracy(obj, pred, target)
            % Function, when fed
            
            
        end
        
    end
end

