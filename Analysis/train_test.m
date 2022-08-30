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
        function obj = train_test(mat, label, ratio, split)
            %NOTE: Any empty input argument is left as a []
            
            
        end
        
        function obj = remove_nan(obj)
            % Function which removes all rows containing NaNs from the
            % feature matrix and their corresponding rows in the label
            % vector
            
            obj.mat = 2;
            
        end
        
        function obj = new_ratio(obj)
            % Takes the given ratio of each label and returns what the
            % sampling should be
            obj.ratio = [NaN];
            
        end
        
        
        
        function obj = untitled(inputArg1,inputArg2)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            obj.Property1 = inputArg1 + inputArg2;
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

