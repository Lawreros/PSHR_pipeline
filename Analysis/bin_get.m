classdef bin_get < handle
    %BIN_GET: Class which serves as a pseudo-generator in MATLAB for
    %sampling data and returning binned vectors    
    
    properties
        mat
        itt = 0
        units = 'second'
        before = 0
        after = 0
    end
    
    methods
        
        % CLASS CONSTRUCTOR
        function obj = bin_get(mat, bin)
            % Load the given vector of values into the class
            obj.mat = mat;
            
            if iscell(bin)
                a = bin{1};
                b = bin{2};
                
                if strcmp(b, 'second')
                    obj.units = b;
                    obj.before = a(1);
                    obj.after = a(2);
                elseif strcmp(b, 'measure')
                    obj.units = b;
                    obj.before = a(1);
                    obj.after = a(2);
                end
            end
            
        end
        
        function obj = itter(obj)
            obj.itt = obj.itt + 1;
        end
        
        function nxt = next(obj)
            % Generator which keeps returning the next bin sample until
            % there are no more
            
            nxt = obj.sample(obj.itt);
            obj.itter()
            
        end
        
        function samp = sample(obj, idx)
            % Returns the binned values for a particular index in the
            % initial mat file
            
            if strcmp(obj.units, 'second')
                
                while (sum(obj.mat(idx-i:idx))/1000) <= obj.before
                    i = i+1;
                    if i == idx
                        break;
                    end
                end
                
                while (sum(obj.mat(idx:idx+j))/1000) <= obj.after
                    j = j+1;
                    if j == idx
                        break;
                    end
                end
                
                samp = obj.mat(idx-i:idx+j);
                
                
            elseif strcmp(obj.units, 'measure')
                samp = obj.mat(idx - obj.before: idx+obj.after);
            end
            
        end
        
    end
end

