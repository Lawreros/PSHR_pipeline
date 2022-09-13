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
        
        % MATLAB classes require you to return the obj if you make changes
        % to it, so I made it a seperate method
        function obj = itter(obj)
            obj.itt = obj.itt + 1;
        end
        
        function nxt = next(obj)
            % Generator which keeps returning the next bin sample until
            % there are no more and a 'false' value is returned. Returns
            % NaNs if various conditions are met.
            
            obj.itter();
            
            if obj.itt > length(obj.mat)
                disp('Class is exhausted')
                nxt = false;
            else
                nxt = obj.sample(obj.itt);
            end
            
        end
        
        function samp = sample(obj, idx)
            % Returns the binned values for a particular index in the
            % initial mat file
            
            % Input:
            %   idx: [int] the index number that you which to take a bin in
            %   reference to
            %
            % Output:
            %   samp: [n-by-1 matrix] 
            
            cap = length(obj.mat);
            
            if strcmp(obj.units, 'second') %If they want to use seconds as their units
                
                i = 0;
              
                % Find how far back to go until the sum is greater than the bin
                % number of seconds
                if obj.before > 0
                    while (sum(obj.mat(idx-i:idx))/1000) <= obj.before
%                         disp(strcat('backsum = ',string(sum(obj.mat(idx-i:idx))/1000)));
                        i = i+1;
                        if i == idx
                            i = i+1;
                            break;
                        end
                    end
                    i = i-1; % the while loop goes until i is one too many, so just subtract 1 for the bounds
                else
                    i = 0;
                end

                
                % Find how far forward to go until the sum is greater than the
                % bin number of seconds
                j = 1; %j starts at 1 so that it looks at the RR-intervals after idx without including it                
                
                if obj.after > 0
                    try
                        while (sum(obj.mat(idx+1:idx+j))/1000) <= obj.after
%                              disp(strcat('frontsum = ',string(sum(obj.mat(idx+1:idx+j))/1000)));
                             j = j+1;
                            if idx+j > cap
                                j = j+1;
                                disp(string(j));
                                break;
                            end
                        end
                        j = j-1;
                    catch
                        j = cap;
                    end
                else
                    j = 0;
                end
                
                
                % If there are NaNs in the bin, then just return NaNs
                % TODO: Maybe add option to ignore the NaNs
                if i < idx && idx+j <= cap && isnan(sum(obj.mat(idx-i:idx+j)))
                    i = false;
                    j = false;
                end
                
                % If there is more than this single data point
                if i < idx && idx+j <= cap %&& (i+j)>1
                    samp = obj.mat(idx-i:idx+j);
                else
                    samp = NaN;
                end
                
                
            elseif strcmp(obj.units, 'measure')
                
                if (idx-obj.before) >= 1 && (idx+obj.after) <= cap
                    samp = obj.mat(idx - obj.before: idx+obj.after);
                else
                    samp = NaN;
                end
            end
            
        end
        
    end
end

