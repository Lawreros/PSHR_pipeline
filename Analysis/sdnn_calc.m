function [ret] = sdnn_calc(mat,bin,band)
    % Calculates the standard deviation for the vector of mesasurements provided.
    %   Inputs:
    %       mat: A [n-by-1] vector which contains the data you want to
    %       calculate standard deviation for
    %       bin: [1-by-2 cell array] Used for creating a vector of the
    %       standard deviations
    %       results from a sliding bin of Y seconds or entries. This takes the
    %       format of {index, 'units'}, so if you want to have a bin of the
    %       last 5 seconds: {5, 'second'} or if you want the last 5 measurements: {5, 'measure'}
    %       If you don't want this, set bin to false.
    %
    %       band: [2 int vector] The range [start, end] of values you want
    %       to calculate the RMSSD of. If false, then analyze the whole
    %       range
    
    %   Returns:
    %       ret: Either a [n-by-1] vector containing the results (if bin is
    %       not false) or an [int] if bin is false.
    
    if band
        r_1 = band(1);
        r_2 = band(2);
    else
        [r_2,c] = size(mat);
        r_1 = 1;
    end
    
    %If they've decided to use the bin values
    if iscell(bin)
        a = bin{1}; % value
        b = bin{2}; % units
        
        ret = NaN(r_2-r_1,1);
        
        if strcmp(b,'second')
            for i = r_1:r_2
                j = 0;
                
                % Check loop backward until you have reached the 'b'
                % seconds in the past through summing
                while (sum(mat(i-j:i,1))/1000) <= a
                    j = j+1;
                    if j == i
                        break;
                    end
                end
                
                if j > 1 && j < i % If there is more than one entry
                    ret(i-r_1+1,1) = std(mat((i-j+1:i),1));
                else
                    ret(i-r_1+1,1) = NaN;
                end
            end
        else % Looking at the past 'a' entries for the calculation
            for i = (a+1):(r_2-r_1+1)
                ret(i,1) = std(mat(i-a:i,1));
            end
        end
        
    else
        % If they just want a single value for the input vector
        ret = std(mat(r_1:r_2,1));
    end

end