function [ret] = rmssd_calc(mat,bin,band)
    % Calculates the root mean square of successive differences (RMSSD) for
    % the vector of mesasurements provided.
    %   Inputs:
    %       mat: A [n-by-1] vector which contains the data you want to
    %       calculate RMSSD for
    %
    %       bin: [1-by-2 cell array] Used for creating a vector of the
    %       RMSSD
    %       results from a sliding bin of Y seconds or entries. This takes the
    %       format of {index, 'units'}, so if you want to have a bin of the
    %       last 5 seconds: {5, 'second'} or if you want the last 5 measurements: {5, 'measure'}
    %       If you don't want this, set bin to false.
    %
    %       band: [2 int vector] The range [start, end] of values you want
    %       to calculate the RMSSD of. If false, then analyze the whole
    %       range
    
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
                summation = 0;
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
                    for  k = (i-j+2):i
                        summation = summation+(mat(k,1)-mat(k-1,1))^2;
                    end
                    ret(i-r_1+1,1) = sqrt((1/(j-2))*summation);
                else
                    ret(i-r_1+1,1) = NaN;
                end
            end
        else % Looking at the past 'a' entries for the calculation
            for i = (a+1):(r_2-r_1+1)
                summation = 0;
                for j = 0:(a-1)
                    summation = summation+(mat(i-j,1)-mat(i-j-1,1))^2;
                end
                ret(i,1) = sqrt((1/a)*summation);
            end
        end
        
    else
        % If they just want a single value for the input vector
        summation = 0;
        for  i = (r_1+1):r_2
            summation = summation+(mat(i,1)-mat(i-1,1))^2;
        end
        ret = sqrt(1/(r_2-r_1)*summation);
    end

end
