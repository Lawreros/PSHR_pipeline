function [ret] = pnnx_calc(mat,diff,bin,band)
% Calculates the percentage of adjacent NN-intervals that differ from
% each other by more than "diff" milliseconds
% Inputs:
%   mat: A [n-by-1] vector which contains the data you want to
%       calculate pNNX for
%
%   diff: [int] The minimum difference in milliseconds between
%       successive NN-intervalse that you want to count
%
%   bin: [1-by-2 cell array] Used for creating a vector of the pNNX
%       results from a sliding bin of Y seconds or entries. This takes the
%       format of {index, 'units'}, so if you want to have a bin of the
%       last 5 seconds: {5, 'second'} or if you want the last 5 measurements: {5, 'measure'}
%       If you don't want this, set bin to false.
%
%   band: [2 int vector] The range [start, end] of values you want
%       to calculate the pnnX of. If false, then analyze the whole
%       range of mat
%
% Returns:
%   ret: Either a [n-by-1] vector containing the results (if bin is
%       not false) or an [int] if bin is false.
    
    if band
        r_1 = band(1);
        r_2 = band(2);
    else
        [r_2,c] = size(mat);
        r_1 = 1;
    end
    
    % If they've decided to use the bin values
    if iscell(bin)
        a = bin{1}; % value
        b = bin{2}; % units
        
        ret = NaN(r_2-r_1,1);
        
        if strcmp(b,'second')
            for i = r_1:r_2
                count = 0;
                j = 0;
                
                % Check loop backward until you have reached the 'b'
                % seconds in the past through summing
                while (sum(mat(i-j:i,1))/1000) <= a
                    j = j+1;
                    if j == i
                        break;
                    end
                end
                
                if j < i && isnan(mat(i-j,1)) % Simple catch to prevent premature binning due to NaNs
                   j = 0;            % The while loop will return false if the sum is NaN,
                end                  % making the binning inconsistent
                
                
                if j > 1 && j < i % If there is more than one entry
                    for k = (i-j+2):i
                        if abs(mat(k,1) - mat(k-1,1))>= diff
                            count = count+1;
                        end
                    end
                    ret(i-r_1+1,1) = count/(j-1);
                else
                    ret(i-r_1+1,1) = NaN;
                end
        
            end
        else % Looking at the past 'b' entries for the calculation
            for i = (a+1):(r_2-r_1+1)
                count = 0;
                for j = 0:(a-1)
                    if abs(mat(i-j,1) - mat(i-j-1,1))>= diff
                        count = count+1;
                    end
                end
                ret(i,1) = count/a;
            end
        end
        
    else
        % If they just want a percentage for a matrix
        count = 0;
        for i = (r_1+1):r_2
            if abs(mat(i,1)-mat(i-1,1))>= diff
                count = count+1;
            end
        end
        ret = count/(r_2-r_1);
    end
end

