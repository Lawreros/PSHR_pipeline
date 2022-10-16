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
%       format of {[before, after], 'units'}, so if you want to have a bin 
%       of 5 seconds before (including current RR-interval) and 3 seconds 
%       after: {[5,3], 'second'} or if you want the 5 entries before and 
%       3 entries after the index: {[5,3], 'measure'}
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
        ret = NaN(r_2-(r_1-1),1);
        bin_gen = bin_get(mat(r_1:r_2), bin);
        
        for i = 1:r_2-(r_1-1)
            count = 0;
                
            dump = bin_gen.next;
                
            if length(dump) > 1
                for k = 2:length(dump)
                    if abs(dump(k) - dump(k-1)) >= diff
                        count = count+1;
                    end
                end
                ret(i) = count/(length(dump)-1);
            else
                % do nothing because NaN if dump = []
            end
        end
    else
        % If they just want a percentage for a matrix
        count = 0;
        rem_ent = 0;
        for i = (r_1+1):r_2
            if isnan(mat(i,1)-mat(i-1,1)) % Check if NaNs are present in the data and skip them
                rem_ent = rem_ent + 1;
            elseif abs(mat(i,1)-mat(i-1,1))>= diff
                count = count+1;
            end
        end
        ret = count/((r_2-r_1)-rem_ent);
    end
end

