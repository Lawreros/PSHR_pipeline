function [ret] = sdnn_calc(mat,bin,band)
% Calculates the standard deviation for the vector of mesasurements provided.
% Inputs:
%   mat: A [n-by-1] vector which contains the data you want to
%       calculate standard deviation for
     
%   bin: [1-by-2 cell array] Used for creating a vector of the
%       standard deviations results from a sliding bin of Y seconds or entries. 
%       This takes the format of {[before, after], 'units'}, so if you want 
%       to have a bin of 5 seconds before (including current RR-interval) 
%       and 3 seconds after: {[5,3], 'second'} or if you want the 5 entries
%       before and 3 entries after the index: {[5,3], 'measure'}
%       If you don't want this, set bin to false.
%
%   band: [2 int vector] The range [start, end] of values you want
%       to calculate the RMSSD of. If false, then analyze the whole range
    
% Returns:
%   ret: Either a [n-by-1] vector containing the results (if bin is not
%       false) or an [int] if bin is false.
    
    if band
        r_1 = band(1);
        r_2 = band(2);
    else
        [r_2,c] = size(mat);
        r_1 = 1;
    end
    
    %If they've decided to use the bin values
    if iscell(bin)
        ret = NaN(r_2-(r_1-1),1);
        bin_gen = bin_get(mat(r_1:r_2), bin);
         
        for i = 1:r_2-(r_1-1)  
            dump = bin_gen.next;
                
            if length(dump) > 1
                ret(i) = std(dump);
            end
        end
       
    else
        % If they just want a single value for the input vector
        ret = nanstd(mat(r_1:r_2,1));
    end

end