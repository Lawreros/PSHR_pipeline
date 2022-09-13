function [ret] = rmssd_calc(mat,bin,band)
% Calculates the root mean square of successive differences (RMSSD) for
% the vector of mesasurements provided.
% Inputs:
%   mat: A [n-by-1] vector which contains the data you want to
%       calculate RMSSD for
%
%   bin: [1-by-2 cell array] Used for creating a vector of the
%       RMSSD results from a sliding bin of Y seconds or entries. This takes the
%       format of {[before, after], 'units'}, so if you want to have a bin 
%       of 5 seconds before (including current RR-interval) and 3 seconds 
%       after: {[5,3], 'second'} or if you want the 5 entries before and 
%       3 entries after the index: {[5,3], 'measure'}
%       If you don't want this, set bin to false.
%
%   band: [2 int vector] The range [start, end] of values you want
%       to calculate the RMSSD of. If false, then analyze the whole
%       range
    
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
    
    %If they've decided to use the bin values
    if iscell(bin)
        ret = NaN(r_2-(r_1-1),1);
        bin_gen = bin_get(mat(r_1:r_2), bin);
            
        for i = 1:r_2-(r_1-1)
            summation = 0;
            dump = bin_gen.next;
                
            if length(dump) > 1
                for k = 2:length(dump)
                    summation = summation+(dump(k,1)-dump(k-1,1))^2;
                end
                ret(i) =  sqrt((1/(length(dump)-1))*summation);
            end
            
        end
        
    else
        % If they just want a single value for the input vector
        summation = 0;
        rem_ent = 0;
        for  i = (r_1+1):r_2
            if isnan(mat(i,1)-mat(i-1,1)) % Check if NaNs are present in the data and skip them
                rem_ent = rem_ent + 1;
            else
                summation = summation+(mat(i,1)-mat(i-1,1))^2;
            end
        end
        ret = sqrt(1/((r_2-rem_ent)-r_1)*summation);
    end

end
