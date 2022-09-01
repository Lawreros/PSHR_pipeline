function [ret] = acar(mat, acar_range, band)
% Applies the Acar filtering method to the data provided. Any
% entry which is outside of the acceptable bounds will be
% replaced with a NaN. The first acar_range number of indexes will also
% be NaN, as you cannot extrapolate data that isn't there
% Inputs:
%   mat: A [n-by-1] vector which contains the data you want to 
%       process.
%   acar_range: [int] the number of entries
%
%   band: [2 int vector] The range [start, end] of values you wish
%       to use the Acar filter on. If false, then analyze the full
%       range
%
% Returns:
%   ret: [band by 1] vector with filtered data
    
    if band
        r_1 = band(1);
        r_2 = band(2);
        
        if r_1 < acar_range
            disp("starting point less than acar range, changing starting point to (acar range)+1");
            r_1 = acar_range+1;
        end
    else
        [r_2, c] = size(mat);
        r_1 = 1;
    end

    ret = nan(1+r_2-r_1,1);
    k = r_1-1;
    for i = (acar_range+1):length(ret)
        a = mean(mat(k+i-acar_range:k+i-1,1));
        
        if abs(a-mat(k+i,1))> (0.2*a)
            ret(i,1) = NaN;
        elseif isnan(a-mat(k+i,1))
            ret(i,1) = NaN;
        else
            ret(i,1) = mat(k+i);
        end
    end
    %ret(1:acar_range,1) = NaN;
    
end

