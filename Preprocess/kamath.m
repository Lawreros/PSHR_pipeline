function [ret] = kamath(mat, band)
% Applies the Kamath filtering method to the data provided. Any
% entry which falls outside the bounds will be replaced with a NaN
%   Inputs:
%       mat: A [n-by-1] vector which contains the data you want to 
%       process.
%       band: [2 int vector] The range [start, end] of values you want
%       to use the Karlsson filter on. If false, then analyse the whole
%       range
%
%   Returns:
%       ret: [band by 1] vector with filtered data
   
    if band
        r_1 = band(1);
        r_2 = band(2);
    else
        [r_2, c] = size(mat);
        r_1 = 1;
    end
    
    %Create copy of matrix to edit
    ret = mat(r_1:r_2,1);
    k = r_1-1;
    for i = 1:(length(ret)-1)
        a = 0.325 * mat(k+i,1);
        b = 0.245 * mat(k+i,1);
        
        c = mat(k+i+1,1) - mat(k+i,1);
        d = mat(k+i,1) - mat(k+i+1,1);
        
        if (0 <= c) && (c <= a)
            %ret(i,1) = NaN;
        elseif (0 <= d) && (d <= b)
            %ret(i,1) = NaN;
        else
            ret(i,1) = NaN;
        end
    end

end

