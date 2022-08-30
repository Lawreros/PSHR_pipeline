function [ret] = karlsson(mat, band)
% Applies the Karlsson filtering method to the data provided. Any
% entry which changes more or less than 20% of the mean of the entries
% directly before and after it will be replaced with a NaN
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
    ret(1,1) = NaN;
    
    for i = 2:(length(r_2)-1)
        a = (ret(i-1,1)+ret(i+1,1))/2;
        
        if abs(a-ret(i,1)) > (0.2*a)
            ret(i,1) = NaN;
        end
    end
end

