function [ret] = bandpass(mat, l_band, u_band, rang)
% Applies a bandpass filtering to vector provided. Any vector
% value outside of the range specified by the lower and upper bounds is
% replaced with a NaN
%   Inputs:
%       mat: [n-by-1 vector] vector of values which are being bandpass
%       filtered
%
%       l_band: [int], the lower bounding value for the bandpass filter
%
%       u_band: [int], the upper bounding value for the bandpass filter
%
%       rang: [2 int vector] The range [start, end] of values you want
%       to use the bandpass on. If false, then analyze the whole
%       range
    
%   Returns:
%       ret: [m-by-1 vector] vector of values that have been bandpass
%       filtered. This should be the length defined by rang

    if rang
        r_1 = rang(1);
        r_2 = rang(2);
    else
        [r_2,c] = size(mat);
        r_1 = 1;
    end

    %Create copy of matrix to edit
    ret = mat(r_1:r_2,1);
    
    for i = 1:length(ret)
        if (ret(i,1)>= u_band) || (ret(i,1) <= l_band)
            ret(i,1) = NaN;
        end
    end 
end