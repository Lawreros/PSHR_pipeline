function [ret] = bandpass(mat, l_band, u_band, rang)
    % Applies a bandpass filtering to vector provided. Any RR-interval
    % value outside of the range specified by the lower and upper bounds is
    % replaced with a NaN
    %   Inputs:
    %       Data: The Data structure
    %       source: [string], Which matrix from the structure you want to
    %       use
    %       l_band: [int], the lower bounding value for the bandpass filter
    %       u_band: [int], the upper bounding value for the bandpass filter
    %       rang: [2 int vector] The range [start, end] of values you want
    %       to use the bandpass on. If false, then analyze the whole
    %       range

    if rang
        r_1 = rang(1);
        r_2 = rang(2);
    else
        [r_2,c] = size(mat);
        r_1 = 1;
    end

    %Create copy of matrix to edit
    ret = mat;
    
    for i = r_1:r_2
        if (ret(i,1)>= u_band) || (ret(i,1) <= l_band)
            ret(i,3) = NaN;
        end
    end 
end