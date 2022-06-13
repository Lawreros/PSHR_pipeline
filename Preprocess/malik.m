function [ret] = malik(mat,band)
% Applies the malik filtering method to the provided data. If an
% entry increases or decreases by more than 20%, it is replaced with
% a NaN
%     Inputs:
%         mat: A [n-by-1] vector which contains the data you want to 
%         process.
%         band: [2 int vector] The range [start, end] of values you want
%         to calculate the pnnX of. If false, then analyze the whole range
% 
%     Returns:
%         ret: [band-by-1] vector of portion of mat defined by band

    if band
        r_1 = band(1);
        r_2 = band(2);
    else
        [r_2, c] = size(mat);
        r_1 = 1;
    end
        
    %Create copy of matrix to edit
    ret = mat(r_1:r_2,1);
    
    for i = 1:(length(ret)-1)
        if abs(ret(i,1) - ret(i+1,1)) > (0.2*ret(i,1))
            ret(i,1) = NaN;
        end
    end
end

