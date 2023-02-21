function [mat, key] = feature_generation(mat, bin, band, meths)
% Function for generating the different features for multiple recording
% sessions

% Inputs:
%   mat: [n-by-m matrix]
%   bin: [1-by-2 cell array] The bin type you want to use
%   band: [1-by-2 matrix]
%   meths: [1-by-4 matrix] Vector of binary values determining which
%       methods you wish to use for feature generation.
%
% Returns:
%   mat:
%
%   key: [1-by-k cell array]

    key{1} = 'n/a';

    if meths(1)
        mat(:,end+1) = rmssd_calc(mat(:,3), bin, band);
        key{end+1} = {'RMSSD'};
    end
    if meths(2)
        mat(:,end+1) = pnnx_calc(mat(:,3),50, bin, band);
        key{end+1} = {'pnn50'};
    end
    if meths(3)
        mat(:,6) = sdnn_calc(mat(:,3),bin,band);
        key{end+1} = {'SDNN'};
    end
    if meths(4)
        mat(:,7) = sdsd_calc(mat(:,3),bin,band);
        key{end+1} = {'SDSD'};
    end
    
    % remove 'n/a' key, because matlab is stupid
    key(1)=[];

end