function [mat, key] = feature_generation(mat, bin, band, meths)
% Function for generating the different features for multiple recording
% sessions

% Inputs:
%   mat: [n-by-m matrix] matrix of values you wish to generate additional
%       features from.
%       NOTE: the function currently assumes that the third
%       column is the column you wish to use when generating additional
%       features.
%
%   bin: [1-by-2 cell array] Used for creating the feature vector
%       results from a sliding bin of Y seconds or entries. This takes the
%       format of {[before, after], 'units'}, so if you want to have a bin 
%       of 5 seconds before (including current RR-interval) and 3 seconds 
%       after: {[5,3], 'second'} or if you want the 5 entries before and 
%       3 entries after the index: {[5,3], 'measure'}
%       If you don't want this, set bin to false.
%
%   band: [1-by-2 matrix] The range [start, end] of values you want
%       to calculate the feature of. If false, then analyze the whole
%       range of mat
%
%   meths: [1-by-4 matrix] Vector of binary values determining which
%       methods you wish to use for feature generation: [RMSSD, pNN50,
%       SDNN, SDSD]. For example, [1,0,0,1] means you want RMSSD and SDSD
%       features calculated.
%
% Returns:
%   mat: [n-by-q matrix] original matrix mat with additional feature
%       columns appended to the end
%
%   key: [1-by-k cell array] list of what features were generated for easy
%       figure labeling

    key{1} = 'n/a';

    if meths(1)
        mat(:,end+1) = rmssd_calc(mat(:,3), bin, band);
        key{end+1} = 'RMSSD';
    end
    if meths(2)
        mat(:,end+1) = pnnx_calc(mat(:,3),50, bin, band);
        key{end+1} = 'pnn50';
    end
    if meths(3)
        mat(:,6) = sdnn_calc(mat(:,3),bin,band);
        key{end+1} = 'SDNN';
    end
    if meths(4)
        mat(:,7) = sdsd_calc(mat(:,3),bin,band);
        key{end+1} = 'SDSD';
    end
    
    % remove 'n/a' key, because matlab is stupid
    key(1)=[];

end