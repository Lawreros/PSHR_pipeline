function [idx] = clustering_pipeline(hr_files,aff_files,varargin)
% Pipeline for the creation of a random forest model
% Required Inputs:
%   hr_files: [1-by-n cell array] list containing the location of the hr
%       files you wish to load and analyze.
%
%   aff_files: [1-by-n cell array] list containging the locations of the
%       coded Affect files that you wish to load and analyze.
%
% Optional Parameters:
%   bin: [1-by-2 cell array] Used for creating the feature_gen vector
%       results from a sliding bin of Y seconds or entries. This takes the
%       format of {[before, after], 'units'}, so if you want to have a bin 
%       of 5 seconds before (including current RR-interval) and 3 seconds 
%       after: {[5,3], 'second'} or if you want the 5 entries before and 
%       3 entries after the index: {[5,3], 'measure'}. Default value is
%       {[5,0], 'second'}.
%
%   plots: [bool] whether to plot a series of 3D scatterplots of the
%       different pairs of features. Default is false.

% Returns:
%   idx: [1-by-n matrix] Matrix containing the cluster id for each of the
%       datapoints.

    p = inputParser;
    addParameter(p,'plots',false, @islogical);
    addParameter(p,'bin', {[5,0], 'second'}, @iscell);
    
    parse(p,varargin{:});
    
    
    % Load in data
     aff_list = {'SIB','ISB','inappropriate face related behavior','polar strap adjustment/removal'...
        'repetitive behaviors','inappropriate movement','crying', 'pulling at pants'};
    
    Data = pshr_load('HR', hr_files, 'Affect', aff_files, 'align', true, 'verbose', false);

    
    %% RR-interval preprocessing
    for i = 1:length(hr_files)
        Data.HR.PP{i} = Data.HR.Raw{i};
        Data.HR.PP{i} = affect_mark(Data.HR.PP{i}, Data.HR.Affect{i},aff_list); %mark the affect locations
        % We'll just work with bandpassing for now...
        Data.HR.PP{i}(:,3) = bandpass(Data.HR.PP{i}(:,3), 300, 1600, false);
    end
    

    %% Generate additional features and call DBSCAN clustering on the data
    
     for i=1:length(Data.HR.PP)
        Data.HR.PP{i} = feature_generation(Data.HR.PP{i}, p.Results.bin, false);
     end
        
     big = vertcat(Data.HR.PP{:});

     [idx] = newFdbscan(big(:,3:end-1), {'RR-interval','RMSSD','pNN50','SDNN','SDSD'}, big(:,end), 50, 10, p.Results.plots);
     
     dats = unique(idx);
     for i = 2:length(dats)
         disp(strcat('Percentage of points belonging to cluster ',string(dats(i)),': ', string(sum(idx(:,1)==dats(i))*100/length(idx)),'%'))
     end
     
     disp(strcat('Percentage of unassigned datapoints: ', string(sum(idx(:,1)==-1)*100/length(idx)),'%'));
end


function [mat] = feature_generation(mat, bin, band)
% Function for generating the different features for multiple recording
% sessions

% Inputs:
%   mat: [n-by-m matrix] where the third column is the data used for
%       feature generation
%   bin: [1-by-2 cell array] The bin type you want to use for the feature
%       calculation
%   band: [1-by-2 matrix] The start and end index you wish to analyze (set
%       this to false to use all available data)

    mat(:,5) = rmssd_calc(mat(:,3), bin, band);
    mat(:,6) = pnnx_calc(mat(:,3),50, bin, band);
    mat(:,7) = sdnn_calc(mat(:,3),bin,band);
    mat(:,8) = sdsd_calc(mat(:,3),bin,band);
    
    %move coding into last column
    mat(:,end+1) = mat(:,4);
    mat(:,4) = [];

end
