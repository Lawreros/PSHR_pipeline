function [idx] = clustering_pipeline(hr_files,aff_files,target,varargin)
% Pipeline for the creation of a random forest model
% Required Inputs:
%   hr_files: [1-by-n cell array] list containing the location of the hr
%       files you wish to load and analyze.
%
%   aff_files: [1-by-n cell array] list containging the locations of the
%       coded Affect files that you wish to load and analyze.
%
%   target: [1-by-t cell] cell array containing the target affects 
%       formatted as {{'A','B'},{'C','D'}}, which can be translated to 
%       "the union of affect A and affect B, intersecting with the union
%       of affect C and D". Clustering will be performed, with the
%       target data being category 1 and control data (i.e. non-target 
%       timepoints randomly selected) will be category 0. This is used to 
%       call Import/table_combo, so see there for more information.
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
%   omit: [1-by-o cell] cell array containing a list of affects to omit
%       from the target data points used for analysis. For
%       example, inputing {{'Q','Z'}} will results in all target
%       affect data points which occur at the same time as either
%       affect Q or affect Z to be labeled as category 0 instead of 
%       category 1. Default is {{'nothing'}}.
% 
%   plots: [bool] whether to plot a series of 3D scatterplots of the
%       different pairs of features. Default is false.

% Returns:
%   idx: [1-by-n matrix] Matrix containing the cluster id for each of the
%       datapoints.

    p = inputParser;
    addParameter(p,'bin', {[5,0], 'second'}, @iscell);
    addParameter(p,'omit', {{'nothing'}}, @iscell);
    addParameter(p,'plots',false, @islogical);
    
    parse(p,varargin{:});
    
    Data = pshr_load('HR', hr_files, 'Affect', aff_files, 'align', true, 'verbose', false);

    
    % Create key for improved figure plotting
    key = [{'RR'}];

    %% Generate additional features and call DBSCAN clustering on the data
    % Iterate through the files and preproces them using bandpass filtering
    % before calculating features (RMSSD, pNN50, etc.)
    for i = 1:length(hr_files)
        % Create new copy of Data.HR.Raw for PreProcessing
        Data.HR.PP{i} = Data.HR.Raw{i};

        % Apply basic bandpass to HR data
        Data.HR.PP{i}(:,3) = bandpass(Data.HR.PP{i}(:,3), 300, 1600, false);
        % Generate new features
        [Data.HR.PP{i}, k] = feature_generation(Data.HR.PP{i}, p.Results.bin, false, [1,1,1,1]);
        if i == 1
            key = [key,k];
        end
    end

    q=1;
    for i = 1:length(Data.HR.Affect)

        if ~isempty(Data.HR.Affect{i})
            % Run table_combo to get the start/stop times for both the
            % target affects and the affects to omit
            new_tabs{q} = table_combo(Data.HR.Affect{i}, target{:}, 'omit', p.Results.omit{:});

            % Create new_dat for data manipulation
            new_dat{q} = Data.HR.PP{i};
                
            new_dat{q} = affect_mark(new_dat{i}, new_tabs{i}(1,:), false);
                
            q=q+1;
        else
            % Skip any data without a corresponding affect file
            disp(strcat('No affect found for :', Data.HR.files{i}));

        end
    end
    
    % Combine all data into one big matrix
    big = vertcat(new_dat{:});

    % Run the clustering function
    [idx] = newFdbscan(big(:,3:end-1), key, big(:,end), 50, 10, p.Results.plots);

    % Calculate the amount of culsters that appeared and how much of the
    % data belongs to each cluster
    dats = unique(idx);
    for i = 2:length(dats)
        disp(strcat('Percentage of points belonging to cluster ',string(dats(i)),': ', string(sum(idx(:,1)==dats(i))*100/length(idx)),'%'))
    end
    
    disp(strcat('Percentage of unassigned datapoints: ', string(sum(idx(:,1)==-1)*100/length(idx)),'%'));
end