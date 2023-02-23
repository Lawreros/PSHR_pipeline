function [] = t_test_pipeline(hr_files, ecg_files, aff_files, target, varargin)
% Pipeline for running t-tests on different aspects of the PHYSMON
% collected data.

% Required Inputs:
%   hr_files: [1-by-n cell array] List of all the files containing
%       RR-intervals for loading using pshr_load
%
%   ecg_files: [1-by-n cell array] List of all the files containing ECG
%       data for loading using pshr_load
%
%   aff_files: [1-by-n cell array] List of all affect files containing
%       coding for the RR-interval/ECG files
%
%   target: [1-by-t cell] cell array containing the target affects 
%       formatted as {{'A','B'},{'C','D'}}, which can be translated to 
%       "the union of affect A and affect B, intersecting with the union
%       of affect C and D". The t-test will be performed, with the
%       target data being category 1 and control data (i.e. non-target 
%       timepoints randomly selected) will be category 0. This is used to 
%       call Import/table_combo, so see there for more information.
%
% Optional Parameters:
%   ecg_features: [bool] Whether to calculate additional features using the
%       ecg data and the relative distance between the different components
%       of the QRS complex. Default is false.
%
%   bin: [1-by-2 cell array] Used for creating a vector of the feature
%       results from a sliding bin of Y seconds or entries. This takes the
%       format of {[before, after], 'units'}, so if you want to have a bin 
%       of 5 seconds before (including current RR-interval) and 3 seconds 
%       after: {[5,3], 'second'} or if you want the 5 entries before and 
%       3 entries after the index: {[5,3], 'measure'}
%       If you don't want this, set bin to false.
%
%   omit: [1-by-o cell] cell array containing a list of affects to omit
%       from both the target and control data points used for analysis. For
%       example, inputing {{'Q','Z'}} will results in all target
%       affect/control data points which occur at the same time as either
%       affect Q or affect Z to be excluded from the training of the
%       t-test. Default is {{'nothing'}}.
% 
%   duration: [bool] Whether to run the t-test using the datapoints from 
%       the duration of the target affects as category 1 and the control
%       datapoints as category 0.
%       Default is true.
% 
%   onset: [bool] Whether to run a t-test using the
%       datapoints from the onset of instances of the target affects. Onset
%       is defined by on_band. Default is true.
% 
%   offset: [bool] Whether to run a t-test using the
%       datapoints from the offset of instances of the target affects as
%       category 1 and the control datapoints as category 0.
%       Offset is defined by off_band. Default is false.
% 
%   onoff: [bool] Whether to run a t-test using the
%       datapoints from both the onset and offset of instances of the target 
%       affects as category 1 and the control datapoints as category 0.
%       Onset is defined by on_band and offset is defind as off_band.
%       Default is false.
% 
%   on_band: [1-by-2 matrix] The bounds, relative to the first timepoint of
%       a target affect instance, which define what counts as an "onset". For
%       example, inputting [3,2] would define onset as the first datapoint
%       of an instance of the target affect, the three datapoints
%       immediatly before the first datapoint, and the two datapoints after
%       the first datapoint. Default is [3,0].
% 
%   off_band: [1-by-2 matrix] The bounds, relative to the first timepoint of
%       a target affect instance, which define what counts as an "offset". For
%       example, inputting [3,2] would define offset as the last datapoint
%       of an instance of the target affect, the three datapoints
%       immediatly before the last datapoint, and the two datapoints after
%       the last datapoint. Default is [0,3].
% 
%   verbose: [bool] Whether to print additional information to the command
%       window and generate additional figures for additional updates and
%       information as the pipeline runs. Default is false.
%
% Returns:

    p = inputParser;
    addParameter(p, 'ecg_features', false, @islogical);
    addParameter(p, 'bin', {[5,0],'second'}, @iscell);
    addParameter(p, 'omit', {{'nothing'}}, @iscell);
    addParameter(p, 'duration', true, @islogical);
    addParameter(p, 'onset', true, @islogical);
    addParameter(p, 'offset', false, @islogical);
    addParameter(p, 'onoff', false, @islogical);
    addParameter(p, 'on_band', [3,0], @ismatrix);
    addParameter(p, 'off_band', [0,3], @ismatrix);
    addParameter(p, 'verbose', false, @islogical);
    parse(p, varargin{:});


    % Change what is loaded depending on whether HR and/or ECG data is
    % present
    if iscell(ecg_files)
        Data = pshr_load('HR', hr_files, 'ECG', ecg_files, 'Affect', aff_files, 'align', true, 'verbose', p.Results.verbose);
    else
        Data = pshr_load('HR', hr_files, 'Affect', aff_files, 'align', true, 'verbose', p.Results.verbose);
    end
    
    % Create Aligned_Metrics field in Data structure to store results of
    % ecg alignment
    if p.Results.ecg_features
        Data.ECG.Aligned_Metrics = {};
    end
    
    % Create key for improved figure plotting
    key = [{'Time'},{'BPM'},{'RR'}];
    
    % Iterate through the files and preproces them using bandpass filtering
    % before calculating features (RMSSD, pNN50, etc.)
    for i = 1:length(hr_files)
        % Create new copy of Data.HR.Raw for PreProcessing
        Data.HR.PP{i} = Data.HR.Raw{i};
        
        % Apply basic bandpass to HR data
        Data.HR.PP{i}(:,3) = bandpass(Data.HR.PP{i}(:,3), 300, 1600, false);
        [Data.HR.PP{i}, k] = feature_generation(Data.HR.PP{i}, p.Results.bin, false, [1,1,0,0]);
        if i == 1
            key = [key,k];
        end        
        
        % If they specified ecg_features, then append them to Data.HR.PP
        % using the 
        if p.Results.ecg_features
            [ecg_aligned, Data.ECG.Aligned_Metrics{i}] = ecg_rr_align(Data.HR.Raw{i}(:,[1,3]), Data.ECG.Raw{i}(:,[1,3]), 130, 'verbose', p.Results.verbose);
            Data.HR.PP{i}(:,end+1:end+3) = ecg_aligned(:,3:5); %Only look at [Q,R,S] complex
            if i == 1
                key = [key, 'Q-R','R','R-S']; % Add to key list for plotting
            end
        end
        
    end 
      
    % If they selected t-test using timepoints during the
    % target affects
    if p.Results.duration
        q=1;
        store_names = {};
        for i = 1:length(Data.HR.Affect)

            if ~isempty(Data.HR.Affect{i})
                % Run table_combo to get the start/stop times for both the
                % target affects and the affects to omit
                new_tabs{q} = table_combo(Data.HR.Affect{i}, target{:}, 'omit', p.Results.omit{:});

                % Create new_dat for data manipulation
                new_dat{q} = Data.HR.PP{i};

                store_names{q} = Data.HR.files{i};
                q=q+1;
            else
                % Skip any data without a corresponding affect file
                disp(strcat('No affect found for :', Data.HR.files{i}));

            end
        end    
        
        % Select data using starts and stops
        dur_mat = [];
        cont_mat = [];
        a = sum(p.Results.bin{1}); % Minimum duration to count as a space outside of both target and omission
        for i = 1:length(new_tabs)
            dur_dump = [];
            cont_dump = [];
            if ~isempty(new_tabs{i}{1,2})
                disp(strcat('Running data from :',store_names{i}))

                % Select data from outside both the omission and target datasets

                aff_vec = affect_mark(zeros(size(new_dat{i},1),1), new_tabs{i}(1,:), false);
                aff_vec(:,1)=[];
                om_vec = affect_mark(zeros(size(new_dat{i},1),1), new_tabs{i}(2,:), false);
                om_vec(:,1)=[];

                % Find times both omit and target affect(s) are not present
                % to use as control. Get start and stop times for these
                % control periods
                inst = find((aff_vec + om_vec)==0);
                vec = inst(2:end)-inst(1:end-1);
                
                % If there is a break in the control period greater than
                % 2 indices, then mark it as a start/stop time
                idx = find(vec > 2);
                starts =[inst(1)]; % Account for the fact that the first timepoint may be part of control
                ends =[];

                for k=1:length(idx)
                    starts = [starts,inst(idx(k)+1)];
                    ends = [ends,inst(idx(k))];
                end
                ends = [ends, inst(end)]; %grab the last ending

                % Add buffer around all of the non-control starts/stops in order to
                % select control timepoints more than X away from affects
                [starts, ends] = dilate(starts, ends, [-a, -a], a);

                 % Select all target data
                for j = 1:length(new_tabs{i}{1,2})
                    dur_dump = [dur_dump; new_dat{i}(new_tabs{i}{1,2}(j):new_tabs{i}{1,3}(j),:)];
                end

                % Select all control data
                for j = 1:length(starts)
                    cont_dump = [cont_dump;new_dat{i}(starts(j):min(ends(j),size(new_dat{i},1)),:)];
                end

                % Create two big matrices for the data from the duration of
                % the target affects and the control timepoints across all
                % provided files
                dur_mat = [dur_mat; dur_dump];
                cont_mat = [cont_mat; cont_dump];
            end

        end
        
        % Run a two sample t-test on the two groups for results, along with
        % confidence interval
        t_stats=zeros(1,size(dur_mat,2));
        
        dur_mat(any(isnan(dur_mat), 2), :) = [];
        cont_mat(any(isnan(cont_mat), 2), :) = [];

        if isempty(dur_mat)
            disp('No instances of target event(s) occuring. This may be due to misspelling of the target(s).');
        else
            disp("------------------------------------------");
            disp("DURATION ANALYSIS RESULTS");
            for q = 2:size(dur_mat,2)
                
                % Calculate the parameters for Cohen's d
                n_1 = length(dur_mat(:,q));
                n_2 = length(cont_mat(:,q));
             
                mu_1 = mean(dur_mat(:,q));
                mu_2 = mean(cont_mat(:,q));
             
                sd_1 = std(dur_mat(:,q));
                sd_2 = std(cont_mat(:,q));
             
                pooled = sqrt(((n_1-1)*sd_1^2 + (n_2-1)*sd_2^2)/(n_1+n_2));
                
                disp("------------------------------------------");
                disp(strcat("Cohen's d for feature : ",key{q}," = ", string((mu_1-mu_2)/pooled)));
                
                % Run ttest with the assumption of unequal variance
                [h, p_, ci, stats] = ttest2(dur_mat(:,q),cont_mat(:,q), 'vartype', 'unequal');
                t_stats(1,q) = p_;
                disp(strcat('two-sample t-test result for feature :', key{q},' = ', string(p_)));
                disp(strcat('target mean: ', string(nanmean(dur_mat(:,q))),' control mean: ', string(nanmean(cont_mat(:,q)))));
                disp('95% CI estimated difference:');
                disp(ci);
                disp('Estimated standard deviation');
                disp(stats.sd);
            end
        end
    end
    
    
    if p.Results.onset || p.Results.offset || p.Results.onoff
        % If they selected t-test using timepoints during the
        % onset of target affects
        q=1;
        store_names = {};
        new_tabs = {};
        
        % Go through each HR file with an associated affect file and create
        % a new start/stop table for the target affects, omitting the
        % "omit" affects
        for i = 1:length(Data.HR.Affect)

            if ~isempty(Data.HR.Affect{i})
                keep_table = table_combo(Data.HR.Affect{i},target{:});
                
                if ~isempty(p.Results.omit)
                    omit_table = table_combo(Data.HR.Affect{i}, p.Results.omit{:});
                    new_tabs{q} = [keep_table(end-1,:); omit_table(end-1,:)];
                else %If they aren't omitting anything, then just use keep_table
                    new_tabs{q} = keep_table;
                end

                % Run the feature_generation function on the HR data
                new_dat{q} = Data.HR.PP{i};%feature_generation(Data.HR.PP{i}, bin, false);

                store_names{q} = Data.HR.files{i};
                q=q+1;
            else
                disp(strcat('No affect found for :', Data.HR.files{i}));

            end
        end


        onset_mat = [];
        offset_mat = [];
        cont_mat = [];

        for i = 1:length(new_tabs)
            
            if ~isempty(new_tabs{i}{1,2})
                disp(strcat('Running data from :',store_names{i}))
                
                % Sample [onset, offset, and control] timpoints using the
                % provided start/stop times from the target/omit affects
                [on_dump, off_dump, cont_dump] = onset_sample(new_dat{i}(:,3:end), new_tabs{i}(1,:), new_tabs{i}(2,:),...
                    'on_band', p.Results.on_band, 'off_band', p.Results.off_band, 'dilate', 10, 'omit_nan', true);
                
                % Create three big matrices to dump these sampled data
                % points
                onset_mat = [onset_mat; on_dump];
                offset_mat = [offset_mat; off_dump];
                cont_mat = [cont_mat; cont_dump];
            end

        end
        
        
        
                
        if p.Results.onoff
            % Run a two sample t-test on the two groups for results, along with
            % confidence interval
            t_stats=zeros(1,size(onset_mat,2));

            if isempty(onset_mat)
                disp('No instances of target event(s) occuring. This may be due to misspelling of the target(s).');
            else
                disp("------------------------------------------");
                disp("ONSET+OFFSET ANALYSIS RESULTS");
                for q = 1:size(onset_mat,2)

                    % Calculate the parameters for Cohen's d
                    n_1 = length([onset_mat(:,q);offset_mat(:,q)]);
                    n_2 = length(cont_mat(:,q));

                    mu_1 = mean([onset_mat(:,q);offset_mat(:,q)]);
                    mu_2 = mean(cont_mat(:,q));

                    sd_1 = std([onset_mat(:,q);offset_mat(:,q)]);
                    sd_2 = std(cont_mat(:,q));

                    pooled = sqrt(((n_1-1)*sd_1^2 + (n_2-1)*sd_2^2)/(n_1+n_2));

                    disp("------------------------------------------");
                    disp(strcat("Cohen's d for feature :", key{q}, '=', string((mu_1-mu_2)/pooled)));

                    % Run ttest with the assumption of unequal variance
                    [h, p_, ci, stats] = ttest2([onset_mat(:,q);offset_mat(:,q)],cont_mat(:,q), 'vartype', 'unequal');
                    t_stats(1,q) = p_;

                    disp(strcat('two-sample t-test result for feature : ', key{q},' = ', string(p_)));
                    disp(strcat('target mean : ', string(nanmean([onset_mat(:,q);offset_mat(:,q)])),' control mean : ', string(nanmean(cont_mat(:,q)))));
                    disp('95% CI estimated difference:');
                    disp(ci);
                    disp('Estimated standard deviation');
                    disp(stats.sd);
                end
            end
        end
        
        if p.Results.onset
            
            % Run a two sample t-test on the two groups for results, along with
            % confidence interval
            t_stats=zeros(1,size(onset_mat,2));

            if isempty(onset_mat)
                disp('No instances of target event(s) occuring. This may be due to misspelling of the target(s).');
            else
                disp("------------------------------------------");
                disp("ONSET ANALYSIS RESULTS");
                for q = 1:size(onset_mat,2)

                    % Calculate the parameters for Cohen's d
                    n_1 = length(onset_mat(:,q));
                    n_2 = length(cont_mat(:,q));

                    mu_1 = mean(onset_mat(:,q));
                    mu_2 = mean(cont_mat(:,q));

                    sd_1 = std(onset_mat(:,q));
                    sd_2 = std(cont_mat(:,q));

                    pooled = sqrt(((n_1-1)*sd_1^2 + (n_2-1)*sd_2^2)/(n_1+n_2));

                    disp("------------------------------------------");
                    disp(strcat("Cohen's d for feature :", key{q}, '=', string((mu_1-mu_2)/pooled)));

                    % Run ttest with the assumption of unequal variance
                    [h, p_, ci, stats] = ttest2(onset_mat(:,q),cont_mat(:,q), 'vartype', 'unequal');
                    t_stats(1,q) = p_;

                    disp(strcat('two-sample t-test result for feature : ', key{q},' = ', string(p_)));
                    disp(strcat('target mean : ', string(nanmean(onset_mat(:,q))),' control mean : ', string(nanmean(cont_mat(:,q)))));
                    disp('95% CI estimated difference:');
                    disp(ci);
                    disp('Estimated standard deviation');
                    disp(stats.sd);
                end
            end
        end
        
        if p.Results.offset
            % Run a two sample t-test on the two groups for results, along with
            % confidence interval
            t_stats=zeros(1,size(offset_mat,2));

            if isempty(offset_mat)
                disp('No instances of target event(s) occuring. This may be due to misspelling of the target(s).');
            else
                disp("------------------------------------------");
                disp("OFFSET ANALYSIS RESULTS");
                for q = 1:size(offset_mat,2)

                    % Calculate the parameters for Cohen's d
                    n_1 = length(offset_mat(:,q));
                    n_2 = length(cont_mat(:,q));

                    mu_1 = mean(offset_mat(:,q));
                    mu_2 = mean(cont_mat(:,q));

                    sd_1 = std(offset_mat(:,q));
                    sd_2 = std(cont_mat(:,q));

                    pooled = sqrt(((n_1-1)*sd_1^2 + (n_2-1)*sd_2^2)/(n_1+n_2));

                    disp("------------------------------------------");
                    disp(strcat("Cohen's d for feature :", key{q}, '=', string((mu_1-mu_2)/pooled)));

                    % Run ttest with the assumption of unequal variance
                    [h, p_, ci, stats] = ttest2(offset_mat(:,q),cont_mat(:,q), 'vartype', 'unequal');
                    t_stats(1,q) = p_;

                    disp(strcat('two-sample t-test result for feature : ', key{q},' = ', string(p_)));
                    disp(strcat('target mean : ', string(nanmean(offset_mat(:,q))),' control mean : ', string(nanmean(cont_mat(:,q)))));
                    disp('95% CI estimated difference:');
                    disp(ci);
                    disp('Estimated standard deviation');
                    disp(stats.sd);
                end
            end
        end
    end
end


function [dstarts, dends] = dilate(starts, ends, amnt, cap)
% Function which dilates the non-zero values of vector vec by the amount
% specified by amnt. In other words, if you want to add 5 seconds before
% and after each affect instance (pair of start/stop timepoints), you call
% this function with amnt = 5
%
% Inputs:
%   starts: [1-by-n matrix] The start times for a given affect
%
%   ends: [1-by-n matrix] The end times for a given affect
%
%   amnt: [1-by-2 matrix] The amount to dialte the start and end times. This value
%       can either be positive or negative, with the first value being
%       subtracted from start times and the second being added to end
%       times.
%
%   cap: [int] The minimum amount of time between the start and stop of two
%       different instances of an affect for them to be considered seperate.
%       For example, if cap = 5, then if instance_1 and instance_2 of a given
%       affect are less than 5 seconds/timepoints apart they are lumped
%       together and considered one instance of the affect.
%
% Returns:
%   dstarts: [1-by-m matrix] The dilated start timepoints
%
%   dends: [1-by-m matrix] The dilated end timepoints

    di_vec = [];

    d_starts = max(1,starts-amnt(1));
    d_ends = ends+amnt(2);
    
    for i = 1:length(d_starts)
        di_vec(d_starts(i):d_ends(i),1)=1;
    end
    
    
    inst = find(di_vec==1);
    vec = inst(2:end)-inst(1:end-1);
    idx = find(vec > cap); % find starts greater in distance then the band being used so you don't
                         % accidently sample from another problematic
                         % behavior
    dstarts = [inst(1)];
    dends = [];
    for k=1:length(idx)
        dstarts = [dstarts,inst(idx(k)+1)];
        dends = [dends,inst(idx(k))];
    end
    dends =[dends, inst(end)];
    
end