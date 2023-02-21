function [] = t_test_pipeline(hr_files, aff_files, target, varargin)
% Pipeline for running t-tests on different aspects of the PHYSMON
% collected data. Effectively any combination of analyses you could want
% will be referenced here.

% Required Inputs:
%   hr_files: [1-by-n cell array] cell array containing the HR files that
%       you wish to run the t-test on
%
%   aff_files: [1-by-n cell array] cell array containing the affect files
%       which you wish to use with the hr_files you provided (in the same
%       order).
%
%   target: [1-by-m cell array] cell array containing the target affects 
%       formatted as {{'A','B'},{'C'}}. This is used to call
%       Import/table_combo, so see there for more information.
%
%
% Optional Parameters:
%   bin: [1-by-2 cell array] Used for creating a vector of the
%       feature results from a sliding bin of Y seconds or entries. This takes the
%       format of {[before, after], 'units'}, so if you want to have a bin 
%       of 5 seconds before (including current RR-interval) and 3 seconds 
%       after: {[5,3], 'second'} or if you want the 5 entries before and 
%       3 entries after the index: {[5,3], 'measure'}
%       If you don't want this, set bin to false.
%
%   omit: [1-by-b cell array] cell array containing the affects to omit, 
%       formatted as {{'A','B'},{'C'}}. This is used to call
%       Import/table_combo, so see there for more information.
%
%   duration: [bool] Whether to run the t-test analysis using data
%       collected from the timepoints during the target affects against
%       comparison to contol data from timepoints outside of the target
%       affects. If false, then this analysis is not run.
%
%   onset: [bool] Whether to run the t-test analysis using data
%       collected from the timepoints during the target affects against
%       comparison to contol data from timepoints outside of the target
%       affects. If false, then this analysis is not run.
%
%   band: [1-by-2 vector] The number of entries before and after the onset
%       to include in the returned matrix. Default is [3,0] for the three entries
%       before the onset (this will results in 4 values being used, with 
%       the last value being the first timepoint of the affect).
%
% Returns:
%   

    p = inputParser;
    addParameter(p, 'bin', {[5,0], 'second'}, @iscell);
    addParameter(p, 'omit', {{'nothing'}}, @iscell);
    addParameter(p, 'duration', true, @islogical);
    addParameter(p, 'onset', true, @islogical);
    addParameter(p, 'band',[3,0], @ismatrix);
    
    parse(p,varargin{:});
    
    omit = p.Results.omit;
    bin = p.Results.bin;
    duration = p.Results.duration;
    onset = p.Results.onset;
    

    Data = pshr_load('HR', hr_files, 'Affect', aff_files, 'align', true, 'verbose', false);
    
     for i = 1:length(hr_files)
        Data.HR.PP{i} = Data.HR.Raw{i};
        % We'll just work with bandpassing for now...
        Data.HR.PP{i}(:,3) = bandpass(Data.HR.PP{i}(:,3), 300, 1600, false);
    end
    
    if duration
        %% Duration sample
        % Sample the HR data using the new affect tables, then run a t-test on the
        % differences between the 

        q=1;
        store_names = {};
        for i = 1:length(Data.HR.Affect)

            if ~isempty(Data.HR.Affect{i})
                new_tabs{q} = table_combo(Data.HR.Affect{i}, target{:}, 'omit', omit{:});

                % Run the feature_generation function on the HR data
                [new_dat{q}, keys] = feature_generation(Data.HR.PP{i}, bin, false, [1,1,0,0]);

                store_names{q} = Data.HR.files{i};
                q=q+1;
            else
                disp(strcat('No affect found for :', Data.HR.files{i}));

            end
        end
        
        %Add to keys
        keys = [{'Timestamp'},{'BPM'},{'RR-interval'}, keys];


        % Select data using starts and stops
        dur_mat = [];
        cont_mat = [];
%         a = 10; % Minimum duration to count as a space outside of both target and omission
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
                idx = find(vec > 2);% a);
                
                starts =[inst(1)]; % Account for the fact that the first timepoint may be part of control
                ends =[];

                for k=1:length(idx)
                    starts = [starts,inst(idx(k)+1)];
                    ends = [ends,inst(idx(k))];
                end
                ends = [ends, inst(end)]; %grab the last ending

                % Add buffer around all of the non-control starts/stops in order to
                % select control timepoints more than -sum(bin{1}) indices away from affects
                [starts, ends] = dilate(starts, ends, -sum(bin{1}), sum(bin{1}));

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
            for q = 1:size(dur_mat,2)
                
                % Calculate the parameters for Cohen's d
                n_1 = length(dur_mat(:,q));
                n_2 = length(cont_mat(:,q));
             
                mu_1 = mean(dur_mat(:,q));
                mu_2 = mean(cont_mat(:,q));
             
                sd_1 = std(dur_mat(:,q));
                sd_2 = std(cont_mat(:,q));
             
                pooled = sqrt(((n_1-1)*sd_1^2 + (n_2-1)*sd_2^2)/(n_1+n_2));
                
                disp("------------------------------------------");
                disp(strcat("Cohen's d for feature : ",keys{q}," = ", string((mu_1-mu_2)/pooled)));
                
                % Run ttest with the assumption of unequal variance
                [h, p_, ci, stats] = ttest2(dur_mat(:,q),cont_mat(:,q), 'vartype', 'unequal');
                t_stats(1,q) = p_;
                disp(strcat('two-sample t-test result for feature :', keys{q},' = ', string(p_)));
                disp(strcat('target mean: ', string(nanmean(dur_mat(:,q))),' control mean: ', string(nanmean(cont_mat(:,q)))));
                disp('95% CI estimated difference:');
                disp(ci);
                disp('Estimated standard deviation');
                disp(stats.sd);
            end
        end
        
        % ADD OPTION FOR CORRELATION/COVARIANCE FOR DIFFERENT FEATURES
        
%          % Check covariance
%         disp('covariance for problematic:');
%         disp(cov([train.cat_1;test.cat_1;unused.cat_1],1));
%         disp('correlation for problematic:');
%         disp(corrcoef([train.cat_1;test.cat_1;unused.cat_1]));
%         disp('covariance for non-problematic:');
%         disp(cov([train.cat_0;test.cat_0;unused.cat_0],1));
%         disp('correlation for non-problematic:');
%         disp(corrcoef([train.cat_0;test.cat_0;unused.cat_0]));
%         
    end

    if onset
        % Onset Sample
        % Pass the HR data and new affect tables into onset_sample in order to get
        % the samples of the onsets and controls

        q=1;
        store_names = {};
        new_tabs = {};
        
        % Go through each HR file with an associated affect file and create
        % a new start/stop table for the target affects, omitting the
        % "omit" affects
        for i = 1:length(Data.HR.Affect)

            if ~isempty(Data.HR.Affect{i})
                keep_table = table_combo(Data.HR.Affect{i}, target{:});
                
                if ~isempty(omit)
                    omit_table = table_combo(Data.HR.Affect{i}, omit{:});
                    new_tabs{q} = [keep_table(end-1,:); omit_table(end-1,:)];
                else %If they aren't omitting anything, then just use keep_table
                    new_tabs{q} = keep_table;
                end

                % Run the feature_generation function on the HR data
                [new_dat{q}, keys] = feature_generation(Data.HR.PP{i}, bin, false, [1,1,0,0]);

                store_names{q} = Data.HR.files{i};
                q=q+1;
            else
                disp(strcat('No affect found for :', Data.HR.files{i}));

            end
        end

        %Add to keys
        keys = [{'Timestamp'},{'BPM'},{'RR-interval'}, keys];


        onset_mat = [];
        offset_mat = [];
        cont_mat = [];
        
        for i = 1:length(new_tabs)
            
            if ~isempty(new_tabs{i}{1,2})
                disp(strcat('Running data from :',store_names{i}))
                
                % Sample [onset, offset, and control] timpoints using the
                % provided start/stop times from the target/omit affects
                [on_dump, off_dump, cont_dump] = onset_sample(new_dat{i}, new_tabs{i}(1,:), new_tabs{i}(2,:), 'band', p.Results.band, 'dilate', 6);
                
                % Create three big matrices to dump these sampled data
                % points
                onset_mat = [onset_mat; on_dump];
                offset_mat = [offset_mat; off_dump];
                cont_mat = [cont_mat; cont_dump];
            end

        end


        % Run a two sample t-test on the two groups for results, along with
        % confidence interval
        t_stats=zeros(1,size(onset_mat,2));
        
        onset_mat(any(isnan(onset_mat), 2), :) = [];
        cont_mat(any(isnan(cont_mat), 2), :) = [];

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
                disp(strcat("Cohen's d for feature :", keys{q}, '=', string((mu_1-mu_2)/pooled)));
                
                % Run ttest with the assumption of unequal variance
                [h, p_, ci, stats] = ttest2(onset_mat(:,q),cont_mat(:,q), 'vartype', 'unequal');
                t_stats(1,q) = p_;
                
                disp(strcat('two-sample t-test result for feature : ', keys{q},' = ', string(p_)));
                disp(strcat('target mean : ', string(nanmean(onset_mat(:,q))),' control mean : ', string(nanmean(cont_mat(:,q)))));
                disp('95% CI estimated difference:');
                disp(ci);
                disp('Estimated standard deviation');
                disp(stats.sd);
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
%   amnt: [int] The amount to dialte the start and end times. This value
%       can either be positive or negative.
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

    d_starts = max(1,starts-amnt);
    d_ends = ends+amnt;
    
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