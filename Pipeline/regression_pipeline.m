function [Data] = regression_pipeline(hr_files, ecg_files, aff_files, target, varargin)
% This function is the default pipeline for regression analysis of
% collected affect data.

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
%       of affect C and D". Logistic regression will be performed, with the
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
%       logistic regression. Default is {{'nothing'}}.
% 
%   duration: [bool] Whether to train and test a logsitic regression model
%       using the datapoints from the duration of the target affects as
%       category 1 and the control datapoints as category 0.
%       Results will be stored in Data.Reg.Off* fields. Default is true.
% 
%   onset: [bool] Whether to train and test a logistic regression using the
%       datapoints from the onset of instances of the target affects. Onset
%       is defined by on_band. Results will be stored in Data.Reg.On*
%       fields. Default is true.
% 
%   offset: [bool] Whether to train and test a logistic regression using the
%       datapoints from the offset of instances of the target affects as
%       category 1 and the control datapoints as category 0.
%       Offset is defined by off_band. Results will be stored in 
%       Data.Reg.Off* fields. Default is false.
% 
%   onoff: [bool] Whether to train and test a logistic regression using the
%       datapoints from both the onset and offset of instances of the target 
%       affects as category 1 and the control datapoints as category 0.
%       Onset is defined by on_band and offset is defind as off_band.
%       Results will be stored in Data.Reg.OO* fields. Default is false.
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
%   iterations: [int] Number of times logistic regressions are run on
%       randomly-sampled data from both the target (duration/onset/offset/onoff
%       depending on selection(s)) and control. The results of each
%       iteration are stored and plotted as histograms to give a better
%       picture of the model performance. Default is 100.
% 
%   verbose: [bool] Whether to print additional information to the command
%       window and generate additional figures for additional updates and
%       information as the pipeline runs. Default is false.
%
% Returns:
%   Data: [struct] Structure containing the loaded data and the results of
%       the regresion analysis, stored in the Data.Reg field as:
%
%       *_Reg: [x-by-y-by-z] matrix containing regression coefficients along
%       with the SE, tStats, and pValue for each coefficient, with the
%       first dimension being the iteration that the results are from.
%
%       *_PPhat: [x-by-3-by-2 matrix] a matrix containing the accuracy results of 
%       the trained logistic regression when classifying the test data. The
%       first dimension x = the number of iterations. The second dimension
%       is the [worst, pred, best], where worst = worst case scenario of
%       the 95% confidence interval, best = best case scenario of the 95%
%       confidence interval, and pred is the unadjusted predictive accuracy
%       of the logsitic regression model. The third dimension is size 2 for
%       the two categories, 1 = non-target, 2 = target.
%
%       *_AUC: [x-by-1 matrix] a vector containing the AUC values for the
%       logistic regression model on the testing data.


    % Input parser to read in optional parameters
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
    addParameter(p,'iterations', 100,@isscalar);
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
    % regression pipeline
    if p.Results.ecg_features
        Data.ECG.Aligned_Metrics = {};
    end
    
    Data.Reg = {};
    
    
    % Create key for improved figure plotting
    key = [{'RR'}];
    
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
      
    % If they selected regression analysis using timepoints during the
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
                new_dat{q} = Data.HR.PP{i};%feature_generation(Data.HR.PP{i}, bin, false);

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
        
        
        i = 1;
        % Run the logistic regression X times in order to get an
        % understanding of the performance across many different random
        % samples
        
        while i <= p.Results.iterations
            [Data.Reg.D_RegPP(i,:,:), Data.Reg.D_PPhat(i,:,:), Data.Reg.D_AUC(i,1)] = gen_regression([cont_mat(:,3:end);dur_mat(:,3:end)],...
                [zeros(size(cont_mat,1),1); ones(size(dur_mat,1),1)],...
                'log', 'verbose', false);
            i = i+1;
        end
        % Plot all of the relevant figures
        results_plotting(Data.Reg.D_RegPP, Data.Reg.D_PPhat, Data.Reg.D_AUC,key,'Duration');
    end
    
    
    if p.Results.onset || p.Results.offset || p.Results.onoff
        % If they selected regression analysis using timepoints during the
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
        
        
        % Run the logistic regression X times in order to get an
        % understanding of the performance across many different random
        % samples
        
        if p.Results.onoff % 'onoff' was set to true, so combine onset and offset
            for i = 1:p.Results.iterations
                [Data.Reg.OO_RegPP(i,:,:), Data.Reg.OO_PPhat(i,:,:), Data.Reg.OO_AUC(i,1)] = gen_regression([cont_mat;onset_mat;offset_mat],...
                    [zeros(size(cont_mat,1),1); ones(size(onset_mat,1),1); ones(size(offset_mat,1),1)],...
                    'log', 'verbose', false);
            end
            % Plot all of the relevant figures
            results_plotting(Data.Reg.OO_RegPP, Data.Reg.OO_PPhat, Data.Reg.OO_AUC,key,'Onset+Offset');
        end
        
        if p.Results.onset % 'onset' set to true
            i = 1;
            for i = 1:p.Results.iterations
                [Data.Reg.On_RegPP(i,:,:), Data.Reg.On_PPhat(i,:,:), Data.Reg.On_AUC(i,1)] = gen_regression([cont_mat;onset_mat],[zeros(size(cont_mat,1),1); ones(size(onset_mat,1),1)],...
                    'log', 'verbose', false);
            end
            % Plot all of the relevant figures
            results_plotting(Data.Reg.On_RegPP, Data.Reg.On_PPhat, Data.Reg.On_AUC,key,'Onset');
        end
        
        if p.Results.offset % 'offset' set to true
            for i = 1:p.Results.iterations
                [Data.Reg.Off_RegPP(i,:,:), Data.Reg.Off_PPhat(i,:,:), Data.Reg.Off_AUC(i,1)] = gen_regression([cont_mat;offset_mat],[zeros(size(cont_mat,1),1); ones(size(offset_mat,1),1)],...
                    'log', 'verbose', false);
            end
            % Plot all of the relevant figures
            results_plotting(Data.Reg.Off_RegPP, Data.Reg.Off_PPhat, Data.Reg.Off_AUC,key,'Offset');
        end
    end
end

function [] = results_plotting(reg, phat, auc, key, source)
% Function to store all of the plotting code
%
% Inputs:
%   reg: [x-by-y-by-z matrix] *_Reg field from Data
%
%   phat: [x-by-3-by-2 matrix] *_PPhat field from Data
%
%   auc: [x-by-1 matrix] *_AUC from Data
%
%   key: [1-by-n cell array] cell array containing the names of each of the
%       features used during the logistic regression. These names are used
%       to correctly label plots.
%
%   source: [string] What label to add to the plots to denote their source
        

    % Plot the histograms for the beta values and the confidence intervals
    row = floor(sqrt(size(reg,2)));
    col = ceil(size(reg,2)/row);
    figure;
    for i =1:size(reg,2)
        subplot(row,col,i)
        histogram(reg(:,i,5));
        hold on;
        histogram(reg(:,i,6));
        histogram(reg(:,i,1));
        legend('lower bound', 'upper bound', 'estimation');
        if i == 1
            title(strcat('estimated Beta coefficient (',source,')'));
        else
            title(strcat('estimated : ', key{i-1},' coefficient'));
        end
    end

    % Plot p-values for the coefficents
    figure;
    for i = 1:size(reg,2)
        subplot(row,col,i)
        histogram(reg(:,i,4));
        if i == 1
            title({'p-value for Beta', strcat('coefficient (',source,')')});
        else
            title({'p-value for : ', strcat(key{i-1},' coefficient')});
        end
    end


    % Plot the accuracy results
    figure;
    titles = {{'Non-Target Accuracy',strcat('(',source,')')},{'Target Accuracy',strcat('(',source,')')}};
    for i = 1:2
        subplot(1,2,i)
        histogram(phat(:,1,i));
        hold on;
        histogram(phat(:,2,i));
        histogram(phat(:,3,i));
        hold off;
        title(titles{i});
        legend('Lower Bound', 'Pred', 'Upper Bound');
    end

    % Plot spread of AUC
    figure;
    histogram(auc);
    title(strcat('AUC Values (',source,')'));
end

function [ret,locs] = ecg_preprocess(mat, amp, cut_bin)
% Function to take the ECG data preprocesses it by removing ECG values
% which fall outside of the accepted amplitude
% Inputs:
%   mat: [n-by-1] vector containing ecg data
%
%   amp: [int], minimum amplitude of accepted ECG values. If a value
%   falls outside of this bound, all values [cut_bin] before and after
%   it are replaced with NaNs.
%
%   cut_bin: [int], the amount of indexes before and after entries
%   which fail [amp] that are replaced with NaNs.
    
% Returns:
%   ret: [n-by-1] matrix containing the ecg data with all removed
%       values replaced by NaNs
%
%   locs: [m-by-1] index of all values which fall outside of the bounds
%       described by 'amp'

    locs = find(abs(mat)>amp);
    ret = mat;
    max_len = length(ret);
    
    for i = 1:length(locs)
    
        if locs(i) <= cut_bin
            ret(1:locs(i)+cut_bin) = NaN;
        elseif locs(i)+cut_bin >= max_len
            ret(locs(i)-cut_bin:end) = NaN;
        else
            ret(locs(i)-cut_bin:locs(i)+cut_bin) = NaN;
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