function [Data] = regression_pipeline(hr_files, ecg_files, ecg_features, aff_files, target, omit, bin, duration, verbose)
% This function is the default pipeline for regression analysis of
% collected affect data.

% Inputs:
%   hr_files: [1-by-n cell array] List of all the files containing
%       RR-intervals for loading using pshr_load
%
%   ecg_files: [1-by-n cell array] List of all the files containing ECG
%       data for loading using pshr_load
%
%   aff_files: [1-by-n cell array] List of all affect files containing
%       coding for the RR-interval/ECG files
%
% Returns:
%   Data: [struct] Structure containing the loaded data and the results of
%       the regresion analysis.


    aff_list = {'SIB','ISB','inappropriate face related behavior','polar strap adjustment/removal'...
        'repetitive behaviors','inappropriate movement','crying', 'pulling at pants', 'meta_chunk'};
    
    if iscell(ecg_files)
        Data = pshr_load('HR', hr_files, 'ECG', ecg_files, 'Affect', aff_files, 'align', true, 'verbose', verbose);
    else
        Data = pshr_load('HR', hr_files, 'Affect', aff_files, 'align', true, 'verbose', false, 'meta', 15);
    end
    
    if ecg_features
        Data.ECG.Aligned_Metrics = {};
    end
      
    %% RR-interval preprocessing
    for i = 1:length(hr_files)
        % Create new copy of Data.HR.Raw for PreProcessing
        Data.HR.PP{i} = Data.HR.Raw{i};
        
        % Apply basic bandpass to HR data
        Data.HR.PP{i}(:,3) = bandpass(Data.HR.PP{i}(:,3), 300, 1600, false);
        new_tabs{q} = table_combo(Data.HR.Affect{i}, target{:}, 'omit', omit{:});

        
        Data.HR.PP{i} = feature_generation(Data.HR.PP{i}, bin, false);
        
        % If they specified ecg_features, then append them to Data.HR.PP
        % using the 
        if ecg_features
            [ecg_aligned, Data.ECG.Aligned_Metrics{i}] = ecg_rr_align(Data.HR.Raw{i}(:,[1,3]), Data.ECG.Raw{i}(:,[1,3]), 130, 'verbose', verbose);
            disp(ecg_files{i})
            Data.HR.PP{i}(:,end+1:end+3) = ecg_aligned(:,3:5); %Only look at [Q,R,S] complex
        end
        
        % Add a column denoting which timepoints are part of aff_list
        % (currently just a binary 0 or 1, hence NumberCategories = false)
%         Data.HR.PP{i} = affect_mark(Data.HR.PP{i}, Data.HR.Affect{i}, aff_list, 'NumberCategories', false); %mark the affect locations
    end 
    
    
    

    %% Function to make sure that the quantity of problematic and nonproblematic behavior datapoints are approximately the same
    
    % Check that no datapoint with NaNs is being included
    on_mat = [];
    off_mat = [];
    un_mat = [];
    
    % Onset calculation
    for i=1:length(Data.HR.PP)
         [on_, off_, un_] = onset_sample(Data.HR.PP{i}(:,3:end), Data.HR.Affect{i}, aff_list,...
                        'band',[5,0], 'omit_nan',true, 'dilate', 8);
         on_mat(end+1:end+length(on_(:,1)),:) = on_(:,1:end);
         off_mat(end+1:end+length(off_(:,1)),:) = off_(:,1:end);
         un_mat(end+1:end+length(un_(:,1)),:) = un_(:,1:end);
    end
    
    % Duration calculation
    
    

    %% Finally get the actual regression portion of the analysis

%     Data.HR.log_reg={};
%     for i=1:length(Data.HR.PP)
%         disp(Data.HR.files{i});
%         Data.HR.log_reg{i} = gen_regression(Data.HR.PP{i}(:,[3:end-1]), Data.HR.PP{i}(:,end),...
%         'log', 'verbose', true);
%     end
    
    % Run log regression on all data concatenated together
%     big = vertcat(Data.HR.PP{:});
    
    i = 1;
    while i < 200
%         [Data.HR.RegPP(i,:,:), Data.HR.PPhat(i,:,:), Data.HR.Ttest(i,:)] = gen_regression(big(:,[3:end-1]),big(:,end), 'log');
        [Data.HR.RegPP(i,:,:), Data.HR.PPhat(i,:,:), Data.HR.AUC(i,1)] = gen_regression([un_mat;on_mat],[zeros(size(un_mat,1),1); ones(size(on_mat,1),1)],...
            'log', 'verbose', false);
        i = i+1;
    end
    
    % Plot the histograms for the beta values and the confidence intervals
    row = floor(sqrt(size(Data.HR.RegPP,2)));
    col = ceil(size(Data.HR.RegPP,2)/row);
    for i =1:size(Data.HR.RegPP,2)
        subplot(row,col,i)
        histogram(Data.HR.RegPP(:,i,5));
        hold on;
        histogram(Data.HR.RegPP(:,i,6));
        histogram(Data.HR.RegPP(:,i,1));
        legend('upper bound', 'lower bound', 'beta');
        title(strcat('estimated value of coefficient: ', string(i)));
        %boxchart([Data.HR.RegPP(:,i,5),Data.HR.RegPP(:,i,6)]);
    end
    
    % Plot p-values for the coefficents
    figure;
    for i = 1:size(Data.HR.RegPP,2)
        subplot(row,col,i)
        histogram(Data.HR.RegPP(:,i,4));
        title(strcat('p-value for coefficient: ',string(i)));
    end
    
    
    % Plot the accuracy results
    figure;
    titles = {'Non-Prob Accuracy','Prob Accuracy'};
    for i = 1:2
        subplot(1,2,i)
        histogram(Data.HR.PPhat(:,1,i));
        hold on;
        histogram(Data.HR.PPhat(:,2,i));
        histogram(Data.HR.PPhat(:,3,i));
        hold off;
        title(titles{i});
        legend('Lower Bound', 'Pred', 'Upper Bound');
    end
    
    % Plot spread of AUC
    figure;
    histogram(Data.HR.AUC)
    title('AUC Values');
    
end


function [ret,locs] = ecg_preprocess(mat, amp, cut_bin)
% Function to take the ECG data preprocesses it by removing ECG values
% which fall outside of the accepted amplitude
    %inputs:
    %   mat: [n-by-1] vector containing ecg data
    %   amp: [int], minimum amplitude of accepted ECG values. If a value
    %   falls outside of this bound, all values [cut_bin] before and after
    %   it are replaced with NaNs.
    %   cut_bin: [int], the amount of indexes before and after entries
    %   which fail [amp] that are replaced with NaNs.
    
    %Returns:
    %   ret: [n-by-1] matrix containing the ecg data with all removed
    %   values replaced by NaNs
    %   locs: [m-by-1] index of all values which fall outside of the bounds
    %   described by [amp]

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

function [mat] = feature_generation(mat, bin, band)
% Function for generating the different features for multiple recording
% sessions

% Inputs:
%   mat: [n-by-m matrix] 
%   bin: [1-by-2 cell array] The bin type you want to use
%   band: [1-by-2 matrix]

    mat(:,end+1) = rmssd_calc(mat(:,3), bin, band);
    mat(:,end+1) = pnnx_calc(mat(:,3),50, bin, band);
%     mat(:,7) = sdnn_calc(mat(:,3),bin,band);
%     mat(:,8) = sdsd_calc(mat(:,3),bin,band);
end
