function [Data] = regression_pipeline(hr_files, ecg_files, aff_files, verbose)
% This function is the default pipeline for regression analysis of
% collected affect data.

% Inputs:
%   hr_files: [1-by-n cell array]
%   ecg_files: [1-by-n cell array]
%   aff_files: [1-by-n cell array]

    aff_list = {'SIB','ISB','inappropriate face related behavior','polar strap adjustment/removal'...
        'repetitive behaviors','inappropriate movement','crying', 'pulling at pants'};
    
    Data = pshr_load('HR', hr_files, 'Affect', aff_files, 'align', true, 'verbose', false);
%   Data = pshr_load('HR', hr_files, 'ECG', ecg_files, 'Affect', aff_files, 'align', false, 'verbose', false);

%     for i = 1:length(ecg_files)
%         [ecg_aligned, aligned_metrics] = ecg_rr_align(Data.HR.Raw{i}(:,[1,3]), Data.ECG.Raw{i}(:,[1,3]), 130, 'verbose', true);
%         disp(ecg_files{i})
%         aligned_metrics.val
%         aligned_metrics.time
%     end
%     
    %% RR-interval preprocessing
    for i = 1:length(hr_files)
        Data.HR.PP{i} = Data.HR.Raw{i};
        Data.HR.PP{i} = affect_mark(Data.HR.PP{i}, Data.HR.Affect{i},aff_list,'NumberCategories',false); %mark the affect locations
        % We'll just work with bandpassing for now...
        %Data.HR.PP{i}(:,3) = bandpass(Data.HR.PP{i}(:,3), 300, 1600, false);
        %Data.HR.PP{i}(:,3) = acar(Data.HR.PP{i}(:,3), 5, false);
        %Data.HR.PP{i}(:,3) = kamath(Data.HR.PP{i}(:,3),false);
        %Data.HR.PP{i}(:,3) = karlsson(Data.HR.PP{i}(:,3),false);
        %Data.HR.PP{i}(:,3) = malik(Data.HR.PP{i}(:,3),false);
    end 
    
    %colored_lineplot(Data.HR.PP{1}(:,3),Data.HR.PP{1}(:,4))


    %% Function to make sure that the quantity of problematic and nonproblematic behavior datapoints are approximately the same
    
    % Check that no datapoint with NaNs is being included
    on_mat = [];
    off_mat = [];
    un_mat = [];
    
    Data.ECG.Aligned_Metrics = {};
    
    for i=1:length(Data.HR.PP)
         Data.HR.PP{i} = feature_generation(Data.HR.PP{i}, {[5,0], 'second'}, false);
         
%          Data.HR.PP{i}(:,end) = []; %get rid of affect labeling
%          [ecg_aligned, Data.ECG.Aligned_Metrics{i}] = ecg_rr_align(Data.HR.PP{i}(:,[1,3]),Data.ECG.Raw{i}(:,[1,3]), 130, 'verbose', false, 'disp_plot', true);
%          Data.HR.PP{i}(:,end+1:end+3) = ecg_aligned(:,3:5); %Only look at [Q,R,S] complex
         
         [on_, off_, un_] = onset_sample(Data.HR.PP{i}(:,3:end), Data.HR.Affect{i}, aff_list,...
                        'band',[4,0], 'omit_nan',true, 'dilate', 10);
         on_mat(end+1:end+length(on_(:,1)),:) = on_(:,1:end);
         off_mat(end+1:end+length(off_(:,1)),:) = off_(:,1:end);
         un_mat(end+1:end+length(un_(:,1)),:) = un_(:,1:end);
    end
    

    %% Finally get the actual regression portion of the analysis

%     Data.HR.log_reg={};
%     for i=1:length(Data.HR.PP)
%         disp(Data.HR.files{i});
%         Data.HR.log_reg{i} = gen_regression(Data.HR.PP{i}(:,[3:end-1]), Data.HR.PP{i}(:,end),...
%         'log', 'verbose', true);
%     end
    
    % Run log regression on all data concatenated together
    big = vertcat(Data.HR.PP{:});
    
%     preprocessing_diagnostic(big(:,3), big(:,end));
    
    
    i = 1;
    while i < 500
%         [Data.HR.RegPP(i,:,:), Data.HR.PPhat(i,:,:), Data.HR.Ttest(i,:)] = gen_regression(big(:,[3:end-1]),big(:,end), 'log');
        [Data.HR.RegPP(i,:,:), Data.HR.PPhat(i,:,:), Data.HR.Ttest(i,:), Data.HR.AUC(i,1)] = gen_regression([on_mat;un_mat],[zeros(size(un_mat,1),1); ones(size(on_mat,1),1)], 'log');
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
        %boxchart([Data.HR.RegPP(:,i,5),Data.HR.RegPP(:,i,6)]);
    end
    
    % Plot p-values for the coefficents
    figure;
    for i = 1:size(Data.HR.RegPP,2)
        subplot(row,col,i)
        histogram(Data.HR.RegPP(:,i,4));
    end
    
    
    % Plot the accuracy results
    figure;
    for i = 1:2
        subplot(1,2,i)
        histogram(Data.HR.PPhat(:,1,i));
        hold on;
        histogram(Data.HR.PPhat(:,2,i));
        histogram(Data.HR.PPhat(:,3,i));
        hold off;
        legend('Lower Bound', 'Pred', 'Upper Bound');
    end
    
    % Plot spread of AUC
    figure;
    histogram(Data.HR.AUC)
    
    disp('done regression');

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

    mat(:,5) = rmssd_calc(mat(:,3), bin, band);
    mat(:,6) = pnnx_calc(mat(:,3),50, bin, band);
%     mat(:,7) = sdnn_calc(mat(:,3),bin,band);
%     mat(:,8) = sdsd_calc(mat(:,3),bin,band);
    
    %move coding into last column
    mat(:,end+1) = mat(:,4);
    mat(:,4) = [];

end
