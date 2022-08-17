function [] = regression_pipeline(hr_files, ecg_files, aff_files, verbose)
% This function is the default pipeline for regression analysis of
% collected affect data.

% Inputs:
%   hr_files: [1-by-n cell array]
%   ecg_files: [1-by-n cell array]
%   aff_files: [1-by-n cell array]

    aff_list = {'SIB','innappropriate face related behavior','polar strap adjustment/removal'...
        'repetitive behaviors','inappropriate movement','crying'};

    Data = pshr_load('HR', hr_files, 'ECG', ecg_files, 'Affect', aff_files, 'align', true, 'verbose', false);

    %% RR-interval preprocessing
    for i = 1:length(hr_files)
        Data.HR.PP{i} = Data.HR.Raw{i};
        Data.HR.PP{i} = affect_mark(Data.HR.PP{i}, Data.HR.Affect{i},aff_list); %mark the affect locations
        % We'll just work with bandpassing for now...
        Data.HR.PP{i}(:,3) = bandpass(Data.HR.PP{i}(:,3), 300, 1200, false);
    end
    
    
    %% ECG preprocessing
    amp = 5000;
    cut_bins = 100;
    for i = 1:length(ecg_files)
        Data.ECG.PP{i} = Data.ECG.Raw{i}; %PP stands for PreProcessed
        Data.ECG.PP{i} = affect_mark(Data.ECG.PP{i}, Data.ECG.Affect{i}, aff_list);
        [ret, locs] = ecg_preprocess(Data.ECG.PP{i}(:,3),amp,cut_bins);
        Data.ECG.PP{i}(:,3) = ret;
    end
    
    %TODO: Test how PQRST handles NaN values
    
    %% Go through combinations of parameters for best RR/ECG alignment

    for i = 1:length(aff_files)
        %[a,b,c] = ecg_rr_alignment(Data.HR.Raw{1}(:,[1,3]), Data.ECG.Raw{1}(:,[1,3]),700,50,130,10,true);
        [Data.ECG.Aligned{i},Data.ECG.Aligned_metrics{i}] = ecg_rr_align(Data.HR.PP{i}(:,[1,3]), Data.ECG.PP{i}(:,[1,3]),...
            130,'subcost',10,'verbose',false);
    end

%     val_iter_results = [NaN, NaN, NaN, NaN, NaN, NaN];
%     time_iter_results = [NaN, NaN, NaN, NaN, NaN, NaN];
    %peak
%     for i=600:50:900
%         % dist
%         for j=40:10:70
%         % subcost
%             for k= 5:5:20
%                 [a,c] = ecg_rr_alignment(Data.HR.Raw{1}(:,[1,3]), Data.ECG.Raw{1}(:,[1,3]),i,j,130,k,false);
%                 val_iter_results(end+1,:) = [sum(~isnan(c.val.diff))/length(c.val.diff), c.val.mean, c.val.std, i,j,k];
%                 time_iter_results(end+1,:) = [sum(~isnan(c.time.diff))/length(c.time.diff), c.time.mean, c.time.std, i,j,k];
%             end
%         end
%     end


    %% Cobble together what features you want to use for analysis


    %% Function to make sure that the quantity of problematic and nonproblematic behavior
    %  datapoints are approximately the same
    
    % Check that no datapoint with NaNs is being included
    for i=1:length(Data.HR.PP)
        Data.HR.PP{i} = feature_generation(Data.HR.PP{i}, {'5', 'second'}, false);
    end
    

    %% Finally get the actual regression portion of the analysis

    if verbose
        for i=1:length(Data.ECG.Aligned_metrics)
            disp(strcat('Aligned val for :',Data.ECG.files{i},' = ', string(Data.ECG.Aligned_metrics{i}.val.mean),' + ', string(Data.ECG.Aligned_metrics{i}.val.std)));
            disp(strcat('Aligned time for :',Data.ECG.files{i},' = ', string(Data.ECG.Aligned_metrics{i}.time.mean),' + ', string(Data.ECG.Aligned_metrics{i}.time.std)));
        end
    end
    
    for i=1:length(Data.HR.PP)
        disp(Data.HR.files{i});
        results = gen_regression(Data.HR.PP{i}(:,[3:end-1]), Data.HR.PP{i}(:,end),...
        'log');
    end
    
    
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
    mat(:,7) = sdnn_calc(mat(:,3),bin,band);
    mat(:,8) = sdsd_calc(mat(:,3),bin,band);
    
    %move coding into last column
    mat(:,end+1) = mat(:,4);
    mat(:,4) = [];

end
