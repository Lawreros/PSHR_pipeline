function [ecg_aligned, aligned_metrics]=ecg_rr_align(rr, ecg, sample_rate, varargin)
% Function which attempts to align the ecg data with the RR-interval data
% and return metrics to judge alignment
% Inputs:
%   rr: [n-by-2 matrix] matrix which contains the timestamps and RR-interval
%       values in the format [timestamp, RR]. This allows for the comparison of
%       timestamps between RR and ECG
%   ecg: [m-by-2 matrix] matrix which contains the timestamps and ECG
%       values in the format [timestamp,ECG].
%       sample_rate: [int or float] the sampling rate used when collecting the
%       ECG data, required for the estimation of RR-intervals from the ECG
%       data.
%   subcost: [int or float] substitution cost coefficient for use in
%       Levenshtein calculation (see levenshtein_align). Default is 10.
%   verbose: [bool] Whether the function should plot the results of the
%       alignment.
%   
% Additional Inputs:
%   This function utilizes both functions ecg_PQRST and levenshtein_align
%   in its calculation. Thus, any parameters used by either function also
%   accepted by this function and passed onward
%
% Returns:
%   ecg_aligned: [n-by-6 matrix] results from feeding the output
%       `time_matrix` through the function levenshtein_align
%   aligned_metrics: [struct] structure containing different comparison
%       results from the alignment process

    p = inputParser;
    p.KeepUnmatched=true;
    validMatrix = @(x) isnumeric(x) && ismatrix(x);
    addRequired(p,'rr', validMatrix);
    addRequired(p,'ecg',validMatrix);
    addRequired(p, 'sample_rate',@isscalar);
    addParameter(p,'subcost', 10, @isscalar);
    addParameter(p, 'verbose', false, @islogical)
    parse(p,rr,ecg,sample_rate,varargin{:});
    
    
    % Fill in any NaN values for the timestamps of RR-itnervals and ECG
    for i=1:size(rr(:,1),1) % Go through each R-wave index and find the closest corresponding timestamp in mat
        j = 0;
        
        while isnan(rr(i-j,1))
            j = j+1;
        end
        rr(i,1) = rr(i-j,1);
    end
    
    for i=1:size(ecg(:,1),1) % Go through each R-wave index and find the closest corresponding timestamp in mat
        j = 0;
        
        while isnan(ecg(i-j,1))
            j = j+1;
        end
        ecg(i,1) = ecg(i-j,1);
    end
    
    
    %Preprocessing of ECG data to remove discrete noise:
    ecg(:,2) = bandpass(ecg(:,2), -3000, 3000, false);
    rr(:,2) = bandpass(rr(:,2), 300, 1600, false);
    
    % Calculate the RR-interval locations/durations from the ECG data
    [ecg_locs, ecg_times] = ecg_PQRST(ecg,'sample_rate',p.Results.sample_rate,varargin{:});
    
    % With the calculated RR data from ECG, align the two matrices
    [ecg_aligned, lev_transform] = levenshtein_align(rr, 2, ecg_times, 4, p.Results.subcost, 10000);
    
    
    aligned_metrics.lev = lev_transform;
    % Calculate difference in time metrics
    aligned_metrics.time={};
    aligned_metrics.time.diff = rr(:,1)-ecg_aligned(:,1);
    aligned_metrics.time.std = nanstd(aligned_metrics.time.diff);
    aligned_metrics.time.mean = nanmean(abs(aligned_metrics.time.diff));
    
    % Calculate difference in RR-interval metrics
    aligned_metrics.val = {};
    aligned_metrics.val.diff = rr(:,2) - ecg_aligned(:,4);
    aligned_metrics.val.std = nanstd(aligned_metrics.val.diff);
    aligned_metrics.val.mean = nanmean(aligned_metrics.val.diff);
    
    if p.Results.verbose
        figure;
        plot(aligned_metrics.time.diff/1000);
        title("(RR-interval timestamp) - (ECG R-peak timestamp)");
        xlabel("RR-interval index");
        ylabel("Difference (sec)");

        figure;
        plot(aligned_metrics.val.diff);
        title("(RR-interval) - (ECG RR estimate)");
        xlabel("RR-interval index");
        ylabel("Difference (millisecond)");
        
        %histogram for time and RR-interval differnces
        figure;
        subplot(1,2,1);
        histogram(aligned_metrics.time.diff/1000,[-1,0,1,2,3,4],'Normalization','probability');
        title("(RR-interval timestamp) - (ECG R-peak timestamp)");
        
        subplot(1,2,2);
        histogram(aligned_metrics.val.diff,[-70:10:70],'Normalization','probability');
        title("(RR-interval) - (ECG RR estimate)");
        
        figure;
        subplot(2,1,1)
        plot(rr(:,1), rr(:,2));
        hold on;
        plot(ecg_aligned(:,1)+aligned_metrics.time.mean, ecg_aligned(:,4));
        xlim([rr(1,1), rr(end,1)]);
        title("RR-interval and RR-interval from ECG");
        xlabel("timestamp (ms)");
        ylabel("RR-interval (ms)");
        legend('RR-interval','RR-interval from ECG');
        
        subplot(2,1,2)
        plot(ecg(1:3:end,1),ecg(1:3:end,2));
        xlim([rr(1,1), rr(end,1)]);
        title("ECG reference data");
        xlabel("timestamp (ms)");
        ylabel("Voltage (uV)");
    end

end