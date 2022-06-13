function [locs] = ecg_rr_conversion(mat, peak, dist, freq)
% Function to take the ECG data and estimate RR-intervals from them
    %inputs:
    %   mat: [n-by-1] vector containing ecg data
    %   peak: [int], minimum peak prominence for peak
    %   dist: [int], minimum distance in index number for peaks
    %   freq: [int], sampling frequency for ECG data in Hz


    [pks, locs] = findpeaks(mat, 'MinPeakProminence', peak, 'MinPeakDistance', dist);

    for i = 2:length(locs)
        %convert to milliseconds assuming a <freq>Hz sampling frequency
        locs(i,2) = (locs(i,1)-locs(i-1,1))*(1000/freq);
    end

end

