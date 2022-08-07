function [locs] = ecg_rr_conversion(mat, peak, dist, freq)
% Function to take the ECG data and estimate RR-intervals from them
    %inputs:
    %   mat: [n-by-2] vector containing ecg data and timestamps in the
    %   format of [timestamp, ecg]
    %   peak: [int], minimum peak prominence for peak
    %   dist: [int], minimum distance in index number for peaks
    %   freq: [int], sampling frequency for ECG data in Hz
    
    %Returns:
    %   locs: [m-by-4] matrix containing the [timestamp, index, voltage, RR-calc]
    %   for the peaks


    [pks, locs] = findpeaks(mat(:,2), 'MinPeakProminence', peak, 'MinPeakDistance', dist);
    
    locs(:,2) = locs(:,1);
    locs(:,1) = NaN;
    locs(:,3) = pks;
    
    for i = 2:length(locs)
        %convert to milliseconds assuming a <freq>Hz sampling frequency
        locs(i,4) = (locs(i,2)-locs(i-1,2))*(1000/freq);
        
        %record timestamp for each peak
        j = 0;
        while isnan(mat(locs(i,2)-j,1))
            j = j+1;
        end
        locs(i,1) = mat(locs(i,2)-j,1);
    end

end

