function [aligned_values, aligned_times, aligned_metrics]=ecg_rr_alignment(rr, ecg, peak, dist, freq, subcost, verbose)
% DEPRICATED: PLEASE USE ecg_rr_align INSTEAD



% Function which attempts to align the ecg data with the RR-interval data
% Inputs:
%   rr: [n-by-2 matrix] matrix which contains the timestamps and
%   RR-interval values in the format [timestamp, RR]. This allows for the
%   comparison of timestamps between RR and ECG
%   ecg: [n-by-2 matrix] matrix which contains the timestamps and
%   ECG values in the format [timestamp, ECG].
%   peak: [int], minimum peak prominence for peak
%   dist: [int], minimum distance in index number for peaks
%   freq: [int], sampling frequency for ECG data in Hz
%   subcost: [int] substitution cost coefficient for use in Levenshtein
%   calculation. The equation for cost is: abs(X - Y)/subcost. Increasing
%   subcost will weight against substitution/deletion. Recomended value is
%   10.
%   verbose: [bool] Whether the function should plot the results of
%   alignment

% Returns:
%   aligned_values: [m-by-2 matrix] of the aligned ECG and RR-interval
%   data, organized as [RR, ECG]
%   aligned_times: [m-by-2 matrix] timestamps corresponding to the entries
%   in aligned_values for the purposes of comparing lag/alignment
%   aligned_metrics: [struct] structure containing different comparison
%   results from the alignment process


    % Get the rr-interval estimates from the ECG data
%     peak = 800;
%     dist = 40;
%     freq = 130;
    ecg_rr = ecg_rr_conversion(ecg, peak, dist, freq);
    ecg_times = ecg_rr(:,1);
    ecg_rr = ecg_rr(:,4);
    
    % Create matrix for timestamp comparison
    rr_times = rr(:,1);
    
    % Replace NaN timestamp values with closest previous timestamp
    for i=1:length(rr_times)
        j = 0;
        while isnan(rr_times(i-j,1))
            j = j+1;
        end
        rr_times(i,1) = rr_times(i-j,1);
    end
    
    rr = rr(:,2);

    % Create matrix for alignment values for levenshtein distance
    lev = zeros(length(rr),length(ecg_rr));
    lev(1,1) = 0;

    for i = 2:length(rr)
        lev(i,1) = i-1;
    end
    for j = 2:length(ecg_rr)
        lev(1,j) = j-1;
    end

    for j = 2:length(ecg_rr)
        for i = 2:length(rr)
            if rr(i) == ecg_rr(j)
                substitutionCost = 0;
            else
                substitutionCost = abs(rr(i) - ecg_rr(j))/subcost;
            end
            lev(i,j) = min([lev(i-1,j)+1, lev(i,j-1)+1, lev(i-1, j-1)+substitutionCost]);
        end
    end

    [Y,I] = min(lev(end,:));

    [r,c] = size(lev);
    c = I;

    move = [];
    q_c = 0;
    q_r = 0;

    while r > 0

        if c>1 && r>1
            q = min([lev(r,c-1),lev(r-1,c-1), lev(r-1,c)]);
        elseif r == 1 && c~=1
            q = lev(r,c-1);
            q_r = 1;
        elseif c == 1 && r ~=1
            q = lev(r-1,c);
            q_c = 1;
        end

        if r == 1
            move(end+1:end+c-1) = 1;
            r= 0;
        elseif c == 1
            move(end+1:end+r-1) = 2;
            r = 0;
        elseif q == lev(r-1,c-1)
            move=[move,0]; %diag
            r=r-1;
            c=c-1;
        elseif q == lev(r,c-1) && q_c ~=1
            move=[move,1]; %left
            c=c-1;
        elseif q == lev(r-1,c) && q_r ~=1
            move=[move,2]; %up
            r = r-1;
        elseif q_c == 1 && q_r == 1
            r=r-1;
            c=c-1;
        end

    end
    move = [move,0];

    % Flip move
    move = flip(move,2);

    clear i j c r I line pks_a q q_c q_r Y substitutionCost;

    % Realign
    j = 1;
    aligned_values=[];
    aligned_times=[];

    for i = 1:length(move)
        if move(1,i) == 0
            aligned_values(end+1,2) = ecg_rr(j);
            aligned_times(end+1,2) = ecg_times(j);
            j = j+1;
        elseif move(1,i) == 1
            j = j+1;
        elseif move(1,i) == 2
            aligned_values(end+1,2) = NaN;
            aligned_times(end+1,2) = NaN;
        end

    end

    
    %Append RR-values to bottom
    aligned_values(:,1) = rr;
    aligned_times(:,1) = rr_times;
    
    % TODO: Provide some analysis metrics for the alignment steps, i.e.
    % metrics for variable `move`
    aligned_metrics = {};
    
    %Calculate difference in time metrics
    aligned_metrics.time = {};
    aligned_metrics.time.diff = aligned_times(:,1)-aligned_times(:,2);
    aligned_metrics.time.std = nanstd(aligned_metrics.time.diff);
    aligned_metrics.time.mean = nanmean(abs(aligned_metrics.time.diff));
                
    %Calculate difference in RR-interval metrics
    aligned_metrics.val = {};
    aligned_metrics.val.diff = aligned_values(:,1)-aligned_values(:,2);
    aligned_metrics.val.std = nanstd(aligned_metrics.val.diff);
    aligned_metrics.val.mean = nanmean(aligned_metrics.val.diff);
    
        
    if verbose % print/plot alignment metrics
        figure(1);
        plot(aligned_metrics.time.diff/1000);
        title("(RR-interval timestamp) - (ECG R-peak timestamp)");
        xlabel("RR-interval index");
        ylabel("Difference (sec)");

        figure(2);
        plot(aligned_metrics.val.diff);
        title("(RR-interval) - (ECG RR estimate)");
        xlabel("RR-interval index");
        ylabel("Difference (millisecond)");
    end
    
    
%editDistance calculates the maximum number of insertions/deletions in
%order to convert the first string into the second string
% TODO: Look into the cost function for editDistance, want to make
% deletions more common than insertions

%d = editDistance(ecg,rr);

end
