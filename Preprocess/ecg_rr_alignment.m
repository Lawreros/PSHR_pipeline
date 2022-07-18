function [aligned_values]=ecg_rr_alignment(rr, ecg)
% Function which attempts to align the ecg data with the RR-interval data
% Inputs:
%   rr: [n-by-2 matrix] matrix which contains the timestamps and
%   RR-interval values in the format [timestamp, RR]. This allows for the
%   comparison of timestamps between RR and ECG
%   ecg: [n-by-2 matrix] matrix which contains the timestamps and
%   ECG values in the format [timestamp, RR].


    % Get the rr-interval estimates from the ECG data
    peak = 800;
    dist = 40;
    freq = 130;
    ecg_rr = ecg_rr_conversion(ecg, peak, dist, freq);
    ecg_times = ecg_rr(:,1);
    ecg_rr = ecg_rr(:,3);
    
    % Create matrix for timestamp comparison
    rr_times = rr(:,1);
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
                substitutionCost = abs(rr(i) - ecg_rr(j))/10;
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
            aligned_values(end+1) = ecg_rr(j);
            aligned_times(end+1) = ecg_times(j);
            j = j+1;
        elseif move(1,i) == 1
            j = j+1;
        elseif move(1,i) == 2
            aligned_values(end+1) = NaN;
            aligned_times(end+1) = ecg_times(j);
        end

    end

    
    %Append RR-values to bottom
    aligned_values(2,:) = flip(rr_times,2);
    
    % TODO: Provide some analysis metrics for the alignment steps, i.e.
    % metrics for variable `move`
    
    
    
%editDistance calculates the maximum number of insertions/deletions in
%order to convert the first string into the second string
% TODO: Look into the cost function for editDistance, want to make
% deletions more common than insertions

%d = editDistance(ecg,rr);

end
