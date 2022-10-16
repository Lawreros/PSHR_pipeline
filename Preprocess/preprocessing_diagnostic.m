function preprocessing_diagnostic(mat, label)
% Function for analyzing the extent that different preprocessing
% pipelines have on the resulting data [WARNING: ONLY WORKS ON BINARY LABLING
% AT THE MOMENT]
    
% Inputs:
%   mat: [n-by-1 matrix] Vector containing the preprocessed data
%   label: [n-by-1 matrix] Vector of labels for each data point


    % Percentage of problematic behavior that is NaNs
    nanmat = isnan(mat)+label;
    
    % Number of removed data points:
    rem = sum(nanmat > 1);
    disp(strcat('Removed number of non-PB datapoints: ',string(sum(isnan(mat))-rem),' of ', string(sum(label==0))));
    disp(strcat('Removed number of PB datapoints: ', string(rem),' of ', string(sum(label==1))));
    disp(string(rem/sum(label==1)));
    
    % Percentage PB's with NaNs at onset
    lab_idx = grp2idx(categorical(label)); % Converts label into categories
    start_idx = find([0;lab_idx]-[lab_idx;0]); % Finds points where the label changes
    % Every odd entry below is the start of an onset
    start_idx(end) = length(lab_idx);
    
    count = 0;
    for i =1:2:length(start_idx)
        if isnan(mat(start_idx(i),1))
            count = count+1;
        end
    end
    disp('PB with NaNs at onset');
    disp(strcat(string(count), ' of ', string(length(start_idx)/2), ' removed'));
    disp(string(count/(length(start_idx)/2)));
    
    % Percentage of PB's with NaNs within X indices of onset
    count = 0;
    skp = 0;
    x = 5;
    for i =1:2:length(start_idx)
        
        if start_idx(i) > x
            if isnan(sum(mat(start_idx(i)-x:start_idx(i),1)))
                count = count+1;
            end
        else
            skp = skp+1;
            %disp('skipped');
        end
    end
    disp(strcat('PB with NaNs within ', string(x), ' indices of onset'));
    disp(strcat(string(count), ' of ', string((length(start_idx)/2)-skp), ' removed'));
    disp(string(count/((length(start_idx)-skp)/2)));

end