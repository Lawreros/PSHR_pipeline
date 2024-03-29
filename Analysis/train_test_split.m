function [train, test, unused] = train_test_split(mat, label, ratio, varargin)
% Function which, when given a vector of different categories, and a ratio
% of different types of each label, splits a matrix into train/test/unused
% groupings.

% Required Inputs:
%   mat: [n-by-k matrix] Matrix consisting of the data you wish to seperate
%       into training and testing parameters
%
%   label: [n-by-1 matrix] Vector of integer labels for all rows of mat,
%       the number will correspond to the field name in the train and test
%       structs.
%
%   ratio: [1-by-m matrix] Matrix containing the ratio of each category
%       that you wish to have in both the training and testing subsets. The
%       number of ratio entries should be equal to the number of unique
%       values in 'label'. For example, inputing [0.3, 0.7] means that the
%       struct 'train' should have 30% category 0 and 70% category 1 (same
%       with struct 'test'). All datapoints omitted in order to make this
%       ratio will be put into struct 'unused', which does not obey this
%       ratio.


% Optional Parameters:
%   split: [float] The train/test split percentage as a decimal. 0.8 will
%       mean that 80% of the data points are put into the training set
%
%   ignore_nan: [bool] Whether to ignore any rows of mat or label with NaNs
%       in them. Default is true.

% Returns:
%   train: [struct] training data, split up into multiple fields by 
%       numerical affect label
%
%   test: [struct] testing data, split up into multiple fields by 
%       numerical affect label
%
%   unused: [struct] extra data not put into either the training or testing
%       sets due to the splitting and ratio process


    % Get the different catagories
        
    
    p = inputParser;
    addParameter(p, 'split', 0.8, @isscalar);
    addParameter(p, 'ignore_nan', true, @islogical);
    addParameter(p, 'verbose', false, @islogical);
    parse(p,varargin{:});
    
    
    if p.Results.ignore_nan % Remove any entries where a NaN is in mat or label
        label(any(isnan(label),2),:)=[];
        mat(any(isnan(label),2),:)=[];
        
                
        label(any(isnan(mat),2),:) = [];
        mat(any(isnan(mat),2),:) = [];
    end
    
    
    [ii, jj, kk] = unique(label); % [categories, locations, sorted output]
    f = histc(kk,1:numel(jj)).'; % number of each category
    if p.Results.verbose
        disp('Category quantity: [Non-problematic   Problematic]');
        disp(f);
    end
    
    % Check that provided ratio matches number of categories
    if length(ratio) ~= length(f)
        disp('ERROR: Provided ratios for categories: ', string(length(ratio)),...
            ' does not match number of categories present: ',string(length(f)));
        return
    end

    
    % Take the smallest category and use that as a reference for the others
    
    % other_categories / smallest = new count
    % Check that works, if not then iterate through subtracting 1 from the
    % smaller category until things work? No
    
    % Calculate max difference between desired ratio and given ratio, fit
    % that then run again
    
    given_ratio = f/sum(f);
    req_change = 2;
    new_f = NaN;
    
    while req_change >= 1
        % find smallest category
        [val, idx] = min(f);

        % find what all other values should be from ratio
        for i =1:length(ratio)
            new_f(i) = (ratio(i)/ratio(idx)) * val;
        end

        [req_change, idz] = max(new_f - f);
        f(idx) = f(idx) - (req_change * ratio(idx)/ratio(idz));
    end
    
    % floor new_f so you don't deal with float issues, there is no way
    % rounding down will require you to access more data
    new_f = floor(new_f);
    
    % Break the collected data into submatrices for the output structure:
    
    train = {};
    test = {};
    unused = {};
    
    for i = 1:length(new_f)
        dump = mat(label==ii(i),:);
        [r,c] = size(dump);
        
        subset = randsample(r,new_f(i));
        train_subset = subset(1:floor(p.Results.split*new_f(i)));
        test_subset = subset(ceil(p.Results.split*new_f(i)):end);
        
        % Get indices of all the datapoints that weren't used for the
        % unused structure
        not_samp = [1:r];
        not_samp(subset) = [];
        
        train.('cat_'+string(ii(i))) = dump(train_subset,:);
        test.('cat_'+string(ii(i))) = dump(test_subset,:);
        unused.('cat_'+string(ii(i))) = dump(not_samp, :);
        
    end
    
        
    % Return struct where it is train.(category) = sub_matrix
    % Also return test.(category) = sub_matrix or NaN if there is nothing
    % left
end