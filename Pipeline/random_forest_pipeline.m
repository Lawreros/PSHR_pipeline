function [conf_mat, feat_imp] = random_forest_pipeline(hr_files,aff_files,varargin)
% Pipeline for the creation of a random forest model
% Required Inputs:
%   hr_files: [1-by-n cell array]
%
%   aff_files: []
%
% Optional Parameters:
%   bin: [1-by-2 cell array] Used for creating a vector of the
%       feature results from a sliding bin of Y seconds or entries. This takes the
%       format of {[before, after], 'units'}, so if you want to have a bin 
%       of 5 seconds before (including current RR-interval) and 3 seconds 
%       after: {[5,3], 'second'} or if you want the 5 entries before and 
%       3 entries after the index: {[5,3], 'measure'}
%       If you don't want this, set bin to false.
%   
%   tree_num: [int] an integer that represents the number of trees the
%       model needs to contain when training. If false, then it is set
%       to default (100 trees)
%
%   hold_out: [float] a decimal between (0,1) that represent the percentage
%       of data that's held out for testing
%
%   duration: [bool] Whether to run the random forest analysis using data
%       collected from the timepoints during the target affects against
%       comparison to contol data from timepoints outside of the target
%       affects. If false, then this analysis is not run.
%
%   onset: [bool]
%
%   band: [1-by-2 vector] The number of entries before and after the onset
%       to include in the returned matrix. Default is [3,0] for the three entries
%       before the onset (this will results in 4 values being used, with 
%       the last value being the first timepoint of the affect).
%
%   iterations: [int]
%

% Returns:
%   conf_mat: [1-by-4 matrix] confusion matrix in the form:
%       [True Positives, False Negatives, False Positives, True Negatives]
%
%   feat_imp: [] 
%

    p = inputParser;
    addParameter(p,'bin',{[5,0],'second'}, @iscell);
    addParameter(p,'tree_num', 100, @isscalar);
    addParameter(p,'hold_out', 0.2, @isscalar);
    addParameter(p,'duration', true, @islogical);
    addParameter(p,'onset', true, @islogical);
    addParameter(p, 'band',[3,0], @ismatrix);
    addParameter(p,'iterations',100,@isscalar);
    addParameter(p,'verbose',true, @islogical);
    
    parse(p,varargin{:});
    
    
    % Load in data
     aff_list = {'SIB','ISB','innappropriate face related behavior','polar strap adjustment/removal'...
        'repetitive behaviors','inappropriate movement','crying', 'pulling at pants'};
    
    Data = pshr_load('HR', hr_files, 'Affect', aff_files, 'align', true, 'verbose', false);

    
    %% RR-interval preprocessing
    for i = 1:length(hr_files)
        Data.HR.PP{i} = Data.HR.Raw{i};
        % We'll just work with bandpassing for now...
        Data.HR.PP{i}(:,3) = bandpass(Data.HR.PP{i}(:,3), 300, 1600, false);
        [Data.HR.PP{i}, labs] = feature_generation(Data.HR.PP{i}, p.Results.bin, false, [1,1,0,0]);
        Data.HR.PP{i} = affect_mark(Data.HR.PP{i}, Data.HR.Affect{i},aff_list); %mark the affect locations
    end
    labs = ['Time', 'BPM', 'RR-interval', labs];

    %% Function to make sure that the quantity of problematic and nonproblematic behavior datapoints are approximately the same
    
    if p.Results.duration %Just fit the data for the entire dataset, making sure to do a train-test split
        
%         for i=1:length(Data.HR.PP)
%             Data.HR.PP{i} = feature_generation(Data.HR.PP{i}, {[5,0], 'second'}, false, [1,1,0,0]);
%             Data.HR.PP{i} = affect_mark(Data.HR.PP{i}, Data.HR.Affect{i},aff_list); %mark the affect locations
%         end
        
        big = vertcat(Data.HR.PP{:});
        
        i = 1;
        acc_mat = [];
        conf_mat = [];
        feat_imp = [];
        while i <= p.Results.iterations
            [train, test, unused] = train_test_split(big(:,[1:end-1]),big(:,end), [0.5 0.5], 'split', 0.8);

            [acc,b,c,d,f] = random_forest([train.cat_0;test.cat_0;train.cat_1;test.cat_1],...
                                    [zeros(size([train.cat_0;test.cat_0],1),1);ones(size([train.cat_1;test.cat_1],1),1)],...
                                    p.Results.tree_num, p.Results.hold_out, true);
            
            acc_mat = [acc_mat,acc];                    
            conf_mat(end+1,:) = d;
            if f
                feat_imp(end+1,:) = f;
            end
            
            if ~rem(i,5)
                disp(strcat('Random forest duration iteration: ', string(i), ' complete'));
            end
            
            i = i+1;
        end
        
        if feat_imp
            si = size(test.cat_0,2)-1;  %since hr_data has is_affect, need to exclude
%             labs =  {'time','BPM','RR-interval','RMSSD','pNN50'}; %feature_names;

            %feature importance
            figure
            for j = 1:si
                subplot(1, si, j);
                histogram(feat_imp(:,j) ./ sum(abs(feat_imp),2));
                title(strcat('Normalized feature importance for: ',labs{j+1},' (Duration)'));
            end
        end
        
        figure
        histogram(acc_mat);
        title("Accuracy of Random Forest model (Duration)");
    end
    
    
    
    if p.Results.onset %If you want to do the onset analysis data for training
        on_mat = [];
        off_mat = [];
        un_mat = [];
        for i=1:length(Data.HR.PP)
%              Data.HR.PP{i} = feature_generation(Data.HR.PP{i}, {[5,0], 'second'}, false, [1,1,0,0]);
%              Data.HR.PP{i} = affect_mark(Data.HR.PP{i}, Data.HR.Affect{i},aff_list); %mark the affect locations

%              Data.HR.PP{i}(:,end) = []; %get rid of affect labeling
             [on_, off_, un_] = onset_sample(Data.HR.PP{i}(:,[1:end-1]), Data.HR.Affect{i}, aff_list,...
                            'band', p.Results.band, 'omit_nan',true, 'dilate', 20);
             on_mat(end+1:end+length(on_(:,1)),:) = on_(:,1:end);
             off_mat(end+1:end+length(off_(:,1)),:) = off_(:,1:end);
             un_mat(end+1:end+length(un_(:,1)),:) = un_(:,1:end);
        end


        i = 1;
        acc_mat = [];
        conf_mat = [];
        feat_imp = [];
        while i <= p.Results.iterations
            [train, test, unused] = train_test_split([on_mat;un_mat], [zeros(size(un_mat,1),1); ones(size(on_mat,1),1)], [0.5 0.5], 'split', 0.8);

            [acc,b,c,d,f] = random_forest([train.cat_0;test.cat_0;train.cat_1;test.cat_1],...
                                        [zeros(size([train.cat_0;test.cat_0],1),1);ones(size([train.cat_1;test.cat_1],1),1)],...
                                        p.Results.tree_num, p.Results.hold_out, true);

            acc_mat = [acc_mat,acc]; 
            conf_mat(end+1,:) = d;
            if f
                feat_imp(end+1,:) = f;
            end
            
            if ~rem(i,5)
                disp(strcat('Random forest onset iteration: ', string(i), ' complete'));
            end
            i = i+1;
        end
        
        if feat_imp
            si = size(test.cat_0,2)-1;  %since hr_data has is_affect, need to exclude
%             labs =  {'time','BPM','RR-interval','RMSSD','pNN50'}; %feature_names;

            %feature importance
            figure
            for j = 1:si
                subplot(1, si, j);
                histogram(feat_imp(:,j) ./ sum(abs(feat_imp),2));
                title(strcat('Normalized feature importance for: ',labs{j+1},' (Onset)'));
            end
        end
        
        figure
        histogram(acc_mat);
        title("Accuracy of Random Forest model (Onset)");
        
        
    end
    
end