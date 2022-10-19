function [] = random_forest_pipeline(hr_files,aff_files,varargin)
% Pipeline for the creation of a random forest model
% Required Inputs:
%   mat: [1-by-n cell array]
%
% Optional Parameters:
%   disp_plot: [bool] Whether to plot the data in mat with the
%       different waves indicated in the figure. Default is false.

% Returns:
%   a

    p = inputParser;
    addParameter(p,'ordinal',false, @islogical);
    addParameter(p,'verbose',true, @islogical);
    
    parse(p,varargin{:});
    
    
    % Load in data
     aff_list = {'SIB','ISB','innappropriate face related behavior','polar strap adjustment/removal'...
        'repetitive behaviors','inappropriate movement','crying', 'pulling at pants'};
    
    Data = pshr_load('HR', hr_files, 'Affect', aff_files, 'align', true, 'verbose', false);

    
    %% RR-interval preprocessing
    for i = 1:length(hr_files)
        Data.HR.PP{i} = Data.HR.Raw{i};
        Data.HR.PP{i} = affect_mark(Data.HR.PP{i}, Data.HR.Affect{i},aff_list); %mark the affect locations
        % We'll just work with bandpassing for now...
        Data.HR.PP{i}(:,3) = bandpass(Data.HR.PP{i}(:,3), 300, 1600, false);
        %Data.HR.PP{i}(:,3) = acar(Data.HR.PP{i}(:,3), 5, false);
        %Data.HR.PP{i}(:,3) = kamath(Data.HR.PP{i}(:,3),false);
        %Data.HR.PP{i}(:,3) = karlsson(Data.HR.PP{i}(:,3),false);
        %Data.HR.PP{i}(:,3) = malik(Data.HR.PP{i}(:,3),false);
    end
    

    %% Function to make sure that the quantity of problematic and nonproblematic behavior datapoints are approximately the same
    
    if onset %If you want to do the onset analysis data for training
    
        on_mat = [];
        off_mat = [];
        un_mat = [];
        for i=1:length(Data.HR.PP)
             Data.HR.PP{i} = feature_generation(Data.HR.PP{i}, {[5,0], 'second'}, false);

             Data.HR.PP{i}(:,end) = []; %get rid of affect labeling
             [on_, off_, un_] = onset_sample(Data.HR.PP{i}(:,[1,3:end]), Data.HR.Affect{i}, aff_list,...
                            'band',[4,0], 'omit_nan',true, 'dilate', 30);
             on_mat(end+1:end+length(on_(:,1)),:) = on_(:,1:end);
             off_mat(end+1:end+length(off_(:,1)),:) = off_(:,1:end);
             un_mat(end+1:end+length(un_(:,1)),:) = un_(:,1:end);
        end


        i = 1;
        conf_mat = [];
        feat_imp = [];
        while i < 400
            [train, test, unused] = train_test_split([on_mat;un_mat], [zeros(size(un_mat,1),1); ones(size(on_mat,1),1)], [0.5 0.5], 'split', 0.8);

            [a,b,c,d,f] = random_forest([train.cat_0;test.cat_0;train.cat_1;test.cat_1],...
                                        [zeros(size([train.cat_0;test.cat_0],1),1);ones(size([train.cat_1;test.cat_1],1),1)],...
                                        100, 0.2,{});%, {'RR-interval','RMSSD','pNN50'});


    %         [a,b,c,d] = random_forest([on_mat;un_mat],[zeros(size(un_mat,1),1); ones(size(on_mat,1),1)], 100, 0.2, {'RR-interval','RMSSD','pNN50'});
            conf_mat(end+1,:) = d;
            feat_imp(end+1,:) = f;
            i = i+1;
        end
    
    else %Just fit the data for the entire dataset, making sure to do a train-test split
        
        for i=1:length(Data.HR.PP)
            Data.HR.PP{i} = feature_generation(Data.HR.PP{i}, {[5,0], 'second'}, false);
        end
        
        big = vertcat(Data.HR.PP{:});
        
        i = 1;
        conf_mat = [];
        feat_imp = [];
        while i < 400
            [train, test, unused] = train_test_split(big(:,[1:end-1]),big(:,end), [0.5 0.5], 'split', 0.8);

            [a,b,c,d,f] = random_forest([train.cat_0;test.cat_0;train.cat_1;test.cat_1],...
                                    [zeros(size([train.cat_0;test.cat_0],1),1);ones(size([train.cat_1;test.cat_1],1),1)],...
                                    100, 0.2,{});
            conf_mat(end+1,:) = d;
            feat_imp(end+1,:) = f;
            i = i+1;
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
