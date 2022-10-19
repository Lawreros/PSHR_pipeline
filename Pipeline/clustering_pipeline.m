function [idx] = clustering_pipeline(hr_files,aff_files,varargin)
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
    
     for i=1:length(Data.HR.PP)
        Data.HR.PP{i} = feature_generation(Data.HR.PP{i}, {[5,0], 'second'}, false);
     end
        
     big = vertcat(Data.HR.PP{:});
        
%      i = 1;
%      while i < 400
%         [train, test, unused] = train_test_split(big(:,[3:end-1]),big(:,end), [0.5 0.5], 'split', 0.8);

     [idx] = newFdbscan(big(:,3:end-1), {'RR-interval','RMSSD','pNN50','SDNN','SDSD'}, big(:,end), 50, 10, true);
     
     dats = unique(idx);
     for i = 2:length(dats)
         disp(strcat('Percentage of points belonging to cluster ',string(dats(i)),': ', string(sum(idx==dats(i))*100/length(idx)),'%'))
     end
     
     disp(strcat('Percentage of unassigned datapoints: ', string(sum(idx==-1)*100/length(idx)),'%'));
%         i = i+1;
%      end
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
