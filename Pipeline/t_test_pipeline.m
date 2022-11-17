hr_files = {'~/Documents/MATLAB/Approved_Data/HR_cropped/HR_03-18-2022_cropped.txt',...
    '~/Documents/MATLAB/Approved_Data/HR_cropped/HR_04-22-2022_cropped.txt',...
    '~/Documents/MATLAB/Approved_Data/HR_cropped/HR_04-25-2022_cropped.txt',...
    '~/Documents/MATLAB/Approved_Data/HR_cropped/HR_05-06-2022_cropped.txt',...
    '~/Documents/MATLAB/Approved_Data/HR_cropped/HR_05-27-2022_cropped.txt',...
    '~/Documents/MATLAB/Approved_Data/HR_cropped/HR_06-03-2022_cropped.txt',...
    '~/Documents/MATLAB/Approved_Data/HR_cropped/HR_06-17-2022_cropped.txt',...
    '~/Documents/MATLAB/Approved_Data/HR_cropped/HR_06-24-2022_cropped.txt',...
    '~/Documents/MATLAB/Approved_Data/HR_cropped/HR_06-23-2022_cropped.txt',...
    '~/Documents/MATLAB/Approved_Data/HR_cropped/HR_07-07-2022_cropped.txt',...
    '~/Documents/MATLAB/Approved_Data/HR_cropped/HR_08-03-2022_cropped.txt',...
    '~/Documents/MATLAB/Approved_Data/HR_cropped/HR_08-12-2022_part1_cropped.txt',...
    '~/Documents/MATLAB/Approved_Data/HR_cropped/HR_08-12-2022_part2_cropped.txt',...
    '~/Documents/MATLAB/Approved_Data/HR_cropped/HR_09-12-2022_cropped.txt',...
    '~/Documents/MATLAB/Approved_Data/HR_cropped/HR_09-19-2022_cropped.txt',...
    '~/Documents/MATLAB/Approved_Data/HR_cropped/HR_09-23-2022_cropped.txt'}; %The name of the HR file(s) you want to analyze (seperated by commas)

%Affect file
aff_files = {'~/Documents/MATLAB/Approved_Data/Coding/Chat&Chew_2022-03-18_1255_V01_Kessler.csv',...
    '~/Documents/MATLAB/Approved_Data/Coding/Chat&Chew_2022-04-22_1255_V01_Baldie.csv',...
    '~/Documents/MATLAB/Approved_Data/Coding/Art_2022-04-25_OGCP_V01_Montanez.csv',...
    '~/Documents/MATLAB/Approved_Data/Coding/Chat&Chew_2022-05-06_1255_V01_Liu.csv',...
    '~/Documents/MATLAB/Approved_Data/Coding/Chat&Chew_2022-05-27_1255_V01_Baldie.csv',...
    '~/Documents/MATLAB/Approved_Data/Coding/Chat&Chew_2022-06-03_1255_V01_Liu.csv',...
    '~/Documents/MATLAB/Approved_Data/Coding/Chat&Chew_2022-06-17_1255_V01_Kessler.csv',...
    '~/Documents/MATLAB/Approved_Data/Coding/Chat&Chew_2022-06-24_1255_V01_Montanez.csv',...
    '~/Documents/MATLAB/Approved_Data/Coding/Speech_2022-06-23_1255_V01_Montanez.csv',...
    '~/Documents/MATLAB/Approved_Data/Coding/Speech_2022-07-07_1255_V01_Kessler.csv',...
    '~/Documents/MATLAB/Approved_Data/Coding/Speech_2022-08-03_1255_V01_Kessler.csv',...
    '~/Documents/MATLAB/Approved_Data/Coding/Speech_2022-08-12_1429_1255_V01_Montanez_Part1.csv',...
    '~/Documents/MATLAB/Approved_Data/Coding/Speech_2022-08-12_1622_1255_V01_Montanez_Part2.csv',...
    NaN,NaN,NaN};


Data = pshr_load('HR', hr_files, 'Affect', aff_files, 'align', true, 'verbose', false);


% For each applicable affect, run table_combo to get just the face_related
% behavior and omission start and ends for problematic behavior
% aff_list = {'SIB','ISB','inappropriate face related behavior','polar strap adjustment/removal'...
%         'repetitive behaviors','inappropriate movement','crying', 'pulling at pants'};
aff_list_2 = {'SIB','ISB','polar strap adjustment/removal'...
        'repetitive behaviors','inappropriate movement','crying', 'pulling at pants'};

q=1;
store_names = {};
for i = 1:length(Data.HR.Affect)

    if ~isempty(Data.HR.Affect{i})
        keep_table = table_combo(Data.HR.Affect{i}, {'inappropriate face related behavior'});
        omit_table = table_combo(Data.HR.Affect{i}, aff_list_2);
        new_tabs{q} = [keep_table(end-1,:); omit_table(end-1,:)];
        
        % Run the feature_generation function on the HR data
        new_dat{q} = feature_generation(Data.HR.Raw{i}, {[5,0], 'second'}, false);
        
        store_names{q} = Data.HR.files{i};
        q=q+1;
    else
        disp(strcat('No affect found for :', Data.HR.files{i}));
        
    end
end

disp('break');

%% Duration sample
% Sample the HR data using the new affect tables, then run a t-test on the
% differences between the 






%% Onset Sample
% Pass the HR data and new affect tables into onset_sample in order to get
% the samples of the onsets and controls

onset_mat = [];
offset_mat = [];
cont_mat = [];

for i = 1:length(new_tabs)

    % format is (data, union, omission, band, offset)
    if ~isempty(new_tabs{i}{1,2})
        disp(strcat('Running data from :',store_names{i}))
        [on_dump, off_dump, cont_dump] = onset_sample_2(new_dat{i}, new_tabs{i}(1,:), new_tabs{i}(2,:), 'band', [5,0], 'offset', true);
        onset_mat = [onset_mat; on_dump];
        offset_mat = [offset_mat; off_dump];
        cont_mat = [cont_mat; cont_dump];
    end

end


% Run a two sample t-test on the two groups for results, along with
% confidence interval
t_stats=zeros(1,size(onset_mat,2));

for q = 1:size(onset_mat,2)
    [h, p_, ci, stats] = ttest2(onset_mat(:,q),cont_mat(:,q), 'Vartype', 'unequal');
    t_stats(1,q) = p_;
    disp(strcat('two-sample t-test result for feature ', string(q),' = ', string(p_)));
    disp(strcat('MEANS: NON-PROB : ', string(mean(onset_mat(:,q))),' PROB : ', string(mean(cont_mat(:,q)))));
    disp('95% estimated difference:');
    disp(ci);
    disp('Esitmated standard deviation');
    disp(stats.sd);
end










%% Extra functions

function [on_mat, off_mat, un_mat] = onset_sample_2(mat, keep_table, omit_table, varargin)
% This function takes a given matrix of values along with the affect start
% and end index table and returns a matrix consisting of a given selection
% of rows before/after each onset.

% Required Inputs:
%   mat:[n-by-m matrix] the given matrix you are sampling rows from
%   aff_table:[a-by-3 cell array] Found in the Data.(type).Affect structure
%       as an output of load_affect.m, contains the index values of the 
%       start and stop points for each affect
%   aff_list: [1-by-x cell array] list of affects from aff_table to use. If
%       false, then all affects present will be marked.

% Optional Inputs:
%   band: [1-by-2 vector] The number of entries before and after the onset
%       to include in the returned matrix. Default is [0,0] for just the 
%       onset value.
%   offset: [bool] Whether to return a matrix of the ends of each affect.
%       Default is true.
%   omit_nan:[bool] Whether to omit rows which contain NaN values from both
%       on_mat and off_mat. This may result in on_mat and off_mat having
%       different numbers of rows. Default is false.

% Output:
%   on_mat: [o-by-m matrix] matrix of the onsets
%   off_mat: [f-by-m matrix] matrix of the offsets. An empty matrix if the
%       optional input 'offset' is false

    p = inputParser;
    addParameter(p, 'band', [0,0], @ismatrix);
    addParameter(p, 'offset', true, @islogical);
    addParameter(p, 'omit_nan', false, @islogical);
    addParameter(p, 'dilate', 0, @isscalar);
    parse(p,varargin{:});
    
    a = sum(p.Results.band);
    
    
    % ALSO THIS SAMPLING IS SUCCEPTABLE TO VERY SHORT AFFECTS, WHERE THE
    % END POINTS MY HAVE BINS WHICH SPAN THE ENTIRE AFFECT, DO NOT USE
    % OFFSETS WITHOUT THINKING THIS THROUGH
    
    aff_vec = affect_mark(zeros(size(mat,1),1), keep_table, false);
    aff_vec(:,1)=[];
    om_vec = affect_mark(zeros(size(mat,1),1), omit_table, false);
    om_vec(:,1)=[];
    inst = find(aff_vec==1);
    vec = inst(2:end)-inst(1:end-1);
    idx = find(vec > a); % find starts greater in distance then the band being used so you don't
                         % accidently sample from another problematic
                         % behavior
    starts =[inst(1)];
    ends =[];
        
    for k=1:length(idx)
        % Check to make sure that you aren't pulling onsets from areas that
        % overlap with the omit times
        if mean(om_vec(inst(idx(k)+1)-a:inst(idx(k)+1))) == 0 && mean(om_vec(inst(idx(k)):inst(idx(k))+a)) == 0
            starts = [starts,inst(idx(k)+1)];
            ends = [ends,inst(idx(k))];
        end
    end
    ends = [ends, inst(end)]; %grab the last ending
    
    [on_mat, off_mat] = samp_(mat, starts, ends, p.Results.band);
    
    
    if length(keep_table{1,2}) ~= length(starts)
        disp('DISPARITY');
    end
    
    
    
    
    % Combine both the target and omits together and dilate from there to
    % get control measurements
    non_mat = [];
    noff_mat = [];
    
    inst = find(aff_vec + om_vec);
    vec = inst(2:end)-inst(1:end-1);
    idx = find(vec > a); 
    starts =[inst(1)];
    ends =[];
        
    for k=1:length(idx)
        starts = [starts,inst(idx(k)+1)];
        ends = [ends,inst(idx(k))];
    end
    ends = [ends, inst(end)]; %grab the last ending
        
    non_mat = [];
    noff_mat = [];
    
    nstarts = starts;
    nends = ends;
    while size(non_mat,1) < size(on_mat,1)
        
        [nstarts, nends] = dilate(nstarts, nends, p.Results.dilate, a);
        [dump_on, dump_off] = samp_(mat, nstarts, nends, p.Results.band);
        
        non_mat = [non_mat;dump_on];
        noff_mat = [noff_mat;dump_off];
    end
    un_mat = [non_mat; noff_mat];
    clear non_mat noff_mat
    

    if p.Results.omit_nan
        on_mat(any(isnan(on_mat),2),:) = [];
        off_mat(any(isnan(off_mat),2),:) = [];
        un_mat(any(isnan(un_mat),2),:) = [];
    end

end

function [on_mat, off_mat] = samp_(mat, starts, ends, band)
    a = sum(band);
    b = band(1);
    c = band(2);
    
    lim = size(mat,1);
    
    on_mat=[];
    off_mat=[];
    
    for i = 1:length(starts)
        if starts(i) >= b && starts(i)+c <= lim
            on_mat(end+1:end+1+a,:) = mat(starts(i)-b:starts(i)+c,:);
        end
    end
    
    if ends
        for i = 1:length(ends)
            if ends(i)+c <= lim
                off_mat(end+1:end+1+a,:) = mat(ends(i)-b:ends(i)+c,:);
            end
        end
    end

end

function [dstarts, dends] = dilate(starts, ends, amnt, cap)
% Function which dilates the non-zero values of vector vec by the amount
% specified by amnt

    di_vec = [];

    d_starts = max(1,starts-amnt);
    d_ends = ends+amnt;
    
    for i = 1:length(d_starts)
        di_vec(d_starts(i):d_ends(i),1)=1;
    end
    
    
    inst = find(di_vec==1);
    vec = inst(2:end)-inst(1:end-1);
    idx = find(vec > cap); % find starts greater in distance then the band being used so you don't
                         % accidently sample from another problematic
                         % behavior
    dstarts = [inst(1)];
    dends = [];
    for k=1:length(idx)
        dstarts = [dstarts,inst(idx(k)+1)];
        dends = [dends,inst(idx(k))];
    end
    dends =[dends, inst(end)];
    
end






function [mat] = feature_generation(mat, bin, band)
% Function for generating the different features for multiple recording
% sessions

% Inputs:
%   mat: [n-by-m matrix]
%   bin: [1-by-2 cell array] The bin type you want to use
%   band: [1-by-2 matrix]

    mat(:,4) = rmssd_calc(mat(:,3), bin, band);
    mat(:,5) = pnnx_calc(mat(:,3),50, bin, band);
%     mat(:,6) = sdnn_calc(mat(:,3),bin,band);
%     mat(:,7) = sdsd_calc(mat(:,3),bin,band);

end


function [] = t_test_pipelin()
% Pipeline for running t-tests on different aspects of the PHYSMON
% collected data. Effectively any combination of analyses you could want
% will be referenced here.

% Required Inputs:
%   mat: [n-by-m matrix]

% Optional Parameters:
%   ordinal: [bool] whether the data you are providing has ordinal target
%       values

% Returns:
%   mdl: [p-by-4 matrix] matrix containing regression coefficients along
%       with the SE, tStats, and pValue for each coefficient


% Things that this function should do:
%   1. Given a list of files, load them and analyze whatever given affects
%   are specified
%       - Call a T-test function which takes any two matrices and runs a
%       t-test on them and returns/plots the results in a significant
%       manner
%   2. Have the option to compare onsets against control time points



end