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
    '~/Documents/MATLAB/Approved_Data/HR_cropped/HR_09-23-2022_cropped.txt',...
    '~/Documents/MATLAB/Approved_Data/HR_cropped/HR_09-26-2022_cropped.txt',...
    '~/Documents/MATLAB/Approved_Data/HR_cropped/HR_09-30-2022_cropped.txt',...
    '~/Documents/MATLAB/Approved_Data/HR_cropped/HR_10-03-2022_cropped.txt',...
    '~/Documents/MATLAB/Approved_Data/HR_cropped/HR_10-17-2022_cropped.txt',...
    '~/Documents/MATLAB/Approved_Data/HR_cropped/HR_10-24-2022_cropped.txt',...
    '~/Documents/MATLAB/Approved_Data/HR_cropped/HR_10-31-2022_cropped.txt',...
    '~/Documents/MATLAB/Approved_Data/HR_cropped/HR_11-07-2022_cropped.txt'}; %The name of the HR file(s) you want to analyze (seperated by commas)
% 
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
    '~/Documents/MATLAB/Approved_Data/Coding/Art_2022-09-12_OGCP_V01_Montanez.csv',...
    '~/Documents/MATLAB/Approved_Data/Coding/Art_2022-09-19_OGCP_V01_Kessler.csv',...
    '~/Documents/MATLAB/Approved_Data/Coding/Chat&Chew_2022-09-23_1255_V01_Montanez.csv',...
    '~/Documents/MATLAB/Approved_Data/Coding/Art_2022-09-26_OGCP_V01_Montanez.csv',...
    '~/Documents/MATLAB/Approved_Data/Coding/Chat&Chew_2022-09-30_1255_V01_Kessler.csv',...
    '~/Documents/MATLAB/Approved_Data/Coding/Art_2022-10-03_OGCP_V01_Kessler.csv',...
    '~/Documents/MATLAB/Approved_Data/Coding/Art_2022-10-17_OGCP_V01_Montanez.csv',...
    '~/Documents/MATLAB/Approved_Data/Coding/Art_2022-10-24_OGCP_V01_Kessler.csv',...
    '~/Documents/MATLAB/Approved_Data/Coding/Art_2022-10-31_OGCP_V01_Montanez.csv',...
    '~/Documents/MATLAB/Approved_Data/Coding/Art_2022-11-07_OGCP_V01_Kessler.csv'};

% hr_files = {'~/Documents/MATLAB/Approved_Data/HR_cropped/HR_09-19-2022_cropped.txt',...
%     '~/Documents/MATLAB/Approved_Data/HR_cropped/HR_09-26-2022_cropped.txt',...
%     '~/Documents/MATLAB/Approved_Data/HR/HR_10-31-2022.txt',...
%     '~/Documents/MATLAB/Approved_Data/HR/HR_11-07-2022.txt'}; %The name of the HR file(s) you want to analyze (seperated by commas)
% 
% aff_files = {'~/Documents/MATLAB/Approved_Data/Timer_coding/Art_2022-09-19_OGCP_V01_Kessler.csv',...
%     '~/Documents/MATLAB/Approved_Data/Timer_coding/Art_2022-09-26_OGCP_V01_Montanez.csv',...
%     '~/Documents/MATLAB/Approved_Data/Coding/Art_2022-10-31_OGCP_V01_Montanez.csv',...
%     '~/Documents/MATLAB/Approved_Data/Timer_coding/Art_2022-11-07_OGCP_V01_Kessler.csv'};




% For each applicable affect, run table_combo to get just the face_related
% behavior and omission start and ends for problematic behavior
 aff_list = {{'SIB','ISB','inappropriate face related behavior','polar strap adjustment/removal'...
         'repetitive behaviors','inappropriate movement','crying', 'pulling at pants'}};
%aff_list_2 = {{'SIB','ISB','polar strap adjustment/removal',...
%        'inappropriate movement','crying', 'pulling at pants','off camera'}};

target = {{'inappropriate face related behavior'}};
    
%t_test_pipelin(hr_files, aff_files, {{'Timer_Used'}}, aff_list, {[5,0], 'second'}, true, true);
%t_test_pipelin(hr_files, aff_files, target, aff_list_2, {[5,0], 'second'}, true, true);

t_test_pipelin(hr_files, aff_files, aff_list, {}, {[5,0],'second'}, false, true);


function [] = t_test_pipelin(hr_files, aff_files, target, omit, bin, duration, onset)
% Pipeline for running t-tests on different aspects of the PHYSMON
% collected data. Effectively any combination of analyses you could want
% will be referenced here.

% Required Inputs:
%   hr_files: [1-by-n cell array] cell array containing the HR files that
%       you wish to run the t-test on
%   aff_files: [1-by-n cell array] cell array containing the affect files
%       which you wish to use with the hr_files you provided (in the same
%       order).
%   target: [1-by-n cell array] cell array containing the targets as {{'A','B'},{'C'}}
%   omit: []
%   duration: [bool] Whether to run duration analysis
%   onset: [bool] Whether to run onset analysis

% Optional Parameters:

% Returns:
%   

    Data = pshr_load('HR', hr_files, 'Affect', aff_files, 'align', true, 'verbose', false);
    
     for i = 1:length(hr_files)
        Data.HR.PP{i} = Data.HR.Raw{i};
        % We'll just work with bandpassing for now...
        Data.HR.PP{i}(:,3) = bandpass(Data.HR.PP{i}(:,3), 300, 1600, false);
    end
    
    
    

    if duration
        %% Duration sample
        % Sample the HR data using the new affect tables, then run a t-test on the
        % differences between the 

        q=1;
        store_names = {};
        for i = 1:length(Data.HR.Affect)

            if ~isempty(Data.HR.Affect{i})
                new_tabs{q} = table_combo(Data.HR.Affect{i}, target{:}, 'omit', omit{:});

                % Run the feature_generation function on the HR data
                new_dat{q} = feature_generation(Data.HR.PP{i}, bin, false);

                store_names{q} = Data.HR.files{i};
                q=q+1;
            else
                disp(strcat('No affect found for :', Data.HR.files{i}));

            end
        end


        % Select data using starts and stops
        dur_mat = [];
        dur_dump = [];
        cont_mat = [];
        cont_dump = [];
        a = 10;
        for i = 1:length(new_tabs)

            if ~isempty(new_tabs{i}{1,2})
                disp(strcat('Running data from :',store_names{i}))

                % Select data from outside both the omission and target datasets

                aff_vec = affect_mark(zeros(size(new_dat{i},1),1), new_tabs{i}(1,:), false);
                aff_vec(:,1)=[];
                om_vec = affect_mark(zeros(size(new_dat{i},1),1), new_tabs{i}(2,:), false);
                om_vec(:,1)=[];

                inst = find((aff_vec + om_vec)==0);
                vec = inst(2:end)-inst(1:end-1);
                idx = find(vec > a);
                starts =[inst(1)];
                ends =[];

                for k=1:length(idx)
                    starts = [starts,inst(idx(k)+1)];
                    ends = [ends,inst(idx(k))];
                end
                ends = [ends, inst(end)]; %grab the last ending

                % Add buffer around all of the non-control starts/stops in order to
                % select control timepoints more than X away from affects
                [starts, ends] = dilate(starts, ends, -10, a);

                 % Select all target data
                for j = 1:length(new_tabs{i}{1,2})
                    dur_dump = [dur_dump; new_dat{i}(new_tabs{i}{1,2}(j):new_tabs{i}{1,3}(j),:)];
                end

                % Select all control data
                for j = 1:length(starts)
                    cont_dump = [cont_dump;new_dat{i}(starts(j):ends(j),:)];
                end

                dur_mat = [dur_mat; dur_dump];
                cont_mat = [cont_mat; cont_dump];
            end

        end

        % Run a two sample t-test on the two groups for results, along with
        % confidence interval
        t_stats=zeros(1,size(dur_mat,2));

        if isempty(dur_mat)
            disp('No instances of target event(s) occuring. This may be due to misspelling of the target(s).');
        else
            for q = 1:size(dur_mat,2)
                [h, p_, ci, stats] = ttest2(dur_mat(:,q),cont_mat(:,q), 'Vartype', 'unequal');
                t_stats(1,q) = p_;
                disp(strcat('two-sample t-test result for feature ', string(q),' = ', string(p_)));
                disp(strcat('MEANS: target : ', string(nanmean(dur_mat(:,q))),' control : ', string(nanmean(cont_mat(:,q)))));
                disp('95% estimated difference:');
                disp(ci);
                disp('Esitmated standard deviation');
                disp(stats.sd);
            end
        end
        
        i=1;
        while i < 500
    %         [Data.HR.RegPP(i,:,:), Data.HR.PPhat(i,:,:), Data.HR.Ttest(i,:)] = gen_regression(big(:,[3:end-1]),big(:,end), 'log');
            [Data.HR.RegPP(i,:,:), Data.HR.PPhat(i,:,:), Data.HR.Ttest(i,:), Data.HR.AUC(i,1)] = gen_regression([dur_mat(:,3:5);cont_mat(:,3:5)],[zeros(size(dur_mat,1),1); ones(size(cont_mat,1),1)], ...
                'log', 'verbose',false);
            i = i+1;
        end

        % Plot the histograms for the beta values and the confidence intervals
        row = floor(sqrt(size(Data.HR.RegPP,2)));
        col = ceil(size(Data.HR.RegPP,2)/row);
        for i =1:size(Data.HR.RegPP,2)
            subplot(row,col,i)
            histogram(Data.HR.RegPP(:,i,5));
            hold on;
            histogram(Data.HR.RegPP(:,i,6));
            histogram(Data.HR.RegPP(:,i,1));
            legend('upper bound', 'lower bound', 'beta');
            %boxchart([Data.HR.RegPP(:,i,5),Data.HR.RegPP(:,i,6)]);
        end
        title('Beta values and 95% confidence intervals');

        % Plot p-values for the coefficents
        figure;
        for i = 1:size(Data.HR.RegPP,2)
            subplot(row,col,i)
            histogram(Data.HR.RegPP(:,i,4));
        end
        title('coefficient p-values');


        % Plot the accuracy results
        figure;
        for i = 1:2
            subplot(1,2,i)
            histogram(Data.HR.PPhat(:,1,i));
            hold on;
            histogram(Data.HR.PPhat(:,2,i));
            histogram(Data.HR.PPhat(:,3,i));
            hold off;
            legend('Lower Bound', 'Pred', 'Upper Bound');
        end
        title('Accuracy results');

        % Plot spread of AUC
        figure;
        histogram(Data.HR.AUC);
        title('AUC Results')
        
        
        

    end

    if onset
        %% Onset Sample
        % Pass the HR data and new affect tables into onset_sample in order to get
        % the samples of the onsets and controls

        q=1;
        store_names = {};
        new_tabs = {};
        for i = 1:length(Data.HR.Affect)

            if ~isempty(Data.HR.Affect{i})
                keep_table = table_combo(Data.HR.Affect{i}, target{:});
                
                if ~isempty(omit)
                    omit_table = table_combo(Data.HR.Affect{i}, omit{:});
                    new_tabs{q} = [keep_table(end-1,:); omit_table(end-1,:)];
                else %If they aren't omitting anything, then just use keep_table
                    new_tabs{q} = keep_table;
                end

                % Run the feature_generation function on the HR data
                new_dat{q} = feature_generation(Data.HR.PP{i}, bin, false);

                store_names{q} = Data.HR.files{i};
                q=q+1;
            else
                disp(strcat('No affect found for :', Data.HR.files{i}));

            end
        end



        onset_mat = [];
        offset_mat = [];
        cont_mat = [];

        for i = 1:length(new_tabs)

            % format is (data, union, omission, band, offset)
            if ~isempty(new_tabs{i}{1,2})
                disp(strcat('Running data from :',store_names{i}))
                [on_dump, off_dump, cont_dump] = onset_sample_2(new_dat{i}, new_tabs{i}(1,:), new_tabs{i}(2,:), 'band', [3,0], 'dilate', 6);
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
            disp(strcat('MEANS: TARGET : ', string(nanmean(onset_mat(:,q))),' CONTROL : ', string(nanmean(cont_mat(:,q)))));
            disp('95% estimated difference:');
            disp(ci);
            disp('Esitmated standard deviation');
            disp(stats.sd);
        end
        
        i=1;
        while i < 500
    %         [Data.HR.RegPP(i,:,:), Data.HR.PPhat(i,:,:), Data.HR.Ttest(i,:)] = gen_regression(big(:,[3:end-1]),big(:,end), 'log');
            [Data.HR.RegPP(i,:,:), Data.HR.PPhat(i,:,:), Data.HR.Ttest(i,:), Data.HR.AUC(i,1)] = gen_regression([onset_mat(:,3:5);cont_mat(:,3:5)],[zeros(size(onset_mat,1),1); ones(size(cont_mat,1),1)], ...
                'log', 'verbose',false);
            i = i+1;
        end

        % Plot the histograms for the beta values and the confidence intervals
        row = floor(sqrt(size(Data.HR.RegPP,2)));
        col = ceil(size(Data.HR.RegPP,2)/row);
        for i =1:size(Data.HR.RegPP,2)
            subplot(row,col,i)
            histogram(Data.HR.RegPP(:,i,5));
            hold on;
            histogram(Data.HR.RegPP(:,i,6));
            histogram(Data.HR.RegPP(:,i,1));
            legend('upper bound', 'lower bound', 'beta');
            %boxchart([Data.HR.RegPP(:,i,5),Data.HR.RegPP(:,i,6)]);
        end
        title('Beta values and 95% confidence intervals');

        % Plot p-values for the coefficents
        figure;
        for i = 1:size(Data.HR.RegPP,2)
            subplot(row,col,i)
            histogram(Data.HR.RegPP(:,i,4));
        end
        title('coefficient p-values');


        % Plot the accuracy results
        figure;
        for i = 1:2
            subplot(1,2,i)
            histogram(Data.HR.PPhat(:,1,i));
            hold on;
            histogram(Data.HR.PPhat(:,2,i));
            histogram(Data.HR.PPhat(:,3,i));
            hold off;
            legend('Lower Bound', 'Pred', 'Upper Bound');
        end
        title('Accuracy results');

        % Plot spread of AUC
        figure;
        histogram(Data.HR.AUC);
        title('AUC Results')
        
    end

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
        if starts(i) > b && starts(i)+c <= lim
            on_mat(end+1:end+1+a,:) = mat(starts(i)-b:starts(i)+c,:);
        end
    end
    
    if ends
        for i = 1:length(ends)
            if ends(i)+c <= lim && ends(i)-b > 0
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