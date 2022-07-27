% PHYSMON Pipeline
% The purpose of this file is to make debugging the functions in the GUI
% easier, as checking analysis values in a linear script is much easier
% that the constantly updating/changing GUI

clear all;
close all;

% Add location of supporting functions to path for use. Doing this without
% saving the new path makes this addition only exist for the current MATLAB
% session

addpath('./Analysis');
addpath('./Preprocess');
addpath('./Export');
addpath('./Import');


%% FILENAME INPUT SECTION
% Input the files you wish to analyze

%HR file
hr_path = './group_HR_analysis/'; %The location of the directory containing the HR file you want to analyze
hr_file = {'HR_05-13-2022.txt','HR_05-09-2022.txt','HR_05-16-2022.txt'}; %The name of the HR file(s) you want to analyze (seperated by commas)

%ECG file
ecg_path = './group_HR_analysis/'; %The location of the directory containing the ECG file you want to analyze
ecg_file = {'ECG_05-09-2022.txt'}; %The name of the ECG file(s) you want to analyze (seperated by commas)

%Affect file
aff_path = "./sample/";
aff_file = {'A_coding.csv'};


%% Load and organize data from iPhone files
Data.HR.Raw{1} = {};
Data.ECG.Raw{1} = {};
Data.Affect.Raw{1} = {};

Data = pshr_load_data(Data, hr_path, hr_file, "HR");
Data = pshr_load_data(Data, ecg_path, ecg_file, "ECG");
%Data = load_affect(Data, aff_path, aff_file);


%% ECG preprocessing

%known issues with ecg data:
%   1. Inverted ECG measurements (lead placement results in upside-down
%   QRS)
%       Solution: The mean uV is not always positive or negative based off
%       of inverted data





%   2. Maxing out of signal/discrete shift
%       Solution: Index where the measured value falls outside of -5000 to 5000
%       then take x indicies before and after it and remove those. This
%       will still work with findpeaks and takes care of the issue with
%       adjustment noise.
mat = Data.ECG.Raw{1}(:,3);
amp = 5000; %Maximum amplitude
cut_bin = 100;
[ret, locs] = ecg_preprocess(mat, amp, cut_bin);

%% Find P,Q,R,S,T

Data.ECG.PP = {};
Data.ECG.PP{1} = Data.ECG.Raw{1};
Data.ECG.PP{1}(:,3) = ret;

%TODO: Go though default settings to find what has the minimum amount of
%NaNs
samp = ecg_PQRST(Data.ECG.PP{1}(:,3));

%% Go through combinations of parameters for best RR/ECG alignment

[a,b,c] = ecg_rr_alignment(Data.HR.Raw{1}(:,[1,3]), Data.ECG.Raw{1}(:,[1,3]),700,50,130,10,true);

val_iter_results = [NaN, NaN, NaN, NaN, NaN, NaN];
time_iter_results = [NaN, NaN, NaN, NaN, NaN, NaN];
%peak
for i=600:50:900
    % dist
    for j=40:10:70
    % subcost
        for k= 5:5:20
            [a,b,c] = ecg_rr_alignment(Data.HR.Raw{1}(:,[1,3]), Data.ECG.Raw{1}(:,[1,3]),i,j,130,k,false);
            val_iter_results(end+1,:) = [sum(~isnan(c.val.diff))/length(c.val.diff), c.val.mean, c.val.std, i,j,k];
            time_iter_results(end+1,:) = [sum(~isnan(c.time.diff))/length(c.time.diff), c.time.mean, c.time.std, i,j,k];
        end
    end
end




%% Plot RR-Interval and ECG data

figure(1);
plot(Data.HR.Raw{1}(:,3));
title("RR-interval Data");
xlabel("Index");
ylabel("Duration (millisecond)");


figure(2);
plot(Data.ECG.Raw{1}(:,3));
title("ECG Data");
xlabel("Index");
ylabel("Voltage");



%% RR-Interval Analysis
% This is where files are exported and saved
Data = group_analysis(Data, "HR", {5,"second"},false);



%% Export Functions

%Export RR-intervals with affects denoted for Richard
aff = {'clothing adjustment/removal','flapping/clapping','loud/rapid humming',...
    'loud/rapid speech','moving at a fast/abrupt pace','polar strap adjustment/removal',...
    'repetitive behaviors','unresponsive/unable to redirect'};


function [] = Richard_export(Data,aff,type,fil_name)
%This function takes the Data information and saves them as .csv files for
%sending to Richard

%input:
%   Data: The data structure
%   aff: Cell array listing the different Affects that should be denoted as
%   problem behaviors
%   type: What type of data you want to export, "HR" or "ECG"
%   fil_name: What you want the created csv file to be called. This will be
%   saved locally unless full path is specified in file name


    %TODO: add functionality for multiple sets of data
    new_mat = Data.(type).Raw(:,[1,3]);
    new_mat(1,3) = 0;

    
    for i = 1:length(Data.(type).Affect)
        
        %only put a 1 in the third column for affects in the list
        if any(strcmp(aff, Data.(type).Affect(i,1))) %if the affect is in the list of aff
            if isempty(Data.(type).Affect(i,2))==1 || isempty(Data.(type).Affect(i,3))==1
                %There are no/missing start and stop indexes for this
                %affect, so do nothing
            else
                %Mark the 3rd column on the new_mat with a 1, denoting that
                %an affect of interest occurs there
                for j = 1:length(Data.(type).Affect{i,2})
                    %this should also take care of the issue of overlapping
                    %affects
                    new_mat(Data.(type).Affect{i,2}(j):Data.(type).Affect{i,3}(j),3) = 1;
                    
                end
            end
        end
    end
    
    writematrix(new_mat,strcat(fil_name,"_List.csv"));
end

function [] = Derek_export(Data,aff,type,fil_name)

% RR interval data with: [Timestamp, RR, pnn50, rmssd, sdnn, problem]


%This function takes the Data information and saves them as .csv files for
%sending to Richard

%input:
%   Data: The data structure
%   aff: Cell array listing the different Affects that should be denoted as
%   problem behaviors
%   type: What type of data you want to export, "HR" or "ECG"
%   fil_name: What you want the created csv file to be called. This will be
%   saved locally unless full path is specified in file name


    %TODO: add functionality for multiple sets of data
    new_mat = Data.(type).Raw(:,[1,3]);
    
    new_mat(:,end+1) = pnnx_calc(new_mat(:,2),50,{10,'second'},false);
    new_mat(:,end+1) = rmssd_calc(new_mat(:,2),{10,'second'},false);
    new_mat(:,end+1) = sdnn_calc(new_mat(:,2),{10,'second'},false);
    new_mat(:,end+1) = sdsd_calc(new_mat(:,2),{10,'second'},false);
     
    new_mat(:,end+1) = zeros(length(new_mat),1);

    
    for i = 1:length(Data.(type).Affect)
        
        %only put a 1 in the third column for affects in the list
        if any(strcmp(aff, Data.(type).Affect(i,1))) %if the affect is in the list of aff
            if isempty(Data.(type).Affect(i,2))==1 || isempty(Data.(type).Affect(i,3))==1
                %There are no/missing start and stop indexes for this
                %affect, so do nothing
            else
                %Mark the 3rd column on the new_mat with a 1, denoting that
                %an affect of interest occurs there
                for j = 1:length(Data.(type).Affect{i,2})
                    %this should also take care of the issue of overlapping
                    %affects
                    new_mat(Data.(type).Affect{i,2}(j):Data.(type).Affect{i,3}(j),end) = 1;
                    
                end
            end
        end
    end
    
    writematrix(new_mat,strcat(fil_name,"_List.csv"));

end

function [Data] = group_analysis(Data, type, bin, band)
% RR-interval analysis for looking at statistics of multiple recording
% sesisons

%inputs:
%   Data: The data structure
%   type: [string] What type of data you want to analyze, "HR" or "ECG"
%   bin: [1-by-2 cell array] The bin type you want to use. If false, no
%   binning is done

    if iscell(bin)
        rname = strcat('RMSSD_',string(bin{1}));
        pname = strcat('pnnx_',string(bin{1}));
    else
        rname = strcat('RMSSD_','nbin');
        pname = strcat('pnnx_','nbin');
    end
    
    if iscell(band)
        rname = strcat(rname, '_band');
        pname = strcat(pname, '_band');
    else
        band = cell(1,length(Data.(type).Raw));
        band(:) = {false};
    end
    
    Data.(type).(rname){1}={};
    Data.(type).(pname){1}={};
    
    rst_tab = strcat(rname,'_stats');
    Data.(type).(rst_tab){1} = {};
    pst_tab = strcat(pname,'_stats');
    Data.(type).(pst_tab){1} = {};
    
    

    for i = 1:length(Data.(type).Raw)
        Data.(type).(rname){i} = rmssd_calc(Data.(type).Raw{i}(:,3),bin,band{i});
        Data.(type).(pname){i} = pnnx_calc(Data.(type).Raw{i}(:,3),50,bin,band{i});
        
        % Record mean and standard deviation for dataset
        index = Data.(type).(rname){i}==Inf;
        
        Data.(type).(rst_tab){i,1} = mean(Data.(type).(rname){i}(index==0));
        Data.(type).(rst_tab){i,2} = std(Data.(type).(rname){i}(index==0));
        
        index = Data.(type).(pname){i}==Inf;
        
        Data.(type).(pst_tab){i,1} = mean(Data.(type).(pname){i}(index==0));
        Data.(type).(pst_tab){i,2} = std(Data.(type).(pname){i}(index==0));
    end
    
    

    disp('done');
end

function [ret,locs] = ecg_preprocess(mat, amp, cut_bin)
% Function to take the ECG data preprocesses it by removing ECG values
% which fall outside of the accepted amplitude
    %inputs:
    %   mat: [n-by-1] vector containing ecg data
    %   amp: [int], minimum amplitude of accepted ECG values. If a value
    %   falls outside of this bound, all values [cut_bin] before and after
    %   it are replaced with NaNs.
    %   cut_bin: [int], the amount of indexes before and after entries
    %   which fail [amp] that are replaced with NaNs.
    
    %Returns:
    %   ret: [n-by-1] matrix containing the ecg data with all removed
    %   values replaced by NaNs
    %   locs: [m-by-1] index of all values which fall outside of the bounds
    %   described by [amp]

    locs = find(abs(mat)>amp);
    ret = mat;
    max_len = length(ret);
    
    for i = 1:length(locs)
    
        if locs(i) <= cut_bin
            ret(1:locs(i)+cut_bin) = NaN;
        elseif locs(i)+cut_bin >= max_len
            ret(locs(i)-cut_bin:end) = NaN;
        else
            ret(locs(i)-cut_bin:locs(i)+cut_bin) = NaN;
        end
    end
end

