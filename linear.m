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


%% Test files:

%HR file
hr_path = './sample/';
hr_file = {'HR_A.txt'};

%ECG file
ecg_path = './sample/';
ecg_file = {'ECG_A.txt'};

%Affect file
aff_path = "./sample/";
aff_file = {'A_coding.csv'};

%realtime = "11:19:15"
%videotime = 728

%% Analysis flags
% As the GUI will most likely chain together preprocessing modules through
% the use of binary triggers, it's good to simulate that with a group of
% True/False statements.


%% Pipeline
Data.HR.Raw{1} = {};
Data.ECG.Raw{1} = {};
Data.Affect.Raw{1} = {};

Data = pshr_load_data(Data, hr_path, hr_file, "HR");
Data = pshr_load_data(Data, ecg_path, ecg_file, "ECG");
Data = load_affect(Data, aff_path, aff_file);

%% RR-Interval Preprocessing


%% Simple Plot of raw RR data

if raw_hr_plot
    fig1 = figure(1);
    plot(Data.HR.Raw(:,3))
    title("Raw RR-interval Data");
end

if poincare
    fig2 = figure(2);
    [SD1, SD2] = poincare_plot(Data.HR.Raw(:,3), fig2);
end
%% RR-Interval Exports
% This is where files are exported and saved
Data = group_analysis(Data, "HR", false,false);

% cuts = {[1,1462],[1,1920],false,[1,1060],false,false,false,false,[1,3180],[1,4900],false};
% Data = group_analysis(Data, 'HR', false, cuts);
% Data = group_analysis(Data, 'HR', {5, 'second'},cuts);



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