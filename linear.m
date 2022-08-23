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
addpath('./Pipeline');


%% FILENAME INPUT SECTION
% Input the files you wish to analyze

%HR file
hr_files = {'~/Documents/MATLAB/Approved_Data/HR/HR_03-18-2022.txt',...
    '~/Documents/MATLAB/Approved_Data/HR/HR_04-22-2022.txt',...
    '~/Documents/MATLAB/Approved_Data/HR/HR_05-06-2022.txt',...
    '~/Documents/MATLAB/Approved_Data/HR/HR_05-27-2022.txt',...
    '~/Documents/MATLAB/Approved_Data/HR/HR_06-03-2022.txt',...
    '~/Documents/MATLAB/Approved_Data/HR/HR_06-17-2022.txt',...
    '~/Documents/MATLAB/Approved_Data/HR/HR_06-24-2022.txt',...
    '~/Documents/MATLAB/Approved_Data/HR/HR_06-23-2022.txt',...
    '~/Documents/MATLAB/Approved_Data/HR/HR_07-07-2022.txt',...
    '~/Documents/MATLAB/Approved_Data/HR/HR_08-03-2022.txt',...
    '~/Documents/MATLAB/Approved_Data/HR/HR_08-12-2022_part1.txt',...
    '~/Documents/MATLAB/Approved_Data/HR/HR_08-12-2022_part2.txt'}; %The name of the HR file(s) you want to analyze (seperated by commas)

%ECG file
ecg_files = {'~/Documents/MATLAB/Approved_Data/ECG/ECG_03-18-2022.txt',...
    '~/Documents/MATLAB/Approved_Data/ECG/ECG_04-22-2022.txt',...
    '~/Documents/MATLAB/Approved_Data/ECG/ECG_05-06-2022.txt',...
    '~/Documents/MATLAB/Approved_Data/ECG/ECG_05-27-2022.txt',...
    '~/Documents/MATLAB/Approved_Data/ECG/ECG_06-03-2022.txt',...
    '~/Documents/MATLAB/Approved_Data/ECG/ECG_06-17-2022.txt',...
    '~/Documents/MATLAB/Approved_Data/ECG/ECG_06-24-2022.txt',...
    '~/Documents/MATLAB/Approved_Data/ECG/ECG_06-23-2022.txt',...
    '~/Documents/MATLAB/Approved_Data/ECG/ECG_07-07-2022.txt',...
    '~/Documents/MATLAB/Approved_Data/ECG/ECG_08-03-2022.txt',...
    '~/Documents/MATLAB/Approved_Data/ECG/ECG_08-12-2022_part1.txt',...
    '~/Documents/MATLAB/Approved_Data/ECG/ECG_08-12-2022_part2.txt'}; %The name of the ECG file(s) you want to analyze (seperated by commas)

%Affect file
aff_files = {'~/Documents/MATLAB/Approved_Data/Coding/Chat&Chew_2022-03-18_1255_V01_Kessler.csv',...
    '~/Documents/MATLAB/Approved_Data/Coding/Chat&Chew_2022-04-22_1255_V01_Baldie.csv',...
    '~/Documents/MATLAB/Approved_Data/Coding/Chat&Chew_2022-05-06_1255_V01_Liu.csv',...
    '~/Documents/MATLAB/Approved_Data/Coding/Chat&Chew_2022-05-27_1255_V01_Baldie.csv',...
    '~/Documents/MATLAB/Approved_Data/Coding/Chat&Chew_2022-06-03_1255_V01_Liu.csv',...
    '~/Documents/MATLAB/Approved_Data/Coding/Chat&Chew_2022-06-17_1255_V01_Kessler.csv',...
    '~/Documents/MATLAB/Approved_Data/Coding/Chat&Chew_2022-06-24_1255_V01_Montanez.csv',...
    '~/Documents/MATLAB/Approved_Data/Coding/Speech_2022-06-23_1255_V01_Montanez.csv',...
    '~/Documents/MATLAB/Approved_Data/Coding/Speech_2022-07-07_1255_V01_Kessler.csv',...
    '~/Documents/MATLAB/Approved_Data/Coding/Speech_2022-08-03_1255_V01_Kessler.csv',...
    '~/Documents/MATLAB/Approved_Data/Coding/Speech_2022-08-12_1429_1255_V01_Montanez_Part1.csv',...
    '~/Documents/MATLAB/Approved_Data/Coding/Speech_2022-08-12_1622_1255_V01_Montanez_Part2.csv'};

%% Load and organize data from iPhone files
Data = pshr_load('HR', hr_files, 'ECG', ecg_files, 'Affect', aff_files,'align',true);

%visualization_pipeline('hr_files', hr_files, 'ecg_files', ecg_files, 'aff_files', aff_files, 'individual_plots', false);
%Data = regression_pipeline(hr_files, ecg_files, aff_files, true);

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