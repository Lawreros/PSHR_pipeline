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
hr_files = {'./sample/HR_A.txt',...
    './sample/HR_B.txt'}; %The name of the HR file(s) you want to analyze (seperated by commas)

%ECG file
ecg_files = {'./sample/ECG_A.txt',...
    './sample/ECG_B.txt'}; %The name of the ECG file(s) you want to analyze (seperated by commas)


%Affect file
aff_files = {'./sample/A_coding.csv',...
    './sample/B_coding.csv'};

%% Load and organize data from iPhone files
Data = pshr_load('HR', hr_files, 'ECG', ecg_files, 'Affect', aff_files,'align',true);


%% ECG and RR alignment
[ecg_aligned, aligned_metrics] = ecg_rr_align(Data.HR.Raw{2}(:,[1,3]), Data.ECG.Raw{2}(:,[1,3]), 130, 'verbose', false);


%% Visualization Pipeline
aff_list = {'SIB','ISB','inappropriate face related behavior','polar strap adjustment/removal'...
        'repetitive behaviors','inappropriate movement','crying', 'pulling at pants'};

visualization_pipeline('hr_files', hr_files, 'aff_files', aff_files, 'ecg_files', ecg_files, 'aff_list', aff_list, 'individual_plots',false);


%% Table Combination
% Set up example variables for processing
% aff_1 = {'a';'b';'c'};
% start_1 = {[1,9];[2,7,14,24];[13,20]};
% end_1 = {[4,11];[2,8,17,26];[15,24]};
% tab_1 = [aff_1, start_1, end_1];
%  
% aff_2 = {'d';'e'};
% start_2 = {[1,6,13];[3,18]};
% end_2 = {[2,8,29];[4,19]};
% tab_2 = [aff_2, start_2, end_2];
% 
% aff_3 = {'f';'g'};
% start_3 = {[1,3,7,13,18,23,28];[2,5,10,15,20,25,29]};
% end_3 = {[1,3,7,13,18,23,28];[2,5,10,15,20,25,29]};
% tab_3 = [aff_3, start_3, end_3];
  
%anser = table_combo(tab_1, {'a'},{'b'}, tab_2, {'d'});


%% Affect Marking
%new_mat = affect_mark([], anser, {'a_U_b'}, 'NumberCategories', true);


%% Onset and Offset



%% Train/Test Split



%% Classification Functions



%% Regression Pipeline
%Data = regression_pipeline(hr_files, ecg_files, aff_files, false);


%% Affect Comparison Pipeline


