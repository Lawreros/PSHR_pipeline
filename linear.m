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
% Data = pshr_load('HR', hr_files, 'ECG', ecg_files, 'Affect', aff_files,'align',true);
visualization_pipeline('hr_files', hr_files, 'aff_files', aff_files, 'ecg_files', ecg_files, 'individual_plots',false);

% Data = regression_pipeline(hr_files, ecg_files, aff_files, true);

