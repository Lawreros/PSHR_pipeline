%Coding comparison
%This script is meant to compare the results between coding done by
%different people
clear all;

aff_1 = "ignore_samples/2022-03-25_Montanez.csv";
aff_2 = "ignore_samples/2022-03-25_Kessler.csv";


tab_1 = readtable(aff_1, 'Format', 'auto');
tab_2 = readtable(aff_2, 'Format', 'auto');

% Change column names so you don't have duplicate issues
tab_2.Properties.VariableNames{2} = 'Affect1_2';
tab_2.Properties.VariableNames{3} = 'Affect2_2';
tab_2.Properties.VariableNames{4} = 'Affect3_2';
tab_2.Properties.VariableNames{25} = 'on_camera_2';
tab_2.Properties.VariableNames{26} = 'problem_yn_2';

%Make table for comparison
comp_tab = tab_1(:,[1:4,25,26]);
comp_tab = [comp_tab, tab_2(:,[2:4,25,26])];

%free up some memory
clear tab_1 tab_2;


%Begin comparison iterations



disp('done');