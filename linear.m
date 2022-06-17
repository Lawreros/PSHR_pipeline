% PHYSMON Pipeline
% The purpose of this file is to make debugging the functions in the GUI
% easier, as checking analysis values in a linear script is much easier
% that the constantly updating/changing GUI.

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
hr_path = '/home/ross/Documents/MATLAB/PSHR_pipeline/sample/';
hr_file = {'HR_A.txt'};

%ECG file
ecg_path = '/home/ross/Documents/MATLAB/PSHR_pipeline/sample/';
ecg_file = {'ECG_A.txt'};

%Affect file
aff_path = "/home/ross/Documents/MATLAB/PSHR_pipeline/sample/";
aff_file = "A_coding.csv";

%realtime = "11:19:15"
%videotime = 728

%Export RR-intervals with affects denoted for Richard
aff = {'clothing adjustment/removal','flapping/clapping','loud/rapid humming',...
    'loud/rapid speech','moving at a fast/abrupt pace','polar strap adjustment/removal',...
    'repetitive behaviors','unresponsive/unable to redirect'};


%% Analysis flags
% As the GUI will most likely chain together preprocessing modules through
% the use of binary triggers, it's good to simulate that with a group of
% True/False statements.

% TODO: Create structure to store flags in compact manner

%RR-interval Preprocess Flags
Bandpass = false;
u_band = 1200;
l_band = 400;

Malik = false;
Kamath = false;
Karlsson = false;

Acar = false;
acar_range = 9;


%RR-interval Analysis Flags
sdsd = false;
sdnn = false;
rmssd = false;
pnnx = false;
poincare = false;

%RR-interval Plotting Flags
raw_hr_plot = false;

%ECG-interval Preprocess Flags
ecg2rr = false;
peak = 800;
dist = 40;

%ECG-interval Analysis Flags


%ECG-interval Plotting Flags



%% Pipeline
Data.HR.Raw{1} = {};
Data.ECG.Raw{1} = {};
Data.Affect.Raw{1} = {};

Data = pshr_load_data(Data, hr_path, hr_file, "HR");
Data = pshr_load_data(Data, ecg_path, ecg_file, "ECG");
Data = LoadAffect(Data, aff_path, aff_file);

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

cuts = {[1,1462],[1,1920],false,[1,1060],false,false,false,false,[1,3180],[1,4900],false};
Data = group_analysis(Data, 'HR', false, cuts);
Data = group_analysis(Data, 'HR', {5, 'second'},cuts);

%Derek_export(Data, aff, "HR", "Derek_HR_output");
%Richard_export(Data, aff, "HR", "test_HR_output");

%% ECG-Preprocessing

    
%% ECG Exports

%Richard_export(Data, aff, "ECG", "test_ECG_output");
disp('done');


%% Affect Loading

function [Data] = LoadAffect(Data, path, file)
% Load in affect file and add it to the structure
Data.Affect.path = path;

%TODO: Add functionality for loading and processing multiple affects at a
%time

    if iscell(file)
        for i = 1:length(file)
            %For some reason you have to set Format to auto in order for
            %readtable to not ignore affects in sparsely filled columns
            %(i.e. if there is only like 4 entries in Affect2, without
            %auto, MATLAB will disregard Affect2 and have the column be all
            %NaNs)
            Data.Affect.Raw{i} = readtable(strcat(path,file{i}),'Format','auto');
        end
        Data.Affect.files = {file};
    else
        Data.Affect.Raw = readtable(strcat(path,file),'Format','auto');
        Data.Affect.files = {file};
    end

    
    %Get list of all unique affects used in the coding (SINGLE)
    aff_list = unique(Data.Affect.Raw.Affect1);
    Data.Affect.Raw.Affect1{1} = "start";
    Data.Affect.Raw.Affect1{end+1} = "end";
    
    if iscell(unique(Data.Affect.Raw.Affect2(1:end-1)))
        ext = unique(Data.Affect.Raw.Affect2(~cellfun(@isempty, Data.Affect.Raw.Affect2)));
%         aff_list = [aff_list; unique(Data.Affect.Raw.Affect2)];
        aff_list = [aff_list; ext];
        Data.Affect.Raw.Affect2{1} = "start";
        Data.Affect.Raw.Affect2{end} = "end";
    else
        disp('No entries in column Affect2');
    end
    
    if iscell(unique(Data.Affect.Raw.Affect3(1:end-1)))
        ext = unique(Data.Affect.Raw.Affect3(~cellfun(@isempty, Data.Affect.Raw.Affect3)));
%         aff_list = [aff_list; unique(Data.Affect.Raw.Affect3)];
        Data.Affect.Raw.Affect3{1} = "start";
        Data.Affect.Raw.Affect3{end} = "end";
    else
        disp('No entries in column Affect3');
    end
    
    aff_list = unique(aff_list); %cell array of all affects used
    
    
    %Generate Start and End times for the Affects
    Data.Affect.Times = {};
    for i = 1:length(aff_list)
        starts = [];
        ends = [];
        for k = 1:3
            col = strcat("Affect",string(k));
        
            buffer = [false, transpose(diff(strcmp(Data.Affect.Raw.(col), aff_list{i}))~=0)];
            buffer = find(buffer);
            
        
            for j = 1:2:length(buffer)
                starts = [starts, Data.Affect.Raw.Time_sec(buffer(j))];%buffer(j)];
                ends = [ends, Data.Affect.Raw.Time_sec(buffer(j+1)-1)];%buffer(j+1)-1];
            end
        end
        Data.Affect.Times{i,1} = aff_list{i};
        Data.Affect.Times{i,2} = starts;
        Data.Affect.Times{i,3} = ends;
    end
    
    %Load and store alignment time
    
    %realtime = "11:19:15"
    %videotime = 728 -1 for lag
    
    pol_time = (((((11*60)+19)*60)+15)*1000);
    vid_time = 727*1000;
    algn = pol_time - vid_time;
    
    %Check if there is any HR or ECG data loaded. If so, then generate the
    %index numbers for the start and stop.
    
    if iscell(Data.HR.Raw) == 1
        disp("HR data found, generating start and stop indexes");
        Data = time_adjust(Data, algn, "HR");
    end
        
    if iscell(Data.ECG.Raw) == 1
        disp("ECG data found, generating start and stop indexes"); 
        Data = time_adjust(Data, algn, "ECG");
    end


end

function [Data] = time_adjust(Data, algn, type)
%Changes the generates indexes for Data.*.Raw that correspond with
%timestamps found in Data.Affect.Times

%inputs:
% Data: main Data structure
% algn: alignment time found with (polar_timestamp - video_time) = corr
% type: What type of Data you are adjusting the time of


    % Iterate through start and stop times
    [r,c] = size(Data.Affect.Times);
    
    % Create place to store new times
    Data.(type).Affect = {};

    
    %Iterate through each affect
    for i = 1:r
        starts = [];
        ends = [];
        
        %Iterate through each start/stop pair
        for j = 1:length(Data.Affect.Times{i,2})
            
            a = find(Data.(type).Raw(:,1) >= (Data.Affect.Times{i,2}(j)*1000+algn));
            %disp(Data.Affect.Times{i,2}(j)*1000+algn);
            
            b = find(Data.(type).Raw(:,1) > (Data.Affect.Times{i,3}(j)*1000+algn));
            %disp(Data.Affect.Times{i,3}(j)*1000+algn);
            
            if isempty(a)==0 && isempty(b)==0
                if a(1) == b(1)
                        %nothing is there
                else
                        starts = [starts, a(1)];
                        ends = [ends, b(1)-1];
                end
            elseif isempty(a)==0 && isempty(b)==1
                fprintf("Affect %s ends after the recording and starts at time %d\n",...
                    Data.Affect.Times{i,1}, Data.Affect.Times{i,2}(j));
                    starts = [starts, a(1)];
                    ends = [ends, length(Data.(type).Raw)];
            else
                %Issue due to the affects occuring outside of collected
                %HR/ECG data
                fprintf('WARNING: Affect %s from video time %d to %d could not be found in data\n',...
                    Data.Affect.Times{i,1}, Data.Affect.Times{i,2}(j), Data.Affect.Times{i,3}(j));
                
            end
        end
            % Store indexes for each affect with their respective data
            % types, instead of keeping everything in the Affect structure.
            % This way, every thing you need for each analysis will be
            % grouped
            Data.(type).Affect{i,1} = Data.Affect.Times{i,1};
            Data.(type).Affect{i,2} = starts;
            Data.(type).Affect{i,3} = ends;
    end


end




%% Export Functions
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