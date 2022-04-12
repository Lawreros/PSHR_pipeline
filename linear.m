% PHYSMON Pipeline
% The purpose of this file is to make debugging the functions in the GUI
% easier, as checking analysis values in a linear script is much easier
% that the constantly updating/changing GUI.

clear all;
close all;

%% Test files:

hr_path = '/home/ross/Downloads/';%"/home/ross/Documents/MATLAB/PSHR_pipeline/sample/";
hr_file = 'HR_03-18-2022.txt';%"HR_03-09-2022.txt";
ecg_path = '/home/ross/Downloads/';%"/home/ross/Documents/MATLAB/PSHR_pipeline/sample/";
ecg_file = 'ECG_03-18-2022.txt';%"ECG_03-09-2022.txt";

aff_path = "/home/ross/Documents/MATLAB/PSHR_pipeline/sample/";
%aff_file = "03-18-code.csv";
aff_file = "2022-03-18_Kessler.csv";

%realtime = "11:19:15"
%videotime = 728

%% Analysis flags
% As the GUI will most likely chain together preprocessing modules through
% the use of binary triggers, it's good to simulate that with a group of
% True/False statements.


%RR-interval
Bandpass = false;
u_band = 1200;
l_band = 400;

Malik = false;
Kamath = false;
Karlsson = false;

Acar = false;
acar_range = 9;


%ECG-interval



%% Pipeline
Data.HR.Raw{1} = {};
Data.ECG.Raw{1} = {};
Data.Affect.Raw{1} = {};

Data = LoadSelected(Data, hr_path, hr_file, "HR");
Data = LoadSelected(Data, ecg_path, ecg_file, "ECG");
Data = LoadAffect(Data, aff_path, aff_file);

%% RR-Interval Preprocessing

%Bandpass Thresholding
if Bandpass
    [r, c] = size(Data.HR.Raw);
    
    %Create Preprocessed matrix of NaNs
    Data.HR.PP = nan(r,c);
    
    for i = 1:r
        if (Data.HR.Raw(i,3)>= u_band) || (Data.HR.Raw(i,3) <= l_band)
            Data.HR.PP(i,3) = NaN;
        else
            Data.HR.PP(i,3) = Data.HR.Raw(i,3);
        end
    end 
end




%Ectopic Heartbeats
    %Malik Method
if Malik
    %input = Data
    [r,c] = size(Data.HR.Raw);
    
    %Create Preprocessed matrix of NaNs
    Data.HR.PP = nan(r,c);
    
    for i = 1:(r-1)
        if abs(Data.HR.Raw(i,3) - Data.HR.Raw(i+1,3)) > (0.2*Data.HR.Raw(i,3))
            Data.HR.PP(i,3) = NaN;
        else
            Data.HR.PP(i,3) = Data.HR.Raw(i,3);
        end
    end
end
    
    
    %Kamath Method
if Kamath
    %input = Data
    [r,c] = size(Data.HR.Raw);
    
    %Create Preprocessed matrix of NaNs
    Data.HR.PP = nan(r,c);
    
    for i = 1:(r-1)
        a = 0.325 * Data.HR.Raw(i,3);
        b = 0.245 * Data.HR.Raw(i,3);
        
        c = Data.HR.Raw(i+1,3) - Data.HR.Raw(i,3);
        d = Data.HR.Raw(i,3) - Data.HR.Raw(i+1,3);
        
        if (0 <= c) && (c <= a)
            Data.HR.PP(i,3) = NaN;
        elseif (0 <= d) && (d <= b)
            Data.HR.PP(i,3) = NaN;
        else
            Data.HR.PP(i,3) = Data.HR.Raw(i,3);
        end
    end
end
    
    %Karlsson Method
if Karlsson
    %input = Data
    [r,c] = size(Data.HR.Raw);
    
    %Create Preprocessed matrix of NaNs
    Data.HR.PP = nan(r,c);
    
    for i = 2:(r-1)
        a = (Data.HR.Raw(i-1,3)+Data.HR.Raw(i+1,3))/2;
        
        if abs(a-Data.HR.Raw(i,3)) > (0.2*a)
            Data.HR.PP(i,3) = NaN;
        else
            Data.HR.PP(i,3) = Data.HR.Raw(i,3);
        end
        
    end
    
end
    
    
    %Acar Method
if Acar
    %input = Data, acar_range
    [r,c] = size(Data.HR.Raw);
    
    Data.HR.PP = nan(r,c);
    
    for i = (acar_range+1):r
        a = sum(Data.HR.Raw(i-acar_range:i-1,3));
        
        if abs(a-Data.HR.Raw(i,3))> (0.2*a)
            Data.HR.PP(i,3) = NaN;
        else
            Data.HR.PP(i,3) = Data.HR.Raw(i,3);
        end
    end
end
    

% Interpolation
    

%% RR-Interval Analysis
    
    


%% RR-Interval Exports
% This is where files are exported and saved

%Export RR-intervals with affects denoted for Richard
aff = {'SIB','not problem'};

Richard_export(Data, aff, "HR");



%% ECG-Preprocessing



%% ECG Analysis
peak = 800;
dist = 40;
ecg_rr = ecg_rr_conversion(Data, peak, dist);



%% ECG Exports

Richard_export(Data, aff, "ECG");
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
    %videotime = 728
    
    pol_time = (((((11*60)+19)*60)+15)*1000);
    vid_time = 728*1000;
    algn = pol_time - vid_time;
    
    %Check if there is any HR or ECG data loaded. If so, then generate the
    %index numbers for the start and stop.
    
    if iscell(Data.HR.Raw) == 0
        disp("HR data found, generating start and stop indexes");
        Data = time_adjust(Data, algn, "HR");
    end
        
    if iscell(Data.ECG.Raw) == 0
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


%% RR-Interval Preprocessing Functions
function [Data] = LoadSelected(Data, path, file, type)
        
        %Clear previously loaded information (add clear back in for gui)
        %clear Data;
        %Data.HR.Raw{1} = {};
        %Data.ECG.Raw{1} = {};
        
        switch type
            case 'HR'
                %[file, path] = uigetfile('*.txt','MultiSelect','on');
                Data.HR.path = path;
                
                if iscell(file)
                    for i = 1:length(file)
                        dump = data_load(strcat(path,file{i}));
                        Data.HR.Raw{i} = vectorize(dump);
                        clear dump;
                    end
                    Data.HR.files = file;
                    %hr_load_list.Items=file;
                else
                    dump = data_load(strcat(path,file));
                    Data.HR.Raw = vectorize(dump);
                    clear dump;
                    %Display the files that are loaded
                    Data.HR.files = {file};
                    %hr_load_list.Items={file};
                end
                
            case 'ECG'
                %[file,path] = uigetfile('*.txt','MultiSelect','on');
                Data.ECG.path = path;
                
                if iscell(file)
                    for i = 1:length(file)
                        dump = data_load(strcat(path,file{i}));
                        Data.ECG.Raw{i} = vectorize(dump);
                        clear dump;                        
                    end
                    Data.ECG.files = file;
                    %ecg_load_list.Items=file;
                else
                    dump = data_load(strcat(path,file));
                    Data.ECG.Raw = vectorize(dump);
                    clear dump;
                    Data.ECG.files = {file};
                    %ecg_load_list.Items={file};
                end
                
        end
        disp(strcat('Loading file: ', path, file));
    end

function [entries] = entry_select(list, target)
    %Given a list of strings and set of target(s), returns the entry number
    %in the list of each of the targets.
    
    %list : cell array of strings
    %target : cell array of strings or single string
    %entries : matrix vector containing the entry numbers
    
    %create cell array if target is single string
    
    disp(target);
    disp(list);
    if iscell(target) == 0
        target = {target};
    end
    
    if iscell(list) == 0
        list = {list};
    end
    
    matches = ismember(list,target);
    entries = find(matches);

end

function [raw_array] = data_load(fpath)
    %Load in file
    fid = fopen(fpath);
    line = fgetl(fid);

    if length(line) < 80
        disp('HR file detected:');
        format = '%f:%f:%f %f %f %f %f';
    else
        disp('ECG file detected:');
        format = strcat('%f:%f:%f %f', repmat(' %f', 1, 57)); %generate 73 ecg entries
    end

    % Read in information, converting app time into milliseconds
    i = 1;
    while line ~= -1
        nline = textscan(line, format, 'Delimiter', '\t');
        
        if isempty(nline{5}) %Add check for ERROR entry/issue and add skip entry
            disp(strcat("ERROR found in row: ", int2str(i)));
            disp("Skipping data packet");
            line = fgetl(fid);
        else
            ntime = ((((nline{1}*60)+nline{2})*60)+nline{3})*1000; %convert time into milliseconds
            raw_array(i,:) = [ntime, nline(1,4:end)];
            line = fgetl(fid);
            i=i+1;
        end
    end
    disp(strcat(fpath, ': LOADED'));
    raw_array = cell2mat(raw_array);
end

function [new_array] = vectorize(matrix_array)
    %Take the raw array and concatonate all of the data entries into a
    %vector for easy manipulation later
    

    [r, c] = size(matrix_array);
    if c < 10 %might not matter, but keeping in case app format changes
        disp('HR cell array detected:');
        skip = 3;
    else
        disp('ECG cell array detected:');
        skip = 3;
    end

    % Convert array into one were the ecg values/RR intervals are all in
    % one column
    %Lines with no RR intervals are automatically removed during this
    %process, as they contain no valuable information.
    entry = 1;
    new_array=[];
    
    for i = 1:r
        new_array(entry,1:2) = matrix_array(i,1:2);
        for j = skip:c
            if matrix_array(i,j) ~= 0
                new_array(entry,3) = matrix_array(i,j);
                entry=entry+1;
            end
        end
    end
    %Add in NaNs for new rows instead of 0's to make time alignment easier
    
    for i = 1:length(new_array)
        if new_array(i,1)==0
            new_array(i,1) = NaN;
        end
    end
    
end

%% ECG Analysis Functions
function [locs] = ecg_rr_conversion(Data, peak, dist)
% Function to take the ECG data and estimate RR-intervals from them

%inputs:
%   Data: Data structure containing Data.ECG and Data.HR
%   peak: minimum peak prominence for peak
%   dist: minimum distance in index number for peaks


[pks, locs] = findpeaks(Data.ECG.Raw(:,3), 'MinPeakProminence', peak, 'MinPeakDistance', dist);

for i = 2:length(locs)
    %convert to milliseconds assuming a 130Hz sampling frequency
    locs(i,2) = (locs(i,1)-locs(i-1,1))*(100/13);
end
    
%15 on locs = 377 on Data.HR.Raw
disp('done');

end


%% Export Functions
function [] = Richard_export(Data,aff,type)
%This function takes the Data information and saves them as .csv files for
%sending to Richard

%input:
%   Data: The data structure
%   aff: Cell array listing the different Affects that should be denoted as
%   problem behaviors
%   type: What type of data you want to export, "HR" or "ECG"


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
    
    writematrix(new_mat,strcat(type,"_List.csv"));
end