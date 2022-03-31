% PHYSMON Pipeline
% The purpose of this file is to make debugging the functions in the GUI
% easier, as checking analysis values in a linear script is much easier
% that the constantly updating/changing GUI.

clear all;
close all;

%% Test files:

hr_path = "/home/ross/Documents/MATLAB/PSHR_pipeline/sample/";
hr_file = "HR_03-09-2022.txt";
ecg_path = "/home/ross/Documents/MATLAB/PSHR_pipeline/sample/";
ecg_file = "ECG_03-09-2022.txt";

%% Analysis flags
% As the GUI will most likely chain together preprocessing modules through
% the use of binary triggers, it's good to simulate that with a group of
% True/False statements.


%RR-interval
Bandpass = false;
u_band = 1200;
l_band = 400;

Malik = true;

Karlsson = false;

%ECG-interval



%% Pipeline
Data.HR.Raw{1} = {};
Data.ECG.Raw{1} = {};

Data = LoadSelected(Data, hr_path, hr_file, "HR");
Data = LoadSelected(Data, ecg_path, ecg_file, "ECG");

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
    
    
    
    %Acar Method


%% RR-Interval Analysis
    
    


%% ECG-Preprocessing



%% ECG Analysis



disp('done');


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
        format = strcat('%f:%f:%f %f', repmat(' %f', 1, 73)); %generate 73 ecg entries
    end

    % Read in information, converting app time into milliseconds
    i = 1;
    while line ~= -1
        nline = textscan(line, format, 'Delimiter', '\t');
        ntime = ((((nline{1}*60)+nline{2})*60)+nline{3})*1000; %convert time into milliseconds
        raw_array(i,:) = [ntime, nline(1,4:end)];
        line = fgetl(fid);
        i=i+1;
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
end