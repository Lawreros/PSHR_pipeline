% Physmon_Pipeline

% This pipeline relies heavily on the struct class in MATLAB, as the order
% of analyses will not alway be consistent

clear all;

%% Hardcoded Parameters (To be added to GUI)

hr_path = "./sample/HR_03-04-2022.txt";

ecg_path = "./sample/ECG_03-04-2022.txt";


%% Hardcoded Parameters (Immutable)



%% Main Pipeline


%a = select_file_gui("choose from these, dummy", "*");

% Load data
data_hr = data_load(hr_path);
data_ecg = data_load(ecg_path);

% Vectorize
data_hr = vectorize(data_hr);
data_ecg = vectorize(data_ecg);

% Preprocess (basic filtration)
data_hr = bandpass(data_hr, 2000, 900, false);

% Plot the loaded data
plot(data_hr(:,3));
plot(data_ecg(:,3));

% ECG and HR alignment


%% Data Loading and Preprocessing

function [file_list] = select_file_gui(prompt, suffix)
    % str prompt : The text displayed at the top of the window to select
    % files
    % str suffix : The file type to display, put '*' if you want to display
    % all file types
    
    [file,path] = uigetfile(strcat('*.', suffix), prompt,'MultiSelect','on');
    if path == 0
       disp('No file selected, stopping program')
       return;
    end
    
    num = length(file);
    file_list = {};
    for i = 1:num
        file_list{i,1} = strcat(path, file{i});
    end
    
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

function [matrix_array] = bandpass(matrix_array, upper, lower, remove)
    % Looks at the third column of a matrix and either removes the rows which
    % contain that entry or replaces the entry with NaN

    [r,c] = size(matrix_array);
    cut = [];
    for i=1:r
        if (matrix_array(i,3) > upper) || (matrix_array(i,3) < lower)
            cut = [cut,i];
        end
    end
    
    if remove
        matrix_array(cut,:) = [];
    else
        matrix_array(cut,3) = NaN;
    end
end

