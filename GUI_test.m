%script file for testing gui configurations
close all;
clear all;
PSHR_GUI;

function PSHR_GUI
    %data structure
    Data.HR.Raw{1} = {};
    Data.ECG.Raw{1} = {};
    
    %figure for displaying current snapshot
    dfig = figure('Position', [0 100 1000 600]);
    dax = axes(dfig);
    set(dfig, 'Name', 'Analysis Output Window');

    %figure for displaying options of analysis
    ofig = uifigure('Name','PSHR Analysis Pipeline','Position', [0 100 1000 600]);
    
    
    %Buttons to load data:
    m = uimenu(ofig,'Text','&Import');
    
    mHR = uimenu(m,'Text', '&Load HR File');
    mHR.MenuSelectedFcn = @(src,event)LoadSelected('hr');
    mECG = uimenu(m,'Text', '&Load ECG File');
    mECG.MenuSelectedFcn = @(src,event)LoadSelected('ecg');

    
    %Analysis list:
    disp('done');
    
    
    
    
    
    function LoadSelected(type)
        switch type
            case 'hr'
                [file, path] = uigetfile('*.txt','MultiSelect','on');
                Data.HR.path = path;
                Data.HR.files = file;
                
                if iscell(file)
                    for i = 1:length(file)
                        dump = data_load(strcat(path,file{i}));
                        Data.HR.Raw{i} = vectorize(dump);
                        clear dump;
                    end
                else
                    dump = data_load(strcat(path,file));
                    Data.HR.Raw = vectorize(dump);
                    clear dump;
                end
                
            case 'ecg'
                [file,path] = uigetfile('*.txt','MultiSelect','on');
                Data.ECG.path = path;
                Data.ECG.files = file;
                
                if iscell(file)
                    for i = 1:length(file)
                        dump = data_load(strcat(path,file{i}));
                        Data.ECG.Raw{i} = vectorize(dump);
                        clear dump;
                    end
                else
                    dump = data_load(strcat(path,file));
                    Data.ECG.Raw = vectorize(dump);
                    clear dump;
                end
                
        end
        disp(strcat('Loading file: ', path, file));
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