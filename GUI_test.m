%script file for testing gui configurations
close all;
clear all;
PSHR_GUI;

function PSHR_GUI
    %data structure
    Data.HR.Raw{1} = {};
    Data.ECG.Raw{1} = {};
    
    %figure for displaying current HR snapshot
%     hfig = figure('Position', [0 100 1000 600]);
%     hax = axes(hfig);
%     set(hfig, 'Name', 'HR Analysis Output Window');

    
    %figure for displaying options of analysis
    ofig = uifigure('Name','PSHR Analysis Pipeline','Position', [0 100 1000 600]);
    
    %Display Loaded files
    uilabel(ofig, 'Position', [20 570 200 20], 'Text', 'Loaded HR Files');
    hr_load_list = uilistbox(ofig, 'Position', [20 460 300 100],...
        'Items', {'None'}, 'Multiselect','on');%,...
        %'ValueChangedFcn', @(src,event));
    
    uilabel(ofig, 'Position', [420 570 200 20], 'Text', 'Loaded ECG Files');
    ecg_load_list = uilistbox(ofig, 'Position', [340 460 300 100],...
        'Items', {'None'}, 'Multiselect','on');
    
    
    uilabel(ofig, 'Position', [720 570 200 20], 'Text', 'Loaded Affect Files');
    aff_load_list = uilistbox(ofig, 'Position', [660 460 300 100],...
        'Items', {'None'}, 'Multiselect', 'on');
    
    
    
    %% Menu Buttons:
    
    %Load Data
    m = uimenu(ofig,'Text','&Import');
    mHR = uimenu(m,'Text', '&Load HR File');
    mHR.MenuSelectedFcn = @(src,event)LoadSelected('hr');
    mECG = uimenu(m,'Text', '&Load ECG File');
    mECG.MenuSelectedFcn = @(src,event)LoadSelected('ecg');
    
    
    %% Script for analysis list
    pipe.HR={'START'};
    pipe.ECG={'START'};
    uilabel(ofig, 'Position', [20 370 125 15],'Text', 'RR Preprocessing');
    phlist = uilistbox(ofig,'Position', [20 300 125 70], ...
        'Items', {'Bandpass','Ectopic','Other'}, ... 
        'ValueChangedFcn', @(src,event)addpipelist(src,event,'hr'));
        
    phtxt = uitextarea(ofig, 'Position', [220 90 100 282],'Value',pipe.HR);
    
    
    uilabel(ofig, 'Position', [320 270 125 15],'Text', 'ECG Preprocessing');
    pelist = uilistbox(ofig,'Position', [320 200 125 70], ...
        'Items', {'ECG', 'Preprocessing','Placeholders'}, ...
        'ValueChangedFcn',@(src,event)addpipelist(src,event,'ecg'));
    
    petxt = uitextarea(ofig, 'Position', [520 90 100 282],'Value',pipe.ECG);
    
    %Buttons for Preprocessing HR
    
    %Bandpass Thresholding
    %Ectopic Heartbeats
        %Malik Method
        %Kamath Method
        %Karlsson Method
        %Acar Method
    
    %Interpolation Methods:
        %Linear Interpolation

    
    %Tab for each sample, with basic analysis metrics
        %SDSD
        %SDNN
        %RMSSD
        %pNNx
    
    
    %Analysis list:
    disp('done');
    
    function addpipelist(src,event,type)
        switch type
            case 'hr'
                %append to processing list
                pipe.HR{end+1}=src.Value;
                phtxt.Value{end+1}=src.Value;
            case 'ecg'
                pipe.ECG{end+1}=src.Value;
                petxt.Value{end+1}=src.Value;
        end
    end

    function removepipelist(src,event,type)
        switch type
            case 'hr'
                
            case 'ecg'
                
        end
    end
    
    %Load Selected Files by User
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
                    hr_load_list.Items=file;
                else
                    dump = data_load(strcat(path,file));
                    Data.HR.Raw = vectorize(dump);
                    clear dump;
                    %Display the files that are loaded
                    hr_load_list.Items={file};
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
                    ecg_load_list.Items=file;
                else
                    dump = data_load(strcat(path,file));
                    Data.ECG.Raw = vectorize(dump);
                    clear dump;
                    ecg_load_list.Items={file};
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