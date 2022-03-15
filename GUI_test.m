%script file for testing gui configurations
close all;
clear all;
PSHR_GUI;

function PSHR_GUI
    %data structure
    Data.HR.Raw{1} = {};
    Data.ECG.Raw{1} = {};
    
    %figure for displaying current HR snapshot
    hfig = figure('Position', [0 100 1000 600]);
    hax = axes(hfig);
    set(hfig, 'Name', 'HR Analysis Output Window');

    %figure for displaying current ECG snapshot
    efig = figure('Position', [0 500 1000 600]);
    eax = axes(efig);
    set(efig, 'Name', 'ECG Analysis Output Window');
    
    %figure for displaying options of analysis
    ofig = uifigure('Name','PSHR Analysis Pipeline','Position', [0 100 1000 600]);
    
    
    %Buttons to load data:
    m = uimenu(ofig,'Text','&Import');
    
    mHR = uimenu(m,'Text', '&Load HR File');
    mHR.MenuSelectedFcn = @(src,event)LoadSelected('hr', hax);
    mECG = uimenu(m,'Text', '&Load ECG File');
    mECG.MenuSelectedFcn = @(src,event)LoadSelected('ecg', eax);
    
    
    %Script for analysis list
    pipe.HR={'START'};
    pipe.ECG={'START'};
    phlist = uilistbox(ofig,'Position', [50 50 125 70], ...
        'Items', {'Bandpass','Ectopic','Other'}, ... 
        'ValueChangedFcn', @(src,event)updatepipelist(src,event,'hr'));
        
    txt = uitextarea(ofig, 'Position', [125 90 100 82],'Value',pipe.HR);
    
    pelist = uilistbox(ofig,'Position', [100 100 100 100], ...
        'Items', {'ECG', 'Preprocessing','Placeholders'}, ...
        'ValueChangedFcn',@(src,event)updatepipelist(src,event,'ecg'));
    
    
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
    
    function updatepipelist(src,event,type)
        switch type
            case 'hr'
                %append to processing list
                pipe.HR{end+1}=src.Value;
                %update display of what's beign done in the pipeline
                txt.Value{end+1}=src.Value;
            case 'ecg'
                pipe.ECG{end+1}=src.Value;
                txt.Value{end+1}=src.Value;
        end
    end
    
    %Load Selected Files by User
    function LoadSelected(type, axis)
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
                        plot(axis, Data.HR.Raw{i}(:,3));
                        hold on;
                    end
                    hold off;
                else
                    dump = data_load(strcat(path,file));
                    Data.HR.Raw = vectorize(dump);
                    plot(axis, Data.HR.Raw(:,3));
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
                        plot(axis, Data.ECG.Raw{i}(:,3));
                        hold on;
                    end
                    hold off;
                else
                    dump = data_load(strcat(path,file));
                    Data.ECG.Raw = vectorize(dump);
                    plot(axis, Data.ECG.Raw(:,3));
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