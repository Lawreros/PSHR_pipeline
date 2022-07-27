%script file for testing gui configurations
close all;
clear all;
PSHR_GUI;

function PSHR_GUI
    %data structure
    Data.HR.Raw{1} = {};
    Data.ECG.Raw{1} = {};

    
    %figure for displaying options of analysis
    ofig = uifigure('Name','PSHR Analysis Pipeline','Position', [0 100 1000 600]);
    
    %Display Loaded files
    viewselect={};
    uilabel(ofig, 'Position', [20 570 200 20], 'Text', 'Loaded HR Files');
    hr_load_list = uilistbox(ofig, 'Position', [20 480 300 80],...
        'Items', {'None'}, 'Multiselect','on',...
        'ValueChangedFcn', @(src,event)disp(src.Value));
    
    uilabel(ofig, 'Position', [420 570 200 20], 'Text', 'Loaded ECG Files');
    ecg_load_list = uilistbox(ofig, 'Position', [340 480 300 80],...
        'Items', {'None'}, 'Multiselect','on',...
        'ValueChangedFcn', @(src,event)disp(src.Value));
    
    
    uilabel(ofig, 'Position', [720 570 200 20], 'Text', 'Loaded Affect Files');
    aff_load_list = uilistbox(ofig, 'Position', [660 480 300 80],...
        'Items', {'None'}, 'Multiselect', 'on');
    
    
    %View Loaded files Raw Data
    uibutton(ofig,'Position',[20 440 100 22],'Text','View',...
         'ButtonPushedFcn', @(src,event)line_plot(hr_load_list.Value,'HR','Raw'));
    
    uibutton(ofig,'Position',[340 440 100 22],'Text','View',...
         'ButtonPushedFcn', @(src,event)line_plot(ecg_load_list.Value,'ECG','Raw'));
    
    
    %% Menu Buttons:
    
    %Load Data
    m = uimenu(ofig,'Text','&Import');
    mHR = uimenu(m,'Text', '&Load HR File');
    mHR.MenuSelectedFcn = @(src,event)LoadSelected('HR');
    mECG = uimenu(m,'Text', '&Load ECG File');
    mECG.MenuSelectedFcn = @(src,event)LoadSelected('ECG');
    
    
    %% Script for analysis list
    pipe.HR={'START'};
    pipe.ECG={'START'};
    uilabel(ofig, 'Position', [20 410 125 15],'Text', 'RR Preprocessing');
    phlist = uilistbox(ofig,'Position', [20 340 125 70], ...
        'Items', {'Bandpass','Ectopic','Other'}, ... 
        'ValueChangedFcn', @(src,event)addpipelist(src,event,'HR'));
    
    
    uilabel(ofig, 'Position', [20 260 125 15], 'Text', 'RR Analysis');
    ahlist = uilistbox(ofig, 'Position', [20 190 125 70], ...
        'Items', {'Place','Holder', 'Stuff'}, ...
        'ValueChangedFcn', @(src,event)addpipelist(src,event,'HR'));
    
    uilabel(ofig, 'Position', [170 410 125 15],'Text', 'RR Pipeline');
    phtxt = uitextarea(ofig, 'Position', [170 130 140 280],'Value',pipe.HR);
    
    
    
    uilabel(ofig, 'Position', [340 410 125 15],'Text', 'ECG Preprocessing');
    pelist = uilistbox(ofig,'Position', [340 340 125 70], ...
        'Items', {'ECG', 'Preprocessing','Placeholders'}, ...
        'ValueChangedFcn',@(src,event)addpipelist(src,event,'ECG'));
    
    uilabel(ofig, 'Position', [340 260 125 15], 'Text', 'ECG Analysis');
    aelist = uilistbox(ofig,'Position', [340 190 125 70], ...
        'Items', {'Place','Holder','Stuff'}, ...
        'ValueChangedFcn', @(src,event)addpipelist(src,event,'ECG'));
    
    
    uilabel(ofig, 'Position', [490 410 125 15], 'Text', 'ECG Pipeline');
    petxt = uitextarea(ofig, 'Position', [490 130 140 280],'Value',pipe.ECG);
    
    %Buttons for Preprocessing HR
    
    %Bandpass Thresholding
    %Ectopic Heartbeats
        %Malik Method
        %Kamath Method
        %Karlsson Method
        %Acar Method
    
    %Interpolation Methods:
        %Linear Interpolation
    
    
    %Analysis list:
    disp('done');
    
    function addpipelist(src,event,type)
        switch type
            case 'HR'
                %append to processing list
                pipe.HR{end+1}=src.Value;
                phtxt.Value{end+1}=src.Value;
            case 'ECG'
                pipe.ECG{end+1}=src.Value;
                petxt.Value{end+1}=src.Value;
        end
    end
    
    %Load Selected Files by User
    function LoadSelected(type)
        
        %Clear previously loaded information
        clear Data;
        Data.HR.Raw{1} = {};
        Data.ECG.Raw{1} = {};
        
        switch type
            case 'HR'
                [file, path] = uigetfile('*.txt','MultiSelect','on');
                Data.HR.path = path;
                
                if iscell(file)
                    for i = 1:length(file)
                        dump = data_load(strcat(path,file{i}));
                        Data.HR.Raw{i} = vectorize(dump);
                        clear dump;
                    end
                    Data.HR.files = file;
                    hr_load_list.Items=file;
                else
                    dump = data_load(strcat(path,file));
                    Data.HR.Raw = vectorize(dump);
                    clear dump;
                    %Display the files that are loaded
                    Data.HR.files = {file};
                    hr_load_list.Items={file};
                end
                
            case 'ECG'
                [file,path] = uigetfile('*.txt','MultiSelect','on');
                Data.ECG.path = path;
                
                if iscell(file)
                    for i = 1:length(file)
                        dump = data_load(strcat(path,file{i}));
                        Data.ECG.Raw{i} = vectorize(dump);
                        clear dump;                        
                    end
                    Data.ECG.files = file;
                    ecg_load_list.Items=file;
                else
                    dump = data_load(strcat(path,file));
                    Data.ECG.Raw = vectorize(dump);
                    clear dump;
                    Data.ECG.files = {file};
                    ecg_load_list.Items={file};
                end
                
        end
        disp(strcat('Loading file: ', path, file));
    end

    function line_plot(files, type, struct)
        %files : txt files that you want to use
        %type : 'HR' or 'ECG'
        %struct : the name of the structure file you want to use
        
        %Plots what is in Data.<type>.<struct>
        
        [entries] = entry_select(Data.(type).files, files);
        if length(entries) > 1
            
            for i=1:length(entries)
                plot(Data.(type).(struct){i}(:,3));
                hold on;
            end
            hold off;
        else
            %Address annoying fact that if only one file is loaded, then
            %structures can stop being cell arrays
            if iscell(Data.(type).(struct))
                plot(Data.(type).(struct){entries}(:,3));
                hold off;
            else
                plot(Data.(type).(struct)(:,3));
            end
        end
    end

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