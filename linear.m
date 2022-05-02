% PHYSMON Pipeline
% The purpose of this file is to make debugging the functions in the GUI
% easier, as checking analysis values in a linear script is much easier
% that the constantly updating/changing GUI.

clear all;
close all;

%% Test files:

%HR file
hr_path = '/home/ross/Downloads/';
hr_file = 'HR_03-18-2022.txt';
%ECG file
ecg_path = '/home/ross/Downloads/';
ecg_file = 'ECG_03-18-2022.txt';
%Affect file
aff_path = "/home/ross/Documents/MATLAB/PSHR_pipeline/sample/";
aff_file = "2022-03-18_Kessler.csv";


%Export RR-intervals with affects denoted for Richard
aff = {'SIB','not problem','Freezing','loud/rapid speech'};

%
%realtime = "11:19:15"
%videotime = 728

%% Analysis flags
% As the GUI will most likely chain together preprocessing modules through
% the use of binary triggers, it's good to simulate that with a group of
% True/False statements.


%RR-interval Preprocess Flags
Bandpass = true;
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

%ECG-interval Preprocess Flags


%ECG-interval Analysis Flags



%% Pipeline
Data.HR.Raw{1} = {};
Data.ECG.Raw{1} = {};
Data.Affect.Raw{1} = {};

Data = LoadSelected(Data, hr_path, hr_file, "HR");
%Data = LoadSelected(Data, ecg_path, ecg_file, "ECG");
Data = LoadAffect(Data, aff_path, aff_file);

%% RR-Interval Preprocessing

%Bandpass Thresholding
if Bandpass
    Data = bandpass(Data,"Raw", l_band, u_band, false);
end

%Ectopic Heartbeats
    %Malik Method
if Malik
    Data = malik(Data, "Raw", false);
end
    
    
    %Kamath Method
if Kamath
    Data = malik(Data, "Raw", false);
end
    
    %Karlsson Method
if Karlsson
    Data = karlsson(Data, "Raw", false);
end
    
    
    %Acar Method
if Acar
    Data = acar(Data, "Raw", acar_range, false);
end
    

% Interpolation
    


%% RR-Interval Analysis

qqq = pnnx_calc_2(Data.HR.Raw(:,3),50,{5,'second'},false);
ppp = rmssd_calc_2(Data.HR.Raw(:,3),{5,'seconds'},false);

% pNN50
if pnnx
    Data = pnnx_calc(Data, "Raw", 50, false);
end

% SDNN
if sdnn
    Data = sdnn_calc(Data, "Raw", false);
end

% SDSD
if sdsd
    Data = sdsd_calc(Data, "Raw", false);
end

% RMSSD
if rmssd
    Data = rmssd_calc(Data, "Raw", false);
end

%% Simple Plot of raw RR data
fig1 = figure(1);
plot(Data.HR.Raw(:,3))
title("Raw RR-interval Data");

%% RR-Interval Exports
% This is where files are exported and saved

Richard_export(Data, aff, "HR", "test_HR_output");

%% ECG-Preprocessing



%% ECG Analysis
peak = 800;
dist = 40;
ecg_rr = ecg_rr_conversion(Data, peak, dist); %This function converts the ECG data into RR-intervals through peak detection


%% Simple Plot of raw ECG data
fig2 = figure(2);
plot(Data.ECG.Raw(:,3));
title("Raw ECG Data");

%% ECG Exports

Richard_export(Data, aff, "ECG", "test_ECG_output");
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

function [Data] = bandpass(Data, source, l_band, u_band, rang)
    % Applies a bandpass filtering to the HR data provided. Any RR-interval
    % value outside of the range specified by the lower and upper bounds is
    % replaced with a NaN
    %   Inputs:
    %       Data: The Data structure
    %       source: [string], Which matrix from the structure you want to
    %       use
    %       l_band: [int], the lower bounding value for the bandpass filter
    %       u_band: [int], the upper bounding value for the bandpass filter
    %       rang: [2 int vector] The range [start, end] of values you want
    %       to use the bandpass on. If false, then analyze the whole
    %       range

    if rang
        r_1 = rang(1);
        r_2 = rang(2);
    else
        [r_2,c] = size(Data.HR.(source));
        r_1 = 1;
    end

    %Create copy of matrix to edit
    Data.HR.PP = Data.HR.(source);
    
    for i = r_1:r_2
        if (Data.HR.(source)(i,3)>= u_band) || (Data.HR.(source)(i,3) <= l_band)
            Data.HR.PP(i,3) = NaN;
        else
            Data.HR.PP(i,3) = Data.HR.(source)(i,3);
        end
    end 
end

function [Data] = malik(Data, source, band)
    % Applies the malik filtering method to the data provided. Any
    % RR-inverval which is outside of the acceptable bounds will be
    % replaced with a NaN
    %   Inputs:
    %       Data: The Data structure
    %       source: [string], Which matrix from the structure you want to
    %       use
    %       band: [2 int vector] The range [start, end] of values you want
    %       to use the malik filter on. If false, then analyze the whole
    %       range
    
    if band
        r_1 = band(1);
        r_2 = band(2);
    else
        [r_2, c] = size(Data.HR.(source));
        r_1 = 1;
    end
        
    %Create copy of matrix to edit
    Data.HR.PP = Data.HR.(source);
    
    for i = r_1:(r_2-1)
        if abs(Data.HR.(source)(i,3) - Data.HR.(source)(i+1,3)) > (0.2*Data.HR.(source)(i,3))
            Data.HR.PP(i,3) = NaN;
        else
            Data.HR.PP(i,3) = Data.HR.(source)(i,3);
        end
    end

end

function [Data] = kamath(Data, source, band)
    % Applies the kamath filtering method to the data provided. Any
    % RR-interval which is outside of the acceptable bounds will be
    % replaced with a NaN
    %   Inputs:
    %       Data: The Data structure
    %       source: [string], Which matrix from the structure you want to
    %       use
    %       band: [2 int vector] The range [start, end] of values you want
    %       to use the kamath filter on. If false, then analyze the whole
    %       range

    
    if band
        r_1 = band(1);
        r_2 = band(2);
    else
        [r_2, c] = size(Data.HR.(source));
        r_1 = 1;
    end
    
    %Create copy of matrix to edit
    Data.HR.PP = Data.HR.(source);
    
    for i = r_1:(r_2-1)
        a = 0.325 * Data.HR.(source)(i,3);
        b = 0.245 * Data.HR.(source)(i,3);
        
        c = Data.HR.(source)(i+1,3) - Data.HR.(source)(i,3);
        d = Data.HR.(source)(i,3) - Data.HR.(source)(i+1,3);
        
        if (0 <= c) && (c <= a)
            Data.HR.PP(i,3) = NaN;
        elseif (0 <= d) && (d <= b)
            Data.HR.PP(i,3) = NaN;
        else
            Data.HR.PP(i,3) = Data.HR.(source)(i,3);
        end
    end

end

function [Data] = karlsson(Data, source, band)
    % Applies the Karlsson filtering method to the data provided. Any
    % RR-interval which is outside of the acceptable bounds will be
    % replaced with a NaN
    %   Inputs:
    %       Data: The Data structure
    %       source: [string], Which matrix from the structure you want to
    %       use
    %       band: [2 int vector] The range [start, end] of values you want
    %       to use the Karlsson filter on. If false, then analyse the whole
    %       range

    if band
        r_1 = band(1);
        r_2 = band(2);
    else
        [r_2, c] = size(Data.HR.(source));
        r_1 = 2;
    end
        
    %Create copy of matrix to edit
    Data.HR.PP = Data.HR.(source);
    
    for i = r_1:(r_2-1)
        a = (Data.HR.(source)(i-1,3)+Data.HR.(source)(i+1,3))/2;
        
        if abs(a-Data.HR.(source)(i,3)) > (0.2*a)
            Data.HR.PP(i,3) = NaN;
        else
            Data.HR.PP(i,3) = Data.HR.(source)(i,3);
        end
        
    end
    
end
        
function [Data] = acar(Data, source, acar_range, band)
    % Applies the Acar filtering method to the data provided. Any
    % RR-interval which is outside of the acceptable bounds will be
    % replaced with a NaN
    %   Inputs:
    %       Data: The Data structure
    %       source: [string], Which matrix from the structure you want to
    %       use
    %       acar_range: [int] 
    %       band: [2 int vector] The range [start, end] of values you wish
    %       to use the Acar filter on. If false, then analyze the full
    %       range


    if band
        r_1 = band(1);
        r_2 = band(2);
    else
        [r_2, c] = size(Data.HR.(source));
        r_1 = acar_range+1;
    end
    
    if r_1 < acar_range
        disp("starting point less than acar range, changing starting point to (acar range)+1");
        r_1 = acar_range+1;
    end
    
    Data.HR.PP = Data.HR.(source);
    
    for i = r_1:r_2
        a = sum(Data.HR.(source)(i-acar_range:i-1,3));
        
        if abs(a-Data.HR.(source)(i,3))> (0.2*a)
            Data.HR.PP(i,3) = NaN;
        else
            Data.HR.PP(i,3) = Data.HR.(source)(i,3);
        end
    end
end


%% RR-Interval Analysis Functions

function [Data] = pnnx_calc(Data,source,diff,band)
    % Calculates the percentage of adjacent NN-intervals that differ from
    % each other by more than "diff" milliseconds
    %   Inputs:
    %       Data: The Data structure
    %       source: [string], Which matrix from the structure you want to
    %       use
    %       diff: [int] The minimum difference in milliseconds between
    %       successive NN-intervalse that you want to count
    %       band: [2 int vector] The range [start, end] of values you want
    %       to calculate the pnnX of. If false, then analyze the whole
    %       range
    
    if band
        r_1 = band(1);
        r_2 = band(2);
    else
        [r_2,c] = size(Data.HR.(source));
        r_1 = 2;
    end
    
    count = 0;
    for i = r_1:r_2
        if abs(Data.HR.(source)(i,3)-Data.HR.(source)(i-1,3))>= diff
            count = count+1;
        end
    end

    Data.HR.(strcat('pNN', string(diff))).(strcat('time_',string(r_1),'_to_',string(r_2))) = count/(r_2-r_1);
    disp(strcat('pNN', string(diff), "=", string(Data.HR.(strcat('pNN', string(diff))).(strcat('time_',string(r_1),'_to_',string(r_2))))));
    
end

function [Data] = sdnn_calc(Data,source,band)
    % Calculates the standard deviation of "NN-intervals" (or preprocessed
    % RR-intervals)
    %   Inputs:
    %       Data: The Data structure
    %       source: [string], Which matrix from the structure you want to
    %       use
    %       band: [2 int vector] The range [start, end] of values you want
    %       to calculade the SDNN of. If false, then analyze the whole
    %       range

    if band
        r_1 = band(1);
        r_2 = band(2);
    else
        [r_2,c] = size(Data.HR.(source));
        r_1 = 1;
    end
    
    Data.HR.SDNN = std(Data.HR.(source)(r_1:r_2,3));
    disp(strcat('SDNN calculated:', string(Data.HR.SDNN)));

end

function [Data] = sdsd_calc(Data,source,band)
    % Calculates the standard deviation of successive differences (SDSD)
    % for the vector of measurements provided
    %   Inputs:
    %       Data: The Data structure
    %       source: [string], Which matrix from the structure you want to use (this is
    %       here to allow you to analyze both Raw and preprocessed)
    %       band: [2 int vector] The range [start, end] of values you want to calculate SDSD for.
    %       If false, then analyze the whole range
    
    if band
        r_1 = band(1);
        r_2 = band(2);
    else
        [r_2,c] = size(Data.HR.(source));
        r_1 = 2;
    end
    
    Data.HR.SDSD = std(Data.HR.(source)(r_1:r_2-1,3)-Data.HR.(source)(r_1+1:r_2,3));
    disp(strcat('SDSD calculated:', string(Data.HR.SDSD)));

end

function [Data] = rmssd_calc(Data, source, band)
    % Calculates the root mean square of successive differences (RMSSD) for
    % the vector of mesasurements provided.
    % Inputs:
    %       Data: The Data structure
    %       source: [string], Which matrix from the structure you want to
    %       use
    %       band: [2 int vector] The range [start, end] of values you want
    %       to calculate RMSSD for. If false, then analyze the whole range

    if band
        r_1 = band(1);
        r_2 = band(2);
    else
        [r_2,c] = size(Data.HR.(source));
        r_1 = 1;
    end
    
    summation = 0;
    for  i = r_1:(r_2-1)
        summation = summation+(Data.HR.(source)(i+1,3)-Data.HR.(source)(i,3))^2;
    end
    
    Data.HR.RMSSD = sqrt(1/((r_2-r_1)-1)*summation);
    disp(strcat('RMSSD calculated:', string(Data.HR.RMSSD)));
    
end


% Change the above analysis functions to not depend on Data structure
function [ret] = pnnx_calc_2(mat,diff,bin,band)
    % Calculates the percentage of adjacent NN-intervals that differ from
    % each other by more than "diff" milliseconds
    %   Inputs:
    %       mat: A [n-by-1] vector which contains the data you want to
    %       calculate pNNX for
    %       diff: [int] The minimum difference in milliseconds between
    %       successive NN-intervalse that you want to count
    %
    %       bin: [1-by-2 cell array] Used for creating a vector of the pNNX
    %       results from a sliding bin of Y seconds or entries. This takes the
    %       format of {index, 'units'}, so if you want to have a bin of the
    %       last 5 seconds: {5, 'second'} or if you want the last 5 measurements: {5, 'measure'}
    %       If you don't want this, set bin to false.
    %
    %       band: [2 int vector] The range [start, end] of values you want
    %       to calculate the pnnX of. If false, then analyze the whole
    %       range
    
    if band
        r_1 = band(1);
        r_2 = band(2);
    else
        [r_2,c] = size(mat);
        r_1 = 2;
    end
    
    % If they've decided to use the bin values
    if iscell(bin)
        a = bin{1}; % value
        b = bin{2}; % units
        
        ret = zeros(r_2-r_1,1);
        
        if strcmp(b,'second')
            for i = r_1:r_2
                count = 0;
                j = 0;
                
                % Check loop backward until you have reached the 'b'
                % seconds in the past through summing
                while (sum(mat(i-j:i,1))/1000) <= a
                    j = j+1;
                    if j == i
                        break;
                    end
                end
                
                if j > 1 && j < i % If there is more than one entry
                    for k = (i-j+2):i
                        if abs(mat(k,1) - mat(k-1,1))>= diff
                            count = count+1;
                        end
                    end
                    ret(i-r_1+1,1) = count/(j-1);
                else
                    ret(i-r_1+1,1) = 0;
                end
        
            end
        else % Looking at the past 'b' entries for the calculation
            for i = (a+1):(r_2-r_1)
                count = 0;
                for j = 1:a
                    if abs(mat(i-j,1) - mat(i-j-1,1))>= diff
                        count = count+1;
                    end
                end
                ret(i,1) = count;
            end
        end
        
    else
        % If they just want a percentage for a matrix
        count = 0;
        for i = r_1:r_2
            if (mat(i,3)-mat(i-1,3))>= diff
                count = count+1;
            end
        end
        ret = count/(r_2-r_1-1);
    end
end

function [ret] = rmssd_calc_2(mat,bin,band)
    % Calculates the root mean square of successive differences (RMSSD) for
    % the vector of mesasurements provided.
    %   Inputs:
    %       mat: A [n-by-1] vector which contains the data you want to
    %       calculate RMSSD for
    %
    %       bin: [1-by-2 cell array] Used for creating a vector of the
    %       RMSSD
    %       results from a sliding bin of Y seconds or entries. This takes the
    %       format of {index, 'units'}, so if you want to have a bin of the
    %       last 5 seconds: {5, 'second'} or if you want the last 5 measurements: {5, 'measure'}
    %       If you don't want this, set bin to false.
    %
    %       band: [2 int vector] The range [start, end] of values you want
    %       to calculate the RMSSD of. If false, then analyze the whole
    %       range
    
    if band
        r_1 = band(1);
        r_2 = band(2);
    else
        [r_2,c] = size(mat);
        r_1 = 2;
    end
    
    %If they've decided to use the bin values
    if iscell(bin)
        a = bin{1}; % value
        b = bin{2}; % units
        
        ret = zeros(r_2-r_1,1);
        
        if strcmp(b,'second')
            for i = r_1:r_2
                summation = 0;
                j = 0;
                
                % Check loop backward until you have reached the 'b'
                % seconds in the past through summing
                while (sum(mat(i-j:i,1))/1000) <= a
                    j = j+1;
                    if j == i
                        break;
                    end
                end
                
                if j > 1 && j < i % If there is more than one entry
%                     for k = (i-j+2):i
%                         if abs(mat(k,1) - mat(k-1,1))>= diff
%                             count = count+1;
%                         end
%                     end
%                    ret(i-r_1+1,1) = count/(j-1);
                    for  k = (i-j+2):i
                        summation = summation+(mat(k,1)-mat(k-1,1))^2;
                    end
                    ret(i-r_1+1,1) = sqrt((1/(j-2))*summation);
                else
                    ret(i-r_1+1,1) = 0;
                end
            end
        else % Looking at the past 'a' entries for the calculation
            for i = (a+1):(r_2-r_1)
%                 count = 0;
%                 for j = 1:a
%                     if abs(mat(i-j,1) - mat(i-j-1,1))>= diff
%                         count = count+1;
%                     end
%                 end
%                 ret(i,1) = count;
                summation = 0;
                for j = 1:a
                    summation = summation+(mat(i-j,1)-mat(i-j-1,1))^2;
                end
                ret(i,1) = sqrt((1/(a-1))*summation);
            end
        end
        
    else
        % If they just want a single value for the input vector
        summation = 0;
        for  i = r_1:(r_2-1)
            summation = summation+(mat(i+1,3)-mat(i,3))^2;
        end
        ret = sqrt(1/((r_2-r_1)-1)*summation);
    end

end

function [ret] = sdnn_calc_2(mat,bin,band)
    % Calculates the standard deviation for the vector of mesasurements provided.
    %   Inputs:
    %       mat: A [n-by-1] vector which contains the data you want to
    %       calculate standard deviation for
    %       bin: [1-by-2 cell array] Used for creating a vector of the
    %       standard deviations
    %       results from a sliding bin of Y seconds or entries. This takes the
    %       format of {index, 'units'}, so if you want to have a bin of the
    %       last 5 seconds: {5, 'second'} or if you want the last 5 measurements: {5, 'measure'}
    %       If you don't want this, set bin to false.
    %
    %       band: [2 int vector] The range [start, end] of values you want
    %       to calculate the RMSSD of. If false, then analyze the whole
    %       range
    
    if band
        r_1 = band(1);
        r_2 = band(2);
    else
        [r_2,c] = size(mat);
        r_1 = 1;
    end
    
    %If they've decided to use the bin values
    if iscell(bin)
        a = bin{1}; % value
        b = bin{2}; % units
        
        ret = zeros(r_2-r_1,1);
        
        if strcmp(b,'second')
            for i = r_1:r_2
                j = 0;
                
                % Check loop backward until you have reached the 'b'
                % seconds in the past through summing
                while (sum(mat(i-j:i,1))/1000) <= a
                    j = j+1;
                    if j == i
                        break;
                    end
                end
                
                if j > 1 && j < i % If there is more than one entry
                    ret(i-r_1+1,1) = std(mat((i-j+1:i),1));
                else
                    ret(i-r_1+1,1) = 0;
                end
            end
        else % Looking at the past 'a' entries for the calculation
            for i = (a+1):(r_2-r_1)
                ret(i,1) = std(mat(i-a:i,1));
            end
        end
        
    else
        % If they just want a single value for the input vector
        for  i = r_1:(r_2-1)
            summation = summation+(mat(i+1,3)-mat(i,3))^2;
        end
        ret = sqrt(1/((r_2-r_1)-1)*summation);
    end

end

function [ret] = sdsd_calc_2(mat,bin,band)
    % Calculates the standard deviation of successive differences (SDSD)
    % for the vector of measurements provided
    %   Inputs:
    %       mat: A [n-by-1] vector which contains the data you want to
    %       calculate SDSD for
    %       bin: [1-by-2 cell array] Used for creating a vector of the
    %       standard deviations
    %       results from a sliding bin of Y seconds or entries. This takes the
    %       format of {index, 'units'}, so if you want to have a bin of the
    %       last 5 seconds: {5, 'second'} or if you want the last 5 measurements: {5, 'measure'}
    %       If you don't want this, set bin to false.
    %       band: [2 int vector] The range [start, end] of values you want
    %       to calculate the RMSSD of. If false, then analyze the whole
    %       range
    
    if band
        r_1 = band(1);
        r_2 = band(2);
    else
        [r_2,c]=size(mat);
        r_1=2;
    end
    
    %If they've decided to use the bin values
    if iscell(bin)
        a = bin{1}; % value
        b = bin{2}; % units
        
        ret = zeros(r_2-r_1,1);
        
        if strcmp(b, 'second')
            for i = r_1:r_2
                j = 0;
                
                while (sum(mat(i-j:i,1))/1000) <= a
                    j = j+1;
                    if j == i
                        break;
                    end
                end
                
                if j > 1 && j < i % If there is more than one entry
                    ret(i-r_1+1,1) = std(mat((i-j:i),1) - mat((i-j+1:i)));
                else
                    ret(i-r_1+1,1) = 0;
                end
            end
        else
            for i = (a+1):(r_2-r_1)
                ret(i,1) = std(mat((i-a:i),1) - mat((i-a+1:i)));
            end
        end
    else
        Data.HR.SDSD = std(mat(r_1:r_2-1,1)-mat(r_1+1:r_2,1));
    end
end

%% RR-Interval Plotting
function [] = poincare_plot(fname, RR, plot_num, max_plot_num)
    % Generates a poincare plot from the data
    % Inputs:
    %       RR: [n-by-1 array], vector containing all RR values
    %       plot_num: [int], plot position on the subplot figure
    %       max_plot_number: [int], maximum intended entries into subplot figure


    cut = [];
    for i=1:length(RR)              % Eliminate NaN's
        if RR(i,1) == 0 || isnan(RR(i,1))
            cut = [cut,i];
        end
    end
    RR(cut) = [];
    RR(:,2) = [RR(2:end,1);0];
    RR(end,:) = [];

    % Calculate SD1 and SD2
    xc = sum(RR(1:end-1,1))/(length(RR)-1);
    yc = sum(RR(2:end,2))/(length(RR)-1);

    SD1 = sqrt((1/length(RR))*nansum(((RR(:,1)-RR(:,2))-nanmean(RR(:,1)-RR(:,2))).^2)/2);
    SD2 = sqrt((1/length(RR))*nansum(((RR(:,1)+RR(:,2))-nanmean(RR(:,1)+RR(:,2))).^2)/2);

    % Making a rotated elipsoid to display SD1 and SD2 https://www.mathworks.com/matlabcentral/answers/342148-how-can-i-rotate-an-ellipse
    x = zeros(1000,1);
    y = zeros(1000,1);
    theta = linspace(0, 2*pi, 1000);
    for k = 1:1000
        x(k) = SD2*cos(theta(k));
        y(k) = SD1*sin(theta(k));
    end
    alpha = pi/4;
    R = [cos(alpha) -sin(alpha); sin(alpha) cos(alpha)];
    rCoords = R*[x'; y'];
    xr = rCoords(1,:)';
    yr = rCoords(2,:)';


    if max_plot_num <= 10    %calculate number of rows and columns for subplots
        if max_plot_num == 1
            p_c = 1;
        else
            p_c = 2;
        end
        p_r = ceil(max_plot_num/p_c);
    else
        p_c = floor(sqrt(max_plot_num));
        p_r = ceil(max_plot_num/p_c);
    end

    min_RR = nanmin(RR(:,1));
    max_RR = nanmax(RR(:,1));


    figure(1);
    % plot onto subplot
    subplot(p_r,p_c,plot_num)
    scatter(RR(:,1), RR(:,2), 15)
    axis([min_RR-50 max_RR+50 min_RR-50 max_RR+50])
    xlabel('RR_n (ms)');
    ylabel('RR_n_+_1 Interval (ms)');
    title(fname)
    hold on;
    plot(xr+xc, yr+yc, 'r');
    plot([0:1600],[0:1600],'r');
    hold off;

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


end