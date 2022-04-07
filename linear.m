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
aff_file = "03-18-code.csv";
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

Karlsson = false;

Acar = false;
acar_range = 9;


%ECG-interval



%% Pipeline
Data.HR.Raw{1} = {};
Data.ECG.Raw{1} = {};
Data.Affect.Raw{1} = {};

%Data = LoadSelected(Data, hr_path, hr_file, "HR");
%Data = LoadSelected(Data, ecg_path, ecg_file, "ECG");
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
    
    


%% ECG-Preprocessing



%% ECG Analysis



disp('done');


%% Affect Loading

function [Data] = LoadAffect(Data, path, file)
% Load in affect file and add it to the structure
Data.Affect.path = path;

    if iscell(file)
        for i = 1:length(file)
            Data.Affect.Raw{i} = readtable(strcat(path,file{i}));
        end
        Data.Affect.files = {file};
    else
        Data.Affect.Raw = readtable(strcat(path,file));
        Data.Affect.files = {file};
    end

    
    %Get list of all unique affects used in the coding (SINGLE)
    aff_list = unique(Data.Affect.Raw.Affect1);
    if isnan(unique(Data.Affect.Raw.Affect2))
        disp('No entries in column Affect2');
    else
        aff_list = [aff_list; unique(Data.Affect.Raw.Affect2)]; 
    end
    
    if isnan(unique(Data.Affect.Raw.Affect3))
        disp('No entries in column Affect3');
    else
        aff_list = [aff_list; unique(Data.Affect.Raw.Affect2)]; 
    end
    
    aff_list = unique(aff_list); %cell array of all affects used
    
    
    %Generate Start and End times for the Affects
    Data.Affect.Times = {};
    for i = 1:length(aff_list)
        starts = [];
        ends = [];
        
        buffer = [true, transpose(diff(strcmp(Data.Affect.Raw.Affect1, aff_list{i}))~=0)];
        for j = 1:2:length(buffer)
            starts = [starts, buffer(j)];
            ends = [ends, buffer(j+1)];
        end
        Data.Affect.Times{i,1} = aff_list{i};
        Data.Affect.Times{i,2} = starts;
        Data.Affect.Times{i,3} = ends;
    end
    
    
    
    %Find and read alignment time
    disp("done");


end



function [Aff_raw] = NO_aff_load_raw(path, file)

%Load in file
fid = fopen(strcat(path,file));
line = fgetl(fid);
line = fgetl(fid);
% Read in information
i = 1;
%while line ~= -1
%    Aff_raw(i,:) = textscan(line,'%d %s%s%s %d%d %s', 'delimiter', ',');
%    line = fgetl(fid);
%    i=i+1;
%end
aa = readtable(strcat(path,file));
aa = table2cell(aa);
Aff_raw = aa(:,[1,2,3,4,25,26,27]);


status = strcat(path, file, ': LOADED');
disp(status);
end

function [Affect] = NO_aff_preprocess(affect)

for i=1:length(affect)
    if isempty(affect{i,7})==0
        a = textscan(affect{i,7}, '%s %s %d:%d:%f','delimiter',' ');
        
        Affect.align_time{1,1} = 'Polar time';
        Affect.align_time{1,2} = 'Video time';
        Affect.align_time{2,1} = ((a{3}*60+a{4})*60+a{5})*1000;
        Affect.align_time{2,2} = i*1000;

    end
        
end

% Alignment times
%Affect.align_time{1,1} = 'Polar time';
%Affect.align_time{1,2} = 'Video time';
%Affect.align_time{2,1} = ((((affect{1,1}*60)+affect{1,2})*60)+affect{1,3})*1000;
%Affect.align_time{2,2} = ((((affect{1,4}*60)+affect{1,5})*60)+affect{1,6})*1000;

%if Affect.align_time{2,2} == 0
    corr = Affect.align_time{2,1};
%else
    %corr = 0;
%end


[r, c] = size(affect);
Affect(1).type = 'camera';
Affect(2).type = 'problem_yn';
for q = 1:2
a=0;
new=1;
for i=1:r
    
    if affect{i,4+q} == 1 && a ==0
        a = i;
    end
    
    if (affect{i,4+q} == 0 && a > 0) || (i==r && a>0)
        b = i;
        if new==1
            Affect(q).start = a*1000;
            Affect(q).end = b*1000;
            new=0;
        else
            Affect(q).start = [Affect(q).start,a*1000];
            Affect(q).end = [Affect(q).end,b*1000];
        end
        a=0;
        
    end
end
end
end

function [Seg_HR,X,Y] = NO_affect_analysis(RR, Affect, pre, name)
% Input parameters:
% RR : Heart rate matrix, consisting of time in the first column and
%     whatever you want to isolate in the second column
% Affect : Preprocessed affect matrix
% pre : number of RR values before the start of an affect to include in the
%      segmentation


% Align video and affect time
vid_time=Affect(1).align_time{2,2};
pol_time=Affect(1).align_time{2,1};

T = pol_time-vid_time;
corr_time=T-RR(1,1);


% Align times and replace any zeros with NaN so they don't show up on graph
for i = 2:length(RR)
    RR(i,3) = 0;
    RR(i,4)=0;
    if isnan(RR(i,1))
    else
        RR(i,1) = RR(i,1) - corr_time;
    end
    %    if RR(i,2) == 0
    %        RR(i,2) = NaN;
    %    end
end


Seg_HR = {};

% locate start and stop of each affect and store the values
for i = 1:length(Affect)
    Seg_HR(i).type = Affect(i).type;
    if isempty(Affect(i).start) ~=1
        
        for j = 1:length(Affect(i).start)
            % find where affect instance starts and ends
            a = find(RR(2:end,1) >= Affect(i).start(j));
            b = find(RR(2:end,1) >= Affect(i).end(j));
            
            if isempty(a)==0 && isempty(b)==0
                if a(1) == b(1)
                    %nothing is there
                    %This has the side effect of adding 0's due
                    % to its entry being skipped ...
                    % so [1,2, nothing, 4,5] => [1,2,0,4,5]
                else
                    Seg_HR(i).truestart(j)=a(1);
                    Seg_HR(i).trueend(j)=b(1);
                end
            end
        end
    end
end

%GLM HACK CODE, START:
GLM = RR;
GLM(:,1)=[]; %Get rid of time values, not necessary
GLM(:,2:20)=0;

affect_key = {'SIB';...
    'freezing';...
    'repetitive behaviors';...
    'moving at a fast/abrupt pace';...
    'contorting face/grimacing';...
    'crying';...
    'jumping';...
    'flapping/clapping';...
    'calm sitting';...
    'finger (s) in mouth';...
    'clothing adjustment/removal';...
    'attentive/responsive';...
    'loud/rapid humming';...
    'loud/rapid speech';...
    'soft melodic humming';...
    'laughing';...
    'smiling';...
    'polar strap adjustment/removal';...
    'unresponsive/unable to redirect'};

%GLM HACK CODE, STOP

% Convert RR to cell array for adding strings
RR=num2cell(RR);


% Add "problem behavior" and "on camera" to RR
for i = 1:length(Seg_HR)
    for j=1:length(Seg_HR(i).truestart)
        taken1=0;
        taken2=0;
        
        if Seg_HR(i).truestart(j)~=0 && Seg_HR(i).trueend(j)~=0
                        
            for q = Seg_HR(i).truestart(j):Seg_HR(i).trueend(j)
                if strcmp(Seg_HR(i).type,'camera')
                    RR{q,3} = 1;
                end
                if strcmp(Seg_HR(i).type,'problem_yn')
                    RR{q,4} = 1;
                end
            end
        end
    end
end





% GLM HACK CODE START:
cut=[];
for i=1:length(GLM)
    if isnan(GLM(i,1))
        cut=[cut,i];
    end
end
GLM(cut,:)=[];

X=GLM(:,2:end);
X(:,end+1)=1; %Add column for constant value (normaly low frequency drift correction)

Y=GLM(:,1);

[r,c]=size(X);
cut=[];
for i=1:c
    if sum(X(:,i))==0
        cut=[cut,i];
    end
end
%X(:,cut)=[]; %Remove nuisance regressors that don't have any values (absent affects)

% ONLY COMMENTED OUT FOR BIG ANALYSIS, UNCOMMENT
%disp(cut);
%B = inv(transpose(X)*X)*transpose(X)*Y;
%for i=1:length(B)
%B(i,2)=B(i,1)/(sum(X(:,i)/sum(sum(X)))); %What percentage affect recording this accounts for
%end

% GLM HACK CODE STOP:


% Cut up processed data for affect analysis
for i = 1:length(Seg_HR)
    chunk = {};
    if isempty(Seg_HR(i).truestart) ~=1
        
        for j = 1:length(Seg_HR(i).truestart)
            if Seg_HR(i).truestart(j) <= pre
                chunk{j,1} = RR(1:Seg_HR(i).trueend(j),2:5);
            else
                chunk{j,1} = RR(Seg_HR(i).truestart(j)-pre:Seg_HR(i).trueend(j),2:5);
            end
            
            % Check purity of chunks
            % Check where the target affect type is in the 'chunk'
            [r,c] = size(chunk{j,1});
            
            if r>=1  %just a quick check whether the chunk contains anything TODO: FIX
                if strcmp(Seg_HR(i).type, chunk{j,1}{1,2})
                    target1=3;
                    target2=4;
                elseif strcmp(Seg_HR(i).type, chunk{j,1}{1,3})
                    target1=2;
                    target2=4;
                elseif strcmp(Seg_HR(i).type, chunk{j,1}{1,4})
                    target1=2;
                    target2=3;
                else
                    disp(strcat('ERROR: Affect value corruption in RR ', num2str(Seg_HR(i).truestart(j)-pre),'to', num2str(Seg_HR(i).trueend(j))))
                end
                
                % Check the percentage of "pure" affect and record it
                count=0;
                for q=1:r
                    if isempty(chunk{j,1}{q,target1}) && isempty(chunk{j,1}{1,target2})
                        count=count+1;
                    end
                end
                chunk{j,2}= count/r;
            end
        end
        Seg_HR(i).HR = chunk;
        
    end
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
        format = strcat('%f:%f:%f %f', repmat(' %f', 1, 73)); %generate 73 ecg entries
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
end