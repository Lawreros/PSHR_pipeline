function [Data] = pshr_load(varargin)
% This function serves to load all HR, ECG, and Affect files that are
% required for analysis and organizes them into a structure that can be
% used for future analysis. This function assumes that the indexes for the
% HR/ECG/Affect files is the same for a given recording session, with NaN's
% input for missing files.
%
% EXAMPLE: 
% output = pshr_load('HR', {H_1, H_2, NaN}, 'ECG', {E_1, NaN, E_3}, 'Affect', {A_1, A_2, A_3})

% Optional Parameters:
%   HR: [cell list of strings] List containing the location of the HR files
%       you wish to load and analyze. Default is an empty cell, representing 
%       no files.

%   ECG: [cell list of strings] List containing the location of the ECG
%       files that you wish to load and analyze. Default is an empty cell,
%       representing no files.

%   Affect: [cell list of strings] List containing the location of the
%       coded Affect files that you wish to load and analyze. Default is an
%       empty cell, representing no files.

%   align: [bool] whether to align the affects from the affect files
%       with the present HR/ECG files. Default is true.

%   lag: [int or float] additional amount of time (in seconds) that you
%       wish to lag the video data for when running alignment. For
%       example, inputing 1 will result in the affect at second 274 of the
%       video to be attributed to second 273, attempting to correct for the
%       delay between the polar strap collecting data and transmitting it.
%       Default is 0 and negative values are accepted.

%   verbose: [bool] whether to display additional information/disp
%       statements during running of the function. Default is true.


%   Returns:
%       Data: [struct] structure containing relevant information laoded
%       from the files you specfy
%
    %https://www.mathworks.com/help/matlab/matlab_prog/suppress-warnings.html
    warning('off','MATLAB:table:RowsAddedExistingVars'); %turn off useless warnings
    warning('off','MATLAB:table:ModifiedAndSavedVarnames');

    % Define input parser for the 
    p = inputParser;
    addParameter(p, 'HR', {}, @iscell);
    addParameter(p, 'ECG', {}, @iscell);
    addParameter(p, 'Affect', {}, @iscell);
    addParameter(p, 'align', true, @islogical);
    addParameter(p, 'lag', 0, @isscalar);
    addParameter(p, 'verbose', true, @islogical);

    parse(p,varargin{:});
    
    % Check the length of the different file lists, if their lengths are
    % greater than 1 but not equal to each other, raise an message to the
    % user
    
    
    
    % Create structure that will hold everything (going with Data for easy
    % implementation of old code)
    
    Data = {};
    
    %% Load the files
    
    % Load the HR data and vectorize it
    if isempty(p.Results.HR)==0
        for i = 1:length(p.Results.HR)
            if isnan(p.Results.HR{i}) % Skip NaN files
                disp('Empty HR file');
                Data.HR.Raw{i} = []; % Put empty matrix for consistency
            else
                dump = data_load(p.Results.HR{i});
                Data.HR.Raw{i} = vectorize(dump);
                disp(strcat('Loading file: ', p.Results.HR{i}));
                clear dump; %just save some memory
            end
        end
        Data.HR.files = p.Results.HR;
    end
    
    
    % Load the ECG data and vectorize it
    if isempty(p.Results.ECG)==0
        for i = 1:length(p.Results.ECG)
            if isnan(p.Results.ECG{i}) % Skip NaN files
                disp('Empty ECG file');
                Data.ECG.Raw{i} = []; % Put empty matrix for consistency
            else
                dump = data_load(p.Results.ECG{i});
                Data.ECG.Raw{i} = vectorize(dump);
                disp(strcat('Loading file: ', p.Results.ECG{i}));
                clear dump; %save memory
            end
        end
        Data.ECG.files = p.Results.ECG;
    end
    
    % Load the Affect data
    if isempty(p.Results.Affect)==0
        Data.Affect.files = p.Results.Affect;
        for i = 1:length(p.Results.Affect)
            if isnan(p.Results.Affect{i})
                disp('Empty Affect file');
                Data.Affect.Raw{i} = {};
            else
                Data.Affect.Raw{i} = readtable(p.Results.Affect{i}, 'Format', 'auto');
                
                
                % Get the alignment time
                % Get all unique entries in the Notes column and their
                % locations
                [strings, indx] = unique(Data.Affect.Raw{i}.Notes);
                % Get the location of the cell containing "iPhone reads"
                % for analysis
                loc = indx(contains(strings, 'iPhone reads'));
                
                if isempty(loc)
                    fprintf(strcat("\n\n WARNING: No alignment time found for Affect file: ", Data.Affect.files{i},'\n'));
                    disp("If this is not true, make sure the alignment time is formatted as 'iPhone reads: #:#:#'");
                    fprintf("Currently, a [] will be input for Data.Affect.align_time\n\n");
                    Data.Affect.align_time{i} = [];
                else
                % Find the first location of the string containing the word
                % "iPhone reads" and its index
                % Save both as part of the Data.Affect structure
                    line = Data.Affect.Raw{i}.Notes{loc};
                    format = '%s%s%d:%d:%d';
                    nline = textscan(line, format);
                    
                    if isempty(nline{3}) || isempty(nline{4}) %Hacky way to check for "iPhone reads " with no actual time
                        fprintf(strcat("\n\n WARNING: No alignment time found for Affect file: ", Data.Affect.files{i},'\n'));
                        disp("If this is not true, make sure the alignment time is formatted as 'iPhone reads: #:#:#'");
                        fprintf("Currently, a [] will be input for Data.Affect.align_time\n\n");
                        Data.Affect.align_time{i} = [];
                    else
                        pol_time = ((((nline{3}*60)+nline{4})*60)+nline{5})*1000;
                        vid_time = Data.Affect.Raw{i}.Time_sec(loc);
                        Data.Affect.align_time{i} = [pol_time, vid_time];
                    end
                end
                
                
        
                % Get list of all unique affects used in the coding
                aff_list = unique(Data.Affect.Raw{i}.Affect1);
                Data.Affect.Raw{i}.Affect1{1} = "start";
                Data.Affect.Raw{i}.Affect1{end+1} = "end";
                
                if iscell(unique(Data.Affect.Raw{i}.Affect2(1:end-1)))
                    ext = unique(Data.Affect.Raw{i}.Affect2(~cellfun(@isempty, Data.Affect.Raw{i}.Affect2)));
                    aff_list = [aff_list; ext];
                    Data.Affect.Raw{i}.Affect2{1} = "start";
                    Data.Affect.Raw{i}.Affect2{end} = "end";
                else
                    disp('No entries in column Affect2');
                end

                if iscell(unique(Data.Affect.Raw{i}.Affect3(1:end-1)))
                    ext = unique(Data.Affect.Raw{i}.Affect3(~cellfun(@isempty, Data.Affect.Raw{i}.Affect3)));
                    Data.Affect.Raw{i}.Affect3{1} = "start";
                    Data.Affect.Raw{i}.Affect3{end} = "end";
                else
                    disp('No entries in column Affect3');
                end

                aff_list = unique(aff_list); %cell array of all affects used


                %Generate Start and End times for the Affects using video time
                Data.Affect.Times{i} = {};
                for k = 1:length(aff_list)
                    starts = [];
                    ends = [];
                    for d = 1:3
                        col = strcat("Affect",string(d));

                        buffer = [false, transpose(diff(strcmp(Data.Affect.Raw{i}.(col), aff_list{k}))~=0)];
                        buffer = find(buffer);


                        for j = 1:2:length(buffer)
                            starts = [starts, Data.Affect.Raw{i}.Time_sec(buffer(j))];%buffer(j)];
                            ends = [ends, Data.Affect.Raw{i}.Time_sec(buffer(j+1)-1)];%buffer(j+1)-1];
                        end
                    end
                    Data.Affect.Times{i}{k,1} = aff_list{k};
                    Data.Affect.Times{i}{k,2} = starts;
                    Data.Affect.Times{i}{k,3} = ends;
                end

            end
                
        end
    end
    
    
    % Load and store alignment time (Automate this for when format is
    % standardized)
    
    if p.Results.align
    
        % realtime = "11:19:15"
        % videotime = 728 - 1 for lag

        %pol_time = (((((11*60)+19)*60)+15)*1000);
        %vid_time = 727*1000;
        %algn = pol_time - vid_time;

        %Check if there is any HR or ECG data loaded. If so, then generate the
        %index numbers for the start and stop.

        for q = 1:length(Data.Affect.Times)
            
            % Get the alignment time
            if isempty(Data.Affect.align_time{q})
                disp(strcat('No alignment time found for Affect: ', Data.Affect.files{q}));
                Data.HR.Affect{q} = {};
                Data.ECG.Affect{q} = {};
            else
                algn = Data.Affect.align_time{q}(1) - ((Data.Affect.align_time{q}(2)-p.Results.lag)*1000); %include lag value
            
                if isfield(Data, 'HR') %Added to see if any HR data was loaded/the HR struct field exists
                    if isempty(Data.Affect.Times{q})==0 && isempty(Data.HR.Raw{q})==0
                        fprintf(strcat("\n\nHR data found, generating start and stop indexes for ", Data.HR.files{q},'\n\n'));
                        Data.HR.Affect{q} = time_adjust(Data.HR.Raw{q}, Data.Affect.Times{q}, algn, p.Results.verbose);
                    else
                        Data.HR.Affect{q} = {};
                    end
                end

                if isfield(Data, 'ECG') %Added to see if any ECG data was loaded/the ECG struct field exists
                    if isempty(Data.Affect.Times{q})==0 && isempty(Data.ECG.Raw{q})==0
                        fprintf(strcat("\n\nECG data found, generating start and stop indexes for ", Data.ECG.files{q},'\n\n'));
                        Data.ECG.Affect{q} = time_adjust(Data.ECG.Raw{q}, Data.Affect.Times{q}, algn, p.Results.verbose);
                    else
                        Data.ECG.Affect{q} = {};
                    end
                end
            end
        end

    end
end


function [raw_array] = data_load(fpath)
    %Load in file
    fid = fopen(fpath);
    line = fgetl(fid);
    
    % Count number of entries in line to determine HR or ECG data
    num = length(strfind(line,'	'));

    if num < 5
        disp('HR file detected:');
        format = '%f:%f:%f %f %f %f %f';
    elseif num < 70
        disp('ECG file detected:');
        format = strcat('%f:%f:%f %f', repmat(' %f', 1, 57));
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
    %Add in NaNs for new rows instead of 0's to make time alignment easier
    
    for i = 1:length(new_array)
        if new_array(i,1)==0
            new_array(i,1) = NaN;
        end
    end
    
end



function [results] = time_adjust(mat,times,algn, verbose)
%Changes the generates indexes for Data.*.Raw that correspond with
%timestamps found in Data.Affect.Times

%inputs:
% mat: [n-by-m matrix] either the Raw HR or ECG data found in the Data
% structure (either Data.HR.Raw or Data.ECG.Raw)
% times: [x-by-3 array] the cell array containing the relevant start and stop times for
% the coded data. This is found in the Data structure as Data.Affect.Times
% algn: [vector of ints] alignment time found with (polar_timestamp - video_time) = corr


    % Iterate through start and stop times
    [r,c] = size(times);
    
    % Create place to store new times
    results = {};

    
    %Iterate through each affect
    for i = 1:r
        starts = [];
        ends = [];
        
        %Iterate through each start/stop pair
        for j = 1:length(times{i,2})
            
            a = find(mat(:,1) >= (times{i,2}(j)*1000+algn));
            %disp(Data.Affect.Times{i,2}(j)*1000+algn);
            
            b = find(mat(:,1) > (times{i,3}(j)*1000+algn));
            %disp(Data.Affect.Times{i,3}(j)*1000+algn);
            
            if isempty(a)==0 && isempty(b)==0
                if a(1) == b(1)
                        %nothing is there
                else
                        starts = [starts, a(1)];
                        ends = [ends, b(1)-1];
                end
            elseif isempty(a)==0 && isempty(b)==1
                if verbose
                    fprintf("Affect %s ends after the recording and starts at time %d\n",...
                    times{i,1}, times{i,2}(j));
                end
                    starts = [starts, a(1)];
                    ends = [ends, length(mat)];
            else
                %Issue due to the affects occuring outside of collected
                %HR/ECG data
                if verbose
                    fprintf('WARNING: Affect %s from video time %d to %d could not be found in data\n',...
                    times{i,1}, times{i,2}(j), times{i,3}(j));
                end
            end
        end
            % Store indexes for each affect with their respective data
            % types, instead of keeping everything in the Affect structure.
            % This way, every thing you need for each analysis will be
            % grouped
            results{i,1} = times{i,1};
            results{i,2} = starts;
            results{i,3} = ends;
    end

end

