function [Data] = pshr_load_data(Data, path, file, type)
    % Loads the data from the provided HR or ECG file and organizes 
    % it into a structure that can be used for future analysis
    %   Inputs:
    %       Data: [cell array] cell array which will be converted into the
    %       structure that is output
    %
    %       path: [string] the path to the folder containing the file you
    %       want to load in
    %
    %       file: [string] name of the file that you want to load
    %
    %       type: [string] the type of the data in the file you are
    %       loading, either "HR" or "ECG"
    
    %   Returns:
    %       Data: [struct] structure containing relevant information loaded
    %       from the file you specify
    

        
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