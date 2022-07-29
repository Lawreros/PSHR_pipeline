function [Data] = load_affect(Data, path, file)
% Loads the data from the provided affect file and organizes 
    % it into a structure that can be used for future analysis
    %   Inputs:
    %       Data: [cell array or struct] cell array which will be converted into the
    %       structure that is output. If HR or ECG data is already
    %       contained in the Data structure, then affect start and stop
    %       times are calculated using the time_adjust function
    %
    %       path: [string] the path to the folder containing the file you
    %       want to load in
    %
    %       file: [string] name of the file that you want to load
    
    %   Returns:
    %       Data: [struct] structure containing relevant information loaded
    %       from the file you specify
    
    Data.Affect.path = path;

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
        Data.Affect.Raw{1} = readtable(strcat(path,file),'Format','auto');
        Data.Affect.files = {file};
    end

    
    for k = 1:length(Data.Affect.Raw)
    %Get list of all unique affects used in the coding
        aff_list = unique(Data.Affect.Raw{k}.Affect1);
        Data.Affect.Raw{k}.Affect1{1} = "start";
        Data.Affect.Raw{k}.Affect1{end+1} = "end";

        if iscell(unique(Data.Affect.Raw{k}.Affect2(1:end-1)))
            ext = unique(Data.Affect.Raw{k}.Affect2(~cellfun(@isempty, Data.Affect.Raw{k}.Affect2)));
    %         aff_list = [aff_list; unique(Data.Affect.Raw.Affect2)];
            aff_list = [aff_list; ext];
            Data.Affect.Raw{k}.Affect2{1} = "start";
            Data.Affect.Raw{k}.Affect2{end} = "end";
        else
            disp('No entries in column Affect2');
        end

        if iscell(unique(Data.Affect.Raw{k}.Affect3(1:end-1)))
            ext = unique(Data.Affect.Raw{k}.Affect3(~cellfun(@isempty, Data.Affect.Raw{k}.Affect3)));
    %         aff_list = [aff_list; unique(Data.Affect.Raw.Affect3)];
            Data.Affect.Raw{k}.Affect3{1} = "start";
            Data.Affect.Raw{k}.Affect3{end} = "end";
        else
            disp('No entries in column Affect3');
        end

        aff_list = unique(aff_list); %cell array of all affects used


        %Generate Start and End times for the Affects using video time
        Data.Affect.Times{k} = {};
        for i = 1:length(aff_list)
            starts = [];
            ends = [];
            for d = 1:3
                col = strcat("Affect",string(d));

                buffer = [false, transpose(diff(strcmp(Data.Affect.Raw{k}.(col), aff_list{i}))~=0)];
                buffer = find(buffer);


                for j = 1:2:length(buffer)
                    starts = [starts, Data.Affect.Raw{k}.Time_sec(buffer(j))];%buffer(j)];
                    ends = [ends, Data.Affect.Raw{k}.Time_sec(buffer(j+1)-1)];%buffer(j+1)-1];
                end
            end
            Data.Affect.Times{k}{i,1} = aff_list{i};
            Data.Affect.Times{k}{i,2} = starts;
            Data.Affect.Times{k}{i,3} = ends;
        end
    
    end
    
    %Load and store alignment time (Automate this for when format is
    %standardized)
    
    %realtime = "11:19:15"
    %videotime = 728 -1 for lag
    
    pol_time = (((((11*60)+19)*60)+15)*1000);
    vid_time = 727*1000;
    algn = pol_time - vid_time;
    
    %Check if there is any HR or ECG data loaded. If so, then generate the
    %index numbers for the start and stop.
    
    if iscell(Data.HR.Raw) == 1
        disp("HR data found, generating start and stop indexes");
        Data = time_adjust(Data, algn, "HR");
    end
        
    if iscell(Data.ECG.Raw) == 1
        disp("ECG data found, generating start and stop indexes"); 
        Data = time_adjust(Data, algn, "ECG");
    end


end

function [Data] = time_adjust(Data, algn, type)
%Changes the generates indexes for Data.*.Raw that correspond with
%timestamps found in Data.Affect.Times

%inputs:
% Data: main Data structure
% algn: [vector of ints] alignment time found with (polar_timestamp - video_time) = corr
% type: What type of Data you are adjusting the time of

for q=1:length(Data.Affect.Times)

    % Iterate through start and stop times
    [r,c] = size(Data.Affect.Times{q});
    
    % Create place to store new times
    Data.(type).Affect{q} = {};

    
    %Iterate through each affect
    for i = 1:r
        starts = [];
        ends = [];
        
        %Iterate through each start/stop pair
        for j = 1:length(Data.Affect.Times{q}{i,2})
            
            a = find(Data.(type).Raw{q}(:,1) >= (Data.Affect.Times{q}{i,2}(j)*1000+algn));
            %disp(Data.Affect.Times{i,2}(j)*1000+algn);
            
            b = find(Data.(type).Raw{q}(:,1) > (Data.Affect.Times{q}{i,3}(j)*1000+algn));
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
                    Data.Affect.Times{q}{i,1}, Data.Affect.Times{q}{i,2}(j));
                    starts = [starts, a(1)];
                    ends = [ends, length(Data.(type).Raw{q})];
            else
                %Issue due to the affects occuring outside of collected
                %HR/ECG data
                fprintf('WARNING: Affect %s from video time %d to %d could not be found in data\n',...
                    Data.Affect.Times{q}{i,1}, Data.Affect.Times{q}{i,2}(j), Data.Affect.Times{q}{i,3}(j));
                
            end
        end
            % Store indexes for each affect with their respective data
            % types, instead of keeping everything in the Affect structure.
            % This way, every thing you need for each analysis will be
            % grouped
            Data.(type).Affect{q}{i,1} = Data.Affect.Times{q}{i,1};
            Data.(type).Affect{q}{i,2} = starts;
            Data.(type).Affect{q}{i,3} = ends;
    end

end
end
