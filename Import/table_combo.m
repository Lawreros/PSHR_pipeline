
% Set up example variables for processing
clear all;

aff_1 = {'crying';'SIB';'face_repeat'};
start_1 = {[1,8,10];[8,10];[1,10]};
end_1 = {[4,9,13];[9,13];[7,18]};
%tab_1 = table(start_1, end_1, 'RowNames', aff_1);
tab_1 = [aff_1, start_1, end_1];

aff_2 = {'crying';'SIB';'face_repeat'};
start_2 = {[11,18,20];[18,20];[13,20]};
end_2 = {[14,19,23];[19,23];[14,23]};
%tab_2 = table(start_2, end_2, 'RowNames', aff_2);
tab_2 = [aff_2, start_2, end_2];

table_comb(tab_1, {'crying','SIB'}, 'intersect', tab_2, {'crying'}, 'omit', {'face_repeat'});


function [] = table_comb(varargin)
% Function which takes N matrices consisting of start and end timepoints
% and combines them


% Parameter:
%   ex_end : whether to apply the exclusions at the end or during

% input the different types of tables are parse the results

% Order is:
% table_1, {columns}, {columns}, 'omit', {columns}, {columns}, 'exclude', {columns} 

% Inputing an empty {} into either the omit or exclude will result in all
% other categories that have not already been explicitly specified being
% placed into those categories.


% First, test how best to read in the parameters and organize

% Find the location of all tables in varargin

    disp('starting');
    % Find the location of the table variables in varargin
    j = 0;
    k = 1;
    omit_flag = false;
    exclude_flag = false;
    ex_end = false;
    un_count = 1;
    om_count = 1;
    ex_count = 1;
    
    for i = 1:length(varargin)
        
        if strcmp(varargin{i}, 'omit')
            omit_flag = true;
            exclude_flag = false;
        elseif strcmp(varargin{i},'exclude')
            omit_flag = false;
            exclude_flag = true;
        elseif strcmp(varargin{i}, 'intersect')
            j=0;
            k=k+1;
            omit_flag = false;
            exclude_flag = false;
        elseif strcmp(varargin{i}, 'ex_end')
            ex_end = true;
            omit_flag = false;
            exclude_flag = false;
        end
        
        
        if length(varargin{i}) == 3 && isnumeric(varargin{i}{1,2})
            org_var.table(j+1,k).raw = varargin{i};
            j = j+1;
            un_count=1;
            omit_flag = false;
            exclude_flag = false;
            
            % Create placeholders for omit and exclude
            org_var.table(j,k).omits{1} = NaN;
            org_var.table(j,k).exclude{1} = NaN;
            
            continue
        elseif iscell(varargin{i}) && ~omit_flag
            org_var.table(j,k).unions{un_count} = varargin{i};
            un_count = un_count+1;
        elseif iscell(varargin{i}) && omit_flag
            org_var.table(j,k).omits{om_count} = varargin{i};
            om_count = om_count+1;
        elseif iscell(varargin{i}) && exclude_flag
            org_var.table(j,k).exclude{ex_count} = varargin{i};
            ex_count = ex_count+1;
        end
        
    end
    
    
    % Create the 'master-times' for each table
    for j = 1:size(org_var.table,1)
        for k = 1:length(org_var.table(j,:))
            org_var.table(j,k).new_tab = create_sub_table(org_var.table(j,k), ex_end);
        end
    end
        
    % Combine each table into the final result
    org_var.final_table = create_master_table(org_var.table);
    
    
    master_table = table();

end

function [com_tab] = create_sub_table(tab, ex_end)

    com_tab = table();
    
    
    % Determine whether they want exclusion to happen during or after
    % sub_table creation

    % Create the union table for each union
    % such that it is 
    % {union_1, start, end;
    %  union_2, start, end;
    %  intersect, start, end;
    %  omit, start, end;
    %  exclude (optional), start, end}
    
    % Go through each union pair and also make unified omits
    com_tab = {};
    usd = []; % Keep track of all categories used
    for i = 1:length(tab.unions)
        grb = [];
        for j=1:length(tab.unions{i})
%             find(strcmp(tab.raw(:,1), tab.unions{i}{j}))
            grb = [grb, find(strcmp(tab.raw(:,1), tab.unions{i}{j}))];
        end
        usd = [usd,grb];
        com_tab(i,:) = gen_union(tab.raw(grb,:));
    end
    
    % If more than one union, then they wanted to find the intersection
    if length(com_tab(:,1)) > 1
        com_tab(end+1,:) = gen_inter(com_tab);
    end
    
    
    
    % If omit is empty cell array, then omit all instances that have not been
    % mentioned in unions
    if isempty(tab.omits{:})
        n_grb = [1:1:length(tab.raw(:,1))];
        n_grb(usd) = [];
        com_tab(end+1,:) = gen_union(tab.raw(n_grb,:));
        
        com_tab{end,1} = 'omit';
        
        skip_flg = false;
        
    elseif iscell(tab.omits{:}) % They specify a union of things to omit
        
        for i = 1:length(tab.omits)
            grb = [];
            for j=1:length(tab.omits{i})
                grb = [grb, find(strcmp(tab.raw(:,1), tab.unions{i}{j}))];
            end
            om_tab(i,:) = gen_union(tab.raw(grb,:));
        end
        
        if length(om_tab(:,1)) > 1 % If they want to omit the intersection between categories
            com_tab(end+1,:) = gen_inter(om_tab);
        else
            com_tab(end+1,:) = om_tab;
        end
        
        com_tab{end,1} = 'omit';
        
        skip_flg = false;
        
    else %They have not mentioned omitting at all
        % Do nothing as they don't want to omit anything
        skip_flg = true;
        com_tab(end+1,:) = {'omit', [], []};
    end
    
    
    % If they have chosen to exclude things
    if iscell(tab.exclude{:})
        % Create union of stuff
        
        % If they want to exclude during, then call the function on each of
        % the individual affects
        
        % If they want to exclude after, then call the function on the
        % resulting intersection, post omit I guess...
        
    else
        com_tab(end+1,:) = {'exclude',[],[]};        
    end
    
    
end


function [mas_tab] = create_master_table(tab)
% Create the final table if multiple tables are being combined

    disp('starting');
    
    % If the tables are on the same z-axis, then find union, if on
    % different z-axis find the intersect
    
    [r,c] = size(tab);
    
    for k = 1:c
        clear dump
        for j = 1:r
            % Generate a union across the same z-axis
            
            % Grab from end-2 which will either be the intersect of the single
            % union that was preposed.
            dump(j,:) = tab(j,k).new_tab(end-2,:);
            om_dump(j,:) = tab(j,k).new_tab(end-1,:);
        end
        
        % Apply omits to master_union_tab, as omits made to the unions between tables
        % will inherently carry over to the intersections 
        % ((A U B)*omit(C)) N (D U E) = ((A U B) N (D U E))*omit(C)
        master_union_tab(k,:) = omit(gen_union(dump), gen_union(om_dump));
        
    end
    
    % Generate a list of intersections
    final_tab = gen_inter(master_union_tab);   
    
    % Apply the excludes to the master category if ex_end is selected

    

end



function [un] = gen_union(tab)
% Takes a table of start and stop times and returns the union of said
% combinations

    dump = zeros(1,max([tab{:,3}]));

    for i = 1:length(tab(:,1))
        for j = 1:length(tab{i,2})
            dump(tab{i,2}(j):tab{i,3}(j)) = 1;
        end
    end
    
    % Find the new start and stop times for the union
    vec = find([0, dump]-[dump, 0]);
    starts = vec(1:2:end);
    ends = vec(2:2:end)-1;
    

    % Create the resulting union cell array
    un = {strjoin(tab(:,1),'_U_'), starts, ends};
    
end

function [intersec] = gen_inter(tab)
% Takes a table of start and stop times and returns an intersection between
% all of the rows

    dump = zeros(length(tab(:,1)),max([tab{:,3}]));
    
    for i = 1:length(tab(:,1))
        for j = 1:length(tab{i,2})
            dump(i,tab{i,2}(j):tab{i,3}(j)) = 1;
        end
    end
    
    % Create a binary vector of 1's where all provided categories overlap
    % with eachother
    dump = [sum(dump,1)==length(tab(:,1))];
    
    % Find the new start and stop times for the intersect
    vec = find([0, dump]-[dump, 0]);
    starts = vec(1:2:end);
    ends = vec(2:2:end)-1;
    

    % Create the resulting intersection cell array
    intersec = {strjoin(tab(:,1),'_N_'), starts, ends};

end

function [om_tab] = omit(tab_1, tab_2)
% Omit the entries of tab_2 from tab_1

    dump = zeros(2,max([tab_1{:,3},tab_2{:,3}]));
    
    for j = 1:length(tab_1{1,2})
        dump(1,tab_1{1,2}(j):tab_1{1,3}(j)) = 1;
    end
    
    for j = 1:length(tab_2{1,2})
        dump(2,tab_2{1,2}(j):tab_2{1,3}(j)) = NaN;
    end
    
    
    dump = sum(dump,1)>0;
    
    % Find the new start and stop times for the omission
    vec = find([0, dump]-[dump, 0]);
    starts = vec(1:2:end);
    ends = vec(2:2:end)-1;
    
    % Create the resulting omission cell array
    om_tab = {strcat(tab_1{1,1},'_omit_',tab_2{1,1}), starts, ends};
    
end


function [ex_tab] = exclude(tab_1, tab_2)
% Exclude instance in tab_1 which overlaps with tab_2

    dump = zeros(1,max([tab_1{:,3},tab_2{:,3}]));
    
    % Make binary vector of tab_1's entries
    for j = 1:length(tab_1{1,2})
        dump(1,tab_1{1,2}(j):tab_1{1,3}(j)) = 1;
    end
    
    % Replace all periods from tab_2 with NaNs
    for j = 1:length(tab_2{1,2})
        dump(1,tab_2{1,2}(j):tab_2{i,3}(j)) = NaN;
    end
    
    % Replace all entries of tab_1 with their mean value (so if any entry
    % is NaN, then all will be NaN)
    for j = 1:length(tab_1{1,2})
        dump(1,tab_1{1,2}(j):tab_1{1,3}(j)) = mean(dump(1,tab_1{1,2}(j):tab_1{1,3}(j)));
    end
    
    % Now find all start and stop times (NaN will not be greater than 0)
    dump = dump > 0;
    
    vec = find([9999, dump]-[dump, 9999]);
    starts = vec(1:2:end);
    ends = vec(2:2:end)-1;
    
    % Create the resulting exclude cell array
    ex_tab = {strcat('(',tab_1{1,1},')_exclude_',tab_2{1,1}), starts, ends};
    
end