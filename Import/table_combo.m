function [final_table] = table_combo(varargin)
% Function which takes N matrices consisting of start and end timepoints
% and combines them

% input the different types of tables are parse the results

% Order is:
% table_1, {affect(s)}, {affect(s)}, 'omit', {affect(s)}, table_2, {affect(s)},
% {affect(s)}, 'omit',...

% Order of events is 1) (Unions -> Intersections -> omits) individual
% tables, then 2) Unions -> Intersections across tables

% Inputing an empty {} into omit will result in all other categories that
% have not already been explicitly specified being omitted.


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
        
        
        if size(varargin{i},2) == 3 && isnumeric(varargin{i}{1,2})
            org_var.table(j+1,k).raw = varargin{i};
            j = j+1;
            un_count=1;
            om_count=1;
            omit_flag = false;
            exclude_flag = false;
            
            % Create placeholders for unions, omits, and exclude
            org_var.table(j,k).unions{1} = NaN;
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
            if ~isempty(org_var.table(j,k).raw) % Because 2D variable table can have empty entries if more rows in one column than others
                org_var.table(j,k).new_tab = create_sub_table(org_var.table(j,k), ex_end);
            end
        end
    end
        
    % Combine each table into the final result
    final_table = create_master_table(org_var.table);
    
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
    
    
    % If unions is empty cell array, then unify all categories available in
    % the table
    if isempty(tab.unions)
        tab.unions = {tab.raw(:,1).'};
    end
    
    usd = []; % Keep track of all categories used
    com_tab = {};
    if iscell([tab.unions{:}])
        % Go through each union pair and also make unified omits
        for i = 1:size(tab.unions,2)
            grb = [];
            for j=1:length(tab.unions{i})
                grb = [grb, find(ismember(tab.raw(:,1), tab.unions{i}{j}))];
            end
            usd = [usd,grb];
            com_tab(i,:) = gen_union(tab.raw(grb,:));
            
            %Add proper union title, in case some of the categories are not
            %found in tab.raw
            com_tab{i,1} = strjoin(tab.unions{i},'_U_');
        end

        % If more than one union, then they wanted to find the intersection
        if length(com_tab(:,1)) > 1
            com_tab(end+1,:) = gen_inter(com_tab);
        end
    
    else % They have not specifed any categories, thus they just want to use a table
        % for omission or exclusion purposes
        
        com_tab = {'',[],[]};
    end
    
    % If omit is empty cell array, then omit all instances that have not been
    % mentioned in unions
    if isempty(tab.omits{:})
        n_grb = [1:1:length(tab.raw(:,1))];
        n_grb(usd) = [];
        com_tab(end+1,:) = gen_union(tab.raw(n_grb,:));
                
        skip_flg = false;
        
    elseif iscell(tab.omits{1}) % They specify a union of things to omit
        
        for i = 1:size(tab.omits,2)
            grb = [];
            for j=1:length(tab.omits{i})
                grb = [grb, find(ismember(tab.raw(:,1), tab.omits{i}{j}))];
            end
            om_tab(i,:) = gen_union(tab.raw(grb,:));
            
            %Add proper union title, in case some of the categories are not
            %found in tab.raw
            om_tab{i,1} = strjoin(tab.omits{i},'_U_');
            
        end
        
        if length(om_tab(:,1)) > 1 % If they want to omit the intersection between categories
            com_tab(end+1,:) = gen_inter(om_tab);
        else
            com_tab(end+1,:) = om_tab;
        end
        
%         com_tab{end,1} = strcat('omit_', com_tab{end,1});
        
        skip_flg = false;
        
    else %They have not mentioned omitting at all
        % Do nothing as they don't want to omit anything
        skip_flg = true;
        com_tab(end+1,:) = {'',[],[]};%{'omit_Blank', [], []};
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
    
    
    % Add parenthesis around all table operations to better seperate them
    for i = 1:length(com_tab(:,1))
        if ~isempty(com_tab{i,1})
            com_tab{i,1} = strcat('(', com_tab{i,1},')');
        end
    end
    
end


function [final_tab] = create_master_table(tab)
% Create the final table if multiple tables are being combined
    
    % If the tables are on the same z-axis, then find union, if on
    % different z-axis find the intersect
    
    [r,c] = size(tab);
    q = 0;
    
    for k = 1:c
        clear dump
        j = 1;
        while j <= r && ~isempty(tab(j,k).new_tab)
            % Generate a union across the same z-axis
            
            % Grab from end-2 which will either be the intersect of the single
            % union that was preposed.
            dump(j,:) = tab(j,k).new_tab(end-2,:);
            om_dump(q+j,:) = tab(j,k).new_tab(end-1,:);
            j = j+1;
        end
        q=q+j-1;
        
        % Apply omits to master_union_tab, as omits made to the unions between tables
        % will inherently carry over to the intersections 
        % ((A U B)*omit(C)) N (D U E) = ((A U B) N (D U E))*omit(C)
        master_union_tab(k,:) = omit(gen_union(dump), gen_union(om_dump(q-(j-2):q,:)));
        
        % TODO: ADD PARENTHESIS HERE
        
    end
    
    % Generate a list of intersections
    final_tab(1,:) = gen_inter(master_union_tab);
    
    % get the collection of omissions so that with onset analysis it is
    % possible to not select from this area
    final_tab(2,:) = gen_union(om_dump);
                                         
    
    % TODO: Apply the excludes to the master category if ex_end is selected
    
    
end



function [un] = gen_union(tab)
% Takes a table of start and stop times and returns the union of said
% combinations

    %TODO: Can save computation time by subtracting the minimum start time
    %from the creation of the zeros and then adding it to the end vector.
    %(also would work in gen_inter)

    dump = zeros(1,max([tab{:,3}]));
    
    non_empty = [];

    for i = 1:length(tab(:,1))
        
        if ~isempty(tab{i,1})
            non_empty = [non_empty,i];
        end
        if ~isempty(tab{i,2})
            for j = 1:length(tab{i,2})
                dump(tab{i,2}(j):tab{i,3}(j)) = 1;
            end
        end
    end
    
    % Find the new start and stop times for the union
    vec = find([0, dump]-[dump, 0]);
    starts = vec(1:2:end);
    ends = vec(2:2:end)-1;
    

    % Create the resulting union cell array
    if ~isempty(non_empty)
        un = {strjoin(tab(non_empty,1),'_U_'), starts, ends};
    else
        un = {'', starts, ends};
    end
    
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

    % Check that there is something to omit from tab_2
    if isempty([tab_2{:,3}])
        om_tab = tab_1; % Just return tab_1
        if ~isempty([tab_2{:,1}])
%             om_tab{1,1} = strcat(tab_1{1,1}(1:end-1),'_omit_', tab_2{1,1}(2:end));
            om_tab{1,1} = strcat(tab_1{1,1},'_omit_', tab_2{1,1});
        end
    else
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
%         om_tab = {strcat(tab_1{1,1}(1:end-1),'_omit_',tab_2{1,1}(2:end)), starts, ends};
        om_tab = {strcat(tab_1{1,1},'_omit_',tab_2{1,1}), starts, ends};
    end
    
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