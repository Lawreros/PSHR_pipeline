
% Set up example variables for processing
disp('test');

aff_1 = ["crying";"SIB";"face_repeat"];
start_1 = {[1,8,10];[8,10];[1,10]};
end_1 = {[4,9,13];[9,13];[4,13]};
tab_1 = table(aff_1, start_1, end_1);


aff_2 = ["crying";"SIB";"face_repeat"];
start_2 = {[1,8,10];[8,10];[1,10]};
end_2 = {[14,19,23];[19,23];[14,23]};
tab_2 = table(aff_2, start_2, end_2);


table_comb(tab_1, {'crying','face_repeat'}, {'SIB'},"omit",{}, tab_2, 'ex_end');


function [] = table_comb(varargin)
% Function which takes N matrices consisting of start and end timepoints
% and combines them


% Parameter:
%   ex_end : whether to apply the exclusions at the end or during

% input the different types of tables are parse the results

% Order is:
% table_1, {columns}, {columns}, "omit", {columns}, {columns}, "invalid", {columns} 

% Inputing an empty {} into either the omit or exclude will result in all
% other categories that have not already been explicitly specified being
% placed into those categories.


% First, test how best to read in the parameters and organize

% Find the location of all tables in varargin

    disp('starting');
    % Find the location of the table variables in varargin
    j = 0;
    omit_flag = false;
    exclude_flag = false;
    ex_end = false;
    un_count = 1;
    om_count = 1;
    ex_count = 1;
    
    for i = 1:length(varargin)
        
        if istable(varargin{i})
            org_var.table(j+1).raw = varargin{i};
            j = j+1;
            un_count=1;
            omit_flag = false;
            exclude_flag = false;
        elseif strcmp(varargin{i}, 'omit')
            omit_flag = true;
            exclude_flag = false;
        elseif strcmp(varargin{i},'exclude')
            omit_flag = false;
            exclude_flag = true;
        elseif strcmp(varargin{i}, 'ex_end')
            ex_end = true;
            omit_flag = false;
            exclude_flag = false;
        end
        
        if iscell(varargin{i}) && ~omit_flag
            org_var.table(j).unions{un_count} = varargin{i};
            un_count = un_count+1;
        elseif iscell(varargin{i}) && omit_flag
            org_var.table(j).omits{om_count} = varargin{i};
            om_count = om_count+1;
        elseif iscell(varargin{i}) && exclude_flag
            org_var.table(j).exclude{ex_count} = varargin{i};
            ex_count = ex_count+1;
        end
        
    end
    
    
    % Create the "master-times" for each table
    output = create_sub_table(org_var.table(1), ex_end);
    
    % Combine each table into the final result
    output = create_master_table(org_var.table(1), ex_end);
    
    
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
    
    

end


function [mas_tab] = create_master_table(tab, ex_end)



end
