function [new_mat] = affect_mark(mat, aff_table, aff_list)
% This funciton takes a matrix of RR or ECG data along with the affect
% start and end index table and returns a matrix with an additional column
% marking a 1 if any of the affects listed in aff_list 

% Input:
%   mat: [n-by-m matrix] which you want the affect data appended to
%   aff_table: [a-by-3 cell array] Found in the Data.(type).Affect structure
%       as an output of load_affect.m, contains the index values of the 
%       start and stop points for each affect
%   aff_list: [1-by-x cell array] list of affects from aff_table to use. If
%   false, then all affects present will be marked.

% Output:
%   new_mat: [n-by-m+1 matrix]

    new_mat = mat;
    new_mat(:,end+1)=0;
    
    if ~iscell(aff_list)
        % put all of aff_table's affects into aff_list, removing common
        % mistakes/errors
        aff_list = {};
        for i = 1:length(aff_table)
            if ~any(strcmp(aff_table{i},{' ', 'not problem', 'off camera'}))
                aff_list{end+1} = aff_table{i};
            end
        end
    end
    
    
    for i = 1:length(aff_table)
        if any(strcmp(aff_list, aff_table(i,1)))
            for j = 1:length(aff_table{i,2})
                % Mark anywhere where 
                new_mat(aff_table{i,2}(j):aff_table{i,3}(j),end) = 1;    
            end
            
            if j
                disp(strcat('Instance(s) found of : ',aff_table{i,1}));
            end
        end
    end
end


