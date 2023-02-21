function [new_mat] = affect_mark(mat, aff_table, aff_list, varargin)
% This funciton takes a matrix of RR or ECG data along with the affect
% start and end index table and returns a matrix with an additional column
% marking a nonzero number (default 1) if any of the affects listed in aff_list 

% Input:
%   mat: [n-by-m matrix] which you want the affect data appended to as a
%       new column (m+1)
%
%   aff_table: [a-by-3 cell array] Found in the Data.(type).Affect structure
%       as an output of load_affect.m, contains the index values of the 
%       start and stop points for each affect.
%
%   aff_list: [1-by-x cell array] list of affects from aff_table to use. If
%       false, then all affects present will be marked, excluding the
%       following affects:
%       ' ' (an empty string, or string with single space)
%       'not problem'
%       'off camera'
%
%   NumberCategories: [bool] whether to assign each of the affects a different
%       number, instead of the binary 0 or 1. The number assigned will be
%       the index number for the affect in aff_list. DOES NOT WORK WITH
%       OVERLAPPING AFFECTS WHICH OCCUR AT THE SAME TIME, THE AFFECT ENTRY
%       WITH THE HIGHER ROW INDEX NUMBER FROM THE TABLE WILL TAKE PRIORITY.

% Output:
%   new_mat: [n-by-m+1 matrix]

    p = inputParser;
    addParameter(p, 'verbose', false, @islogical);
    addParameter(p, 'NumberCategories', false, @islogical);
    parse(p,varargin{:});
    
    new_mat = mat;
    new_mat(:,end+1)=0;
    
    if ~iscell(aff_list)
        % put all of aff_table's affects into aff_list, removing common
        % mistakes/errors
        aff_list = {};
        for i = 1:size(aff_table,1)
            if ~any(strcmp(aff_table{i},{' ', 'not problem', 'off camera'}))
                aff_list{end+1} = aff_table{i};
            end
        end
    end
    
    
    for i = 1:size(aff_table,1)
        if any(strcmp(aff_list, aff_table(i,1)))
            for j = 1:length(aff_table{i,2})
                % Mark anywhere where the affect occurs
                if p.Results.NumberCategories
                    new_mat(aff_table{i,2}(j):aff_table{i,3}(j),end) = find(strcmp(aff_table(i,1),aff_list));
                else
                    new_mat(aff_table{i,2}(j):aff_table{i,3}(j),end) = 1;
                end
            end
            
            if p.Results.verbose
                disp(strcat('Instance(s) found of : ',aff_table{i,1}));
            end
        end
    end
end


