function [on_mat, off_mat] = onset_sample(mat, aff_table, aff_list, varargin)
% This function takes a given matrix of values along with the affect start
% and end index table and returns a matrix consisting of a given selection
% of rows before/after each onset.

% Required Inputs:
%   mat:[n-by-m matrix] the given matrix you are sampling rows from
%   aff_table:[a-by-3 cell array] Found in the Data.(type).Affect structure
%       as an output of load_affect.m, contains the index values of the 
%       start and stop points for each affect
%   aff_list: [1-by-x cell array] list of affects from aff_table to use. If
%       false, then all affects present will be marked.

% Optional Inputs:
%   bound: [1-by-2 vector] The number of entries before and after the onset
%       to include in the returned matrix. Default is [0,0] for just the 
%       onset value.
%   offset: [bool] Whether to return a matrix of the ends of each affect.
%       Default is true.
%   omit_nan:[bool] Whether to omit rows which contain NaN values from both
%   on_mat and off_mat. This may result in on_mat and off_mat having
%   different numbers of rows. Default is false.

% Output:
%   on_mat: [o-by-m matrix] matrix of the onsets
%   off_mat: [o-by-m matrix] matrix of the offsets. An empty matrix if the
%       optional input 'offset' is false

    p = inputParser;
    addParameter(p, 'bond', [0,0], @ismatrix);
    addParameter(p, 'offset', true, @islogical);
    addParameter(p, 'omit_nan', false, @islogical);
    parse(p,varargin{:});
    
    on_mat = [];
    off_mat = [];
    a = sum(p.Results.bond);
    b = p.Results.bond(1);
    c = p.Results.bond(2);
    

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
                % Append data to onset table
                on_mat(end+1:end+1+a,:) = mat(aff_table{i,2}(j)-b:aff_table{i,2}(j)+c,:);
                
                % If they want offset data
                if p.Results.offset
                    off_mat(end+1:end+1+a,:) = mat(aff_table{i,3}(j)-b:aff_table{i,3}(j)+c,:);
                end
                
                %new_mat(aff_table{i,2}(j):aff_table{i,3}(j),end) = 1;    
            end
        end
    end
    if p.Results.omit_nan
        on_mat(any(isnan(on_mat),2),:) = [];
        off_mat(any(isnan(off_mat),2),:) = [];
    end

end