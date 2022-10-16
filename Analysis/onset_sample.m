function [on_mat, off_mat, un_mat] = onset_sample(mat, aff_table, aff_list, varargin)
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
%   band: [1-by-2 vector] The number of entries before and after the onset
%       to include in the returned matrix. Default is [0,0] for just the 
%       onset value.
%   offset: [bool] Whether to return a matrix of the ends of each affect.
%       Default is true.
%   omit_nan:[bool] Whether to omit rows which contain NaN values from both
%       on_mat and off_mat. This may result in on_mat and off_mat having
%       different numbers of rows. Default is false.

% Output:
%   on_mat: [o-by-m matrix] matrix of the onsets
%   off_mat: [f-by-m matrix] matrix of the offsets. An empty matrix if the
%       optional input 'offset' is false

    p = inputParser;
    addParameter(p, 'band', [0,0], @ismatrix);
    addParameter(p, 'offset', true, @islogical);
    addParameter(p, 'omit_nan', false, @islogical);
    addParameter(p, 'dilate', 0, @isscalar);
    parse(p,varargin{:});
    
    on_mat = [];
    off_mat = [];
    a = sum(p.Results.band);
    b = p.Results.band(1);
    c = p.Results.band(2);
    

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
    
    % TODO: Make vector of 1's or 0's then find the onsets, the current way
    % is not compatible with the "meta_chunk" added affect
    
    % ALSO THIS SAMPLING IS SUCCEPTABLE TO VERY SHORT AFFECTS, WHERE THE
    % END POINTS MY HAVE BINS WHICH SPAN THE ENTIRE AFFECT, DO NOT USE
    % OFFSETS WITHOUT THINKING THIS THROUGH
    
    aff_vec = affect_mark([], aff_table, aff_list);
    aff_vec = aff_vec;
    inst = find(aff_vec==1);
    vec = inst(2:end)-inst(1:end-1);
    idx = find(vec > a); % find starts greater in distance then the band being used so you don't
                         % accidently sample from another problematic
                         % behavior
    starts =[inst(1)];
    ends =[];
        
    for k=1:length(idx)
        starts = [starts,inst(idx(k)+1)];
        ends = [ends,inst(idx(k))];
    end
    ends = [ends, inst(end)]; %grab the last ending
    
    [on_mat, off_mat] = samp_(mat, starts, ends, p.Results.band);
    
    non_mat = [];
    noff_mat = [];
    
    nstarts = starts;
    nends = ends;
    while size(non_mat,1) < size(on_mat,1)
        
        [nstarts, nends] = dilate(nstarts, nends, p.Results.dilate, a);
        [dump_on, dump_off] = samp_(mat, nstarts, nends, p.Results.band);
        
        non_mat = [non_mat;dump_on];
        noff_mat = [noff_mat;dump_off];
    end
    un_mat = [non_mat; noff_mat];
    clear non_mat noff_mat

    %%
%     for i = 1:length(aff_table)
%         if any(strcmp(aff_list, aff_table(i,1)))
%             for j = 1:length(aff_table{i,2})
%                 % Append data to onset table
%                 on_mat(end+1:end+1+a,:) = mat(aff_table{i,2}(j)-b:aff_table{i,2}(j)+c,:);
%                 
%                 % If they want offset data
%                 if p.Results.offset
%                     off_mat(end+1:end+1+a,:) = mat(aff_table{i,3}(j)-b:aff_table{i,3}(j)+c,:);
%                 end
%                 
%                 %new_mat(aff_table{i,2}(j):aff_table{i,3}(j),end) = 1;    
%             end
%         end
%     end
    %%
    if p.Results.omit_nan
        on_mat(any(isnan(on_mat),2),:) = [];
        off_mat(any(isnan(off_mat),2),:) = [];
        un_mat(any(isnan(un_mat),2),:) = [];
    end

end

function [on_mat, off_mat] = samp_(mat, starts, ends, band)
    a = sum(band);
    b = band(1);
    c = band(2);
    
    lim = size(mat,1);
    
    on_mat=[];
    off_mat=[];
    
    for i = 1:length(starts)
        if starts(i) >= b && starts(i)+c <= lim
            on_mat(end+1:end+1+a,:) = mat(starts(i)-b:starts(i)+c,:);
        end
    end
    
    if ends
        for i = 1:length(ends)
            if ends(i)+c <= lim
                off_mat(end+1:end+1+a,:) = mat(ends(i)-b:ends(i)+c,:);
            end
        end
    end

end

function [dstarts, dends] = dilate(starts, ends, amnt, cap)
% Function which dilates the non-zero values of vector vec by the amount
% specified by amnt

    di_vec = [];

    d_starts = max(1,starts-amnt);
    d_ends = ends+amnt;
    
    for i = 1:length(d_starts)
        di_vec(d_starts(i):d_ends(i),1)=1;
    end
    
    
    inst = find(di_vec==1);
    vec = inst(2:end)-inst(1:end-1);
    idx = find(vec > cap); % find starts greater in distance then the band being used so you don't
                         % accidently sample from another problematic
                         % behavior
    dstarts = [inst(1)];
    dends = [];
    for k=1:length(idx)
        dstarts = [dstarts,inst(idx(k)+1)];
        dends = [dends,inst(idx(k))];
    end
    dends =[dends, inst(end)];
    
end
