function [on_mat, off_mat, un_mat] = onset_sample(mat, keep_table, omit_table, varargin)
% This function takes a given matrix of values along with the affect start
% and end index table and returns a matrix consisting of a given selection
% of rows before/after each onset.

% Required Inputs:
%   mat:[n-by-m matrix] the given matrix you are sampling rows from
%
%   aff_table:[a-by-3 cell array] Found in the Data.(type).Affect structure
%       as an output of load_affect.m, contains the index values of the 
%       start and stop points for each affect
%
%   aff_list: [1-by-x cell array] list of affects from aff_table to use. If
%       false, then all affects present will be marked.

% Optional Inputs:
%   on_band: [1-by-2 matrix] The bounds, relative to the first timepoint of
%       a target affect instance, which define what counts as an "onset". For
%       example, inputting [3,2] would define onset as the first datapoint
%       of an instance of the target affect, the three datapoints
%       immediatly before the first datapoint, and the two datapoints after
%       the first datapoint. Default is [0,0].
% 
%   off_band: [1-by-2 matrix] The bounds, relative to the first timepoint of
%       a target affect instance, which define what counts as an "offset". For
%       example, inputting [3,2] would define offset as the last datapoint
%       of an instance of the target affect, the three datapoints
%       immediatly before the last datapoint, and the two datapoints after
%       the last datapoint. Default is [0,0].
%
%   omit_nan: [bool] Whether to omit rows which contain NaN values from both
%       on_mat and off_mat. This may result in on_mat and off_mat having
%       different numbers of rows. Default is false.
%
%   dilate: [int] Amount to dilate affect instance start/end times when
%       sampling for control values. Default is 5.

% Output:
%   on_mat: [o-by-m matrix] matrix of the onsets
%   off_mat: [f-by-m matrix] matrix of the offsets

    p = inputParser;
    addParameter(p, 'on_band', [0,0], @ismatrix);
    addParameter(p, 'off_band', [0,0], @ismatrix);
    addParameter(p, 'omit_nan', false, @islogical);
    addParameter(p, 'dilate', 5, @isscalar);
    parse(p,varargin{:});
    
    a = sum(p.Results.on_band);
    b = sum(p.Results.off_band);
    
    
    % ALSO THIS SAMPLING IS SUCCEPTABLE TO VERY SHORT AFFECTS, WHERE THE
    % END POINTS MY HAVE BINS WHICH SPAN THE ENTIRE AFFECT, DO NOT USE
    % OFFSETS WITHOUT THINKING THIS THROUGH
    
    aff_vec = affect_mark(zeros(size(mat,1),1), keep_table, false);
    aff_vec(:,1)=[];
    om_vec = affect_mark(zeros(size(mat,1),1), omit_table, false);
    om_vec(:,1)=[];
    inst = find(aff_vec==1);
    vec = inst(2:end)-inst(1:end-1);
    idx = find(vec > a); % find starts greater in distance then the band being used so you don't
                         % accidently sample from another problematic
                         % behavior
    starts =[];%;[inst(1)];
    ends =[];
    
    idx = [1;idx;length(inst)];
        
    for k=1:length(idx)-1
        % Check to make sure that you aren't pulling onsets from areas that
        % overlap with the omit times at any point
        if mean(om_vec(max(1,inst(idx(k)+1)-a):inst(idx(k)+1))) == 0 && ...
                mean(om_vec(inst(idx(k+1)):min(length(aff_vec),inst(idx(k+1))+b))) == 0
            
            if k == 1 %clunky fix for the first start value
                starts = [starts,inst(idx(k))];
            else
                starts = [starts,inst(idx(k)+1)];
            end
            
            ends = [ends,inst(idx(k+1))];
        end
    end
%     ends = [ends, inst(end)]; %grab the last ending
    
    [on_mat, off_mat] = samp_(mat, starts, ends, p.Results.on_band, p.Results.off_band);
    
    
    
    % Combine both the target and omits together and dilate from there to
    % get control measurements
    non_mat = [];
    noff_mat = [];
    
    inst = find(aff_vec + om_vec);
    vec = inst(2:end)-inst(1:end-1);
    idx = find(vec > a); 
    starts =[inst(1)];
    ends =[];
        
    for k=1:length(idx)
        starts = [starts,inst(idx(k)+1)];
        ends = [ends,inst(idx(k))];
    end
    ends = [ends, inst(end)]; %grab the last ending
    
    % "not-onset" and "not-offset" matrices to store the control values
    non_mat = [];
    noff_mat = [];
    
    nstarts = starts;
    nends = ends;
    
    
    if p.Results.omit_nan
        on_mat(any(isnan(on_mat),2),:) = [];
        off_mat(any(isnan(off_mat),2),:) = [];
    end
    
        
    % Calculate the amount of onset and/or offset measurements
    q = size(on_mat,1)+size(off_mat,1);
    cap = size(mat,1); % have maximum value to stop while loop from going forever
    
    while (size(non_mat,1) < q) && (nends(1) <= cap)
        
        [nstarts, nends] = dilate(nstarts, nends, [p.Results.dilate + randi(3),p.Results.dilate + randi(3)], a); % Add some randomness to dilation/sampling
        [dump_on, dump_off] = samp_(mat, nstarts, nends, p.Results.on_band, p.Results.off_band);
        
        non_mat = [non_mat;dump_on];
        noff_mat = [noff_mat;dump_off];
    end
    un_mat = [non_mat; noff_mat];
    clear non_mat noff_mat
    

    if p.Results.omit_nan
        un_mat(any(isnan(un_mat),2),:) = [];
    end

end

function [on_mat, off_mat] = samp_(mat, starts, ends, on_band, off_band)
    a = sum(on_band);
    b = on_band(1);
    c = on_band(2);
    
    lim = size(mat,1);
    
    on_mat=[];
    off_mat=[];
    
    for i = 1:length(starts)
        if starts(i) > b && starts(i)+c <= lim
            on_mat(end+1:end+1+a,:) = mat(starts(i)-b:starts(i)+c,:);
        end
    end
    
    
    a = sum(off_band);
    b = off_band(1);
    c = off_band(2);
    
    if ends
        for i = 1:length(ends)
            if ends(i)+c <= lim && ends(i)-b > 0
                off_mat(end+1:end+1+a,:) = mat(ends(i)-b:ends(i)+c,:);
            end
        end
    end

end

function [dstarts, dends] = dilate(starts, ends, amnt, cap)
% Function which dilates the non-zero values of vector vec by the amount
% specified by amnt. In other words, if you want to add 5 seconds before
% and after each affect instance (pair of start/stop timepoints), you call
% this function with amnt = 5
%
% Inputs:
%   starts: [1-by-n matrix] The start times for a given affect
%
%   ends: [1-by-n matrix] The end times for a given affect
%
%   amnt: [1-by-2 matrix] The amount to dialte the start and end times. This value
%       can either be positive or negative, with the first value being
%       subtracted from start times and the second being added to end
%       times.
%
%   cap: [int] The minimum amount of time between the start and stop of two
%       different instances of an affect for them to be considered seperate.
%       For example, if cap = 5, then if instance_1 and instance_2 of a given
%       affect are less than 5 seconds/timepoints apart they are lumped
%       together and considered one instance of the affect.
%
% Returns:
%   dstarts: [1-by-m matrix] The dilated start timepoints
%
%   dends: [1-by-m matrix] The dilated end timepoints

    di_vec = [];

    d_starts = max(1,starts-amnt(1));
    d_ends = ends+amnt(2);
    
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