function [aligned, move] = levenshtein_align(mat_1, col_1, mat_2, col_2, subcost, time_pen)
% Functions which attempts to align to matrices by the specified column
% using the Levenshtein distance metric:
% Inputs:
%   mat_1: [n-by-m matrix] The matrix to which you are aligning mat_2. This
%       matrix will not be altered in any way during the alignment process, and
%       thus is not returned by this function.
%
%   col_1: [int] What column to use in mat_1 for the alignment process,
%       this will be aligned to mat_2(:,col_2).
%
%   mat_2: [a-by-b matrix] The matrix which you are aligning to mat_1. This
%       matrix will have rows deleted or new rows of NaN values created as a
%       result of this alignment process.
%
%   col_2: [int] What column to use in mat_2 for the alignment process,
%       this will be aligned to mat_1(:,col_1).
%
%   subcost: [int] substitution cost coefficient for the use in Levenshtein
%       calculation. The equation of cost is: abs(X-Y)/subcost. Increasing
%       subcost will weight against substitution/deletion. Recommended value is
%       10.
%
%   time_pen: [int] denominator to divide the absolute difference between
%       the two columns being aligned as an added weight to subcost. The equation
%       for cost becomes: abs(X-Y)/subcost + abs(X-Y)/time_pen. If 0 then
%       this penality method will not be included with the subcost calculation.

% Returns:
%   aligned: [n-by-b matrix] the resulting contents of mat_2 after being 
%       aligned to mat_1
%   move: [n-by-1 matrix] a vector of the decided additions/deletions for
%       mat_2 in order for it to be aligned with mat_1


    % Isolate columns of interest so I don't have to keep typing out
    % (:,col*)
    va = mat_1;%mat_1(:,col_1);
    vb = mat_2(:,[1,4]);%mat_2(:,col_2);
    
    % Create matrix for alignment values for levenshtein distance
    lev = zeros(size(va,1),size(vb,1));
    lev(1,1) = 0;

    for i = 2:size(va,1)
        lev(i,1) = i-1;
    end
    for j = 2:size(vb,1)
        lev(1,j) = j-1;
    end

    % Go through the combinations of values and calculate the
    % substitution costs between values
    if time_pen
        for j = 2:size(vb,1)
            for i = 2:size(va,1)
                if va(i,2) == vb(j,2)
                    substitutionCost = 0;
                else
                    substitutionCost = (abs(va(i,2) - vb(j,2))/subcost)+(abs(va(i,1)-vb(j,1))/time_pen);
                end
                lev(i,j) = min([lev(i-1,j)+1, lev(i,j-1)+1, lev(i-1, j-1)+substitutionCost]);
            end
        end
    else
        for j = 2:size(vb,1)
            for i = 2:size(va,1)
                if va(i,2) == vb(j,2)
                    substitutionCost = 0;
                else
                    substitutionCost = (abs(va(i,2) - vb(j,2))/subcost);
                end
                lev(i,j) = min([lev(i-1,j)+1, lev(i,j-1)+1, lev(i-1, j-1)+substitutionCost]);
            end
        end
    end

    [Y,I] = min(lev(end,:));

    [r,c] = size(lev);
    c = I;

    % Create the vector move which will store the movements that the
    % computer follows along the lev matrix
    move = [];
    q_c = 0;
    q_r = 0;

    while r > 0

        if c>1 && r>1
            q = min([lev(r,c-1),lev(r-1,c-1), lev(r-1,c)]);
        elseif r == 1 && c~=1
            q = lev(r,c-1);
            q_r = 1;
        elseif c == 1 && r ~=1
            q = lev(r-1,c);
            q_c = 1;
        end

        if r == 1
            move(end+1:end+c-1) = 1;
            r= 0;
        elseif c == 1
            move(end+1:end+r-1) = 2;
            r = 0;
        elseif q == lev(r-1,c-1)
            move=[move,0]; %diag
            r=r-1;
            c=c-1;
        elseif q == lev(r,c-1) && q_c ~=1
            move=[move,1]; %left
            c=c-1;
        elseif q == lev(r-1,c) && q_r ~=1
            move=[move,2]; %up
            r = r-1;
        elseif q_c == 1 && q_r == 1
            r=r-1;
            c=c-1;
        end

    end
    move = [move,0];

    % Flip move for path to maximize alignment
    move = flip(move,2);

    % Clean up a bunch of variables that were made during this process
    clear i j c r I line pks_a q q_c q_r Y substitutionCost;

    % Realign the two matrices 
    j = 1;
    aligned=[];
    
    for i = 1:length(move)
        if move(1,i) == 0
            aligned(end+1,:) = mat_2(j,:);
            j = j+1;
        elseif move(1,i) == 1
            j = j+1;
        elseif move(1,i) == 2
            aligned(end+1,:) = NaN;
        end

    end
    
end