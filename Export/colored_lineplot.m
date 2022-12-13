function [] = colored_lineplot(mat, label,varargin)
% Function which takes a given matrix and creates a 2D plot where each
% category of label is given a different color and displayed as a shaded
% region over a lineplot.

% Required Inputs:
%   mat: [n-by-1 or n-by-2 matrix] 

% Optional Parameters:
%   label: [n-by-1 vector] either a vector of numbers or a cell array of
%       strings, labeling each datapoint. This function will take the input
%       values and convert them into categories, where it then asigns each
%       category a color

%   title: [string] Title of the figure you wish to create. Default is
%       'Shaded Plot'.
%
%   NumCategories: [int] The number of unique categories present across all
%       of the recordings you wish to plot. This is so there is consistent
%       coloring across shading. If this is 0 then the number of unique
%       categories is calculated from the provided label vector. Default is
%       0.
%
%   xlabel: [string] Label put on the x-axis of the figure you wish to
%       create. Default is 'x-axis'.

%   ylabel: [string] Label put on the y-axis of the figure you wish to
%       create. Default is 'y-axis'.
%
%   legend: [bool] Super hacky, but writes text onto the figure to label
%       what each of the shaded regions represents. Default is true.
%
%   fig_gen: [bool] Whether to plot the results as a new figure. If false,
%       then this function can work with subplot. Default is true.

    %Input parser
    p = inputParser;
    addParameter(p,'title','Shaded Plot',@ischar);
    addParameter(p,'NumCategories',0,@isnumeric);
    addParameter(p,'xlabel','x-axis', @ischar);
    addParameter(p,'ylabel','y-axis', @ischar);
    addParameter(p,'legend',true,@islogical);
    addParameter(p,'fig_gen',true,@islogical);
    
    parse(p,varargin{:});


    % Check the dimensions of mat to see if the x-axis has unique values
    [r,c] = size(mat);

    if c == 1 % Generate index values for the x-axis
        mat = [[1:r].',mat];
    end
    
    % Plot the data in mat
    if p.Results.fig_gen
        figure;
    end
    plot(mat(:,1),mat(:,2));
    title(p.Results.title);
    xlabel(p.Results.xlabel);
    ylabel(p.Results.ylabel);
    hold on;


    % Convert label vector into usable values
    lab_idx = grp2idx(categorical(label)); % Converts label into categories
    start_idx = find([99;lab_idx]-[lab_idx;99]); % Finds points where the label changes
    start_idx(end) = length(lab_idx);
    
    % Check that, if the last affect goes unil the end, the last value is
    % not NaN
    q = 0;
    while isnan(mat(start_idx(end)-q,1))
        q = q + 1;
    end
    start_idx(end) = start_idx(end)-q;
    
    % Convert the unique value locations into coordinates
    patch_x = [mat(start_idx(1:end-1),1),mat(start_idx(1:end-1),1),...
        mat(start_idx(2:end),1), mat(start_idx(2:end),1)].';
    patch_y = zeros(size(patch_x))-500;
    patch_y(2:3,:) = 4000;
    
    % Construct color grid for plotting:
    face_colors = zeros(length(start_idx)-1,1,3);
    
    if p.Results.NumCategories == 0
        color_key = turbo(length(unique(lab_idx))); %Generate a color key for the different categories
    else
        color_key = turbo(p.Results.NumCategories + 1); % 1 added for the coloring of 0
    end
    
    q = 1;
    k_ = unique(label);
    for i = double(lab_idx(start_idx(1:end-1))).' %for loop trips over itself when iterating through column
        face_colors(q,:,:) = color_key(k_(i)+1,:); % use specified color_key for each category
        q = q+1;
    end
    
    pa = patch(patch_x, patch_y,'w');
    pa.CData = face_colors; % n-by-1-by-3 array of colors for each patch
    pa.FaceColor = 'flat';
    pa.FaceAlpha = .2;
    
    % Somewhat hacky way of getting a custom legend into the figure, but I
    % really didn't have much of a choice because of how much MATLAB's
    % legend function sucks :(
    % https://www.mathworks.com/matlabcentral/answers/231531-legend-for-a-patches-object
    
    if p.Results.legend
        label = num2str(k_);
        %hold on;
        for il = 1:length(k_)
            hl(il) = patch(NaN,NaN, color_key(k_(il)+1,:),'FaceAlpha',.2);
        end
        legend(hl,label);
    end


    hold off;
    
end