function [] = colored_lineplot(mat, label,varargin)
% Function which takes a given matrix and creates a 2D plot where each
% category of label is given a different color and displayed as a shaded
% region over a lineplot.

% Inputs:
%   mat: [n-by-1 or n-by-2 matrix] 

%   label: [n-by-1 vector] either a vector of numbers or a cell array of
%   strings, labeling each datapoint. This function will take the input
%   values and convert them into categories, where it then asigns each
%   category a color

%   title: [string] Title of the figure you wish to create. Default is
%   'Shaded Plot'.
%   xlabel: [string] Label put on the x-axis of the figure you wish to
%   create. Default is 'x-axis'.
%   ylabel: [string] Label put on the y-axis of the figure you wish to
%   create. Default is 'y-axis'.
%
%   legend: [bool] Super hacky, but writes text onto the figure to label
%   what each of the shaded regions represents. Default is true.
%
%   fig_gen: [bool] Whether to plot the results as a new figure. If false,
%   then this function can work with subplot. Default is true.

    %Input parser
    p = inputParser;
    addParameter(p,'title','Shaded Plot',@ischar);
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
    start_idx = find([0;lab_idx]-[lab_idx;0]); % Finds points where the label changes
    start_idx(end) = length(lab_idx);
    
    % Convert the unique value locations into coordinates
    patch_x = [mat(start_idx(1:end-1),1),mat(start_idx(1:end-1),1),...
        mat(start_idx(2:end),1), mat(start_idx(2:end),1)].';
    patch_y = zeros(size(patch_x));
    patch_y(2:3,:) = 2000;
    
    % Construct color grid for plotting:
    face_colors = zeros(length(start_idx)-1,1,3);
    
    color_key = [0,0,1; % Color key for the different categories
        0,1,0;
        0,0,1;
        0,1,1;
        1,1,0];
    
    q = 1;
    for i = double(lab_idx(start_idx(1:end-1))).' %for loop trips over itself when iterating through column
        face_colors(q,:,:) = color_key(i,:); % use specified color_key for each category
        q = q+1;
    end
    
    pa = patch(patch_x, patch_y,'w');
    pa.CData = face_colors; % n-by-1-by-3 array of colors for each patch
    pa.FaceColor = 'flat';
    pa.FaceAlpha = .2;
    
%     if true % Whether to make a legend for the plot
%         for i = 1:length(unique(categorical(label)))
%             % For each unique label category, color one entry on legend
%             leg{i} = plot(NaN, color_key(i,:), 'Linewidth', 8);
%         end
%         %https://www.mathworks.com/matlabcentral/answers/1626265-create-a-custom-legend
%         %https://www.mathworks.com/matlabcentral/answers/388575-convert-string-array-into-cell-array
%         legend([leg{:}],cellstr(string(unique(label))),'location','best');
%     end


    % Somewhat hacky way of getting a custom legend into the figure, but I
    % really didn't have much of a choice because of how much MATLAB's
    % legend function sucks :(
    % https://www.delftstack.com/howto/matlab/matlab-custom-legend/#:~:text=Add%20Custom%20Legends%20Using%20the%20text()%20Function%20in%20MATLAB,-We%20can%20also&text=Simply%20plot%20the%20variable%20and,able%20to%20see%20the%20text.
    
    if p.Results.legend
        q = string(unique(categorical(label)));
        for i = 1:length(q)
            text(0,i*50,q(i),'Fontsize', 18, 'Color', color_key(i,:));
        end
    end


    hold off;
    
end