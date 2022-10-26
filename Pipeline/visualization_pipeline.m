function [] = visualization_pipeline(varargin)
% Pipeline which focuses on data quality assurance through visual
% inspection. This pipeline serves to quickly and clearly display the
% collected data in a series of plots and simple metrics.

% Required Inputs:


% Optional Parameters:
%   hr_files: [1-by-n cell array] List of all the HR files you wish to
%       visualize. Default is {}.
%
%   ecg_files: [1-by-n cell array] List of all the ECG files you wish to
%       visualize. Default is {}.
%
%   aff_files: [1-by-n cell array] Providing affect files will results in
%       plots shaded by when certain affects occur. Default is {}.
%
%   aff_list: [1-by-n cell array] List of affects you wish to color on the
%       generated plots. Currently this function only works with binary
%       (problematic/non-problematic) coloring. If false, then all affects
%       (excluding `off camera` and `not problem` will be colored. Default is
%       false.
%
%   ignore_timestamps: [bool] Used when plotting data that isn't aligned,
%       due to the use of NaNs as times for the flattened data. When true
%       and aff_list = {}, then timestamps are replaced by index numbers.
%       Only relevant for RR-interval plots, as ECG has 1 non-nan timestamp
%       for every 50+ entries.
%       Default is false.
%
%   individual_plots: [bool] Whether to generate individual plots for each
%       file or to make figures containing multiple plots. Default is false.

%   break_points: [bool] Whether to put in breaks between each round
%       of figure generation in order to prevent too many plots appearing on
%       the screen at once. This prevents confusion about what plots belong together.
%       Default is false.

    aff_list = {'SIB','ISB','innappropriate face related behavior','polar strap adjustment/removal'...
        'repetitive behaviors','inappropriate movement','crying', 'pulling at pants'};

    p = inputParser;
    addParameter(p, 'hr_files', {}, @iscell);
    addParameter(p, 'ecg_files', {}, @iscell);
    addParameter(p, 'aff_files', {}, @iscell);
    addParameter(p, 'aff_list', false, @iscell);
    addParameter(p, 'ignore_timestamps', false, @islogical);
    addParameter(p, 'individual_plots', false, @islogical);
    parse(p,varargin{:});

    %If they provided affect files, assume that they want to align the
    %codings with the HR and ECG data when plotting
    align=false;
    if isempty(p.Results.aff_files)==0
        align = true;
    end
        

    %Load in data
    Data = pshr_load('HR', p.Results.hr_files, 'ECG', p.Results.ecg_files,...
        'Affect', p.Results.aff_files, 'align', align);
    
    % Depending on whether they want unique plots or a multiplot, generate
    % figure accordingly
    if p.Results.individual_plots == false
        % Calculate dimensions of subplot
        row = floor(sqrt(length(Data.HR.Raw)));
        col = ceil(length(Data.HR.Raw)/row);
    end
    
    
    % Flag what data is present for plotting
    if isfield(Data, 'HR')
        % Plot HR data
        for i = 1:length(Data.HR.Raw)
            if p.Results.individual_plots == 0
                subplot(row, col, i);
            else
                figure;
            end
            
            [dump, f_nam] = fileparts(Data.HR.files{i});

            if align && ~isempty(Data.HR.Affect{i})
                Data.HR.Raw{i} = affect_mark(Data.HR.Raw{i}, Data.HR.Affect{i}, false);
                colored_lineplot(Data.HR.Raw{i}(:,[1,3]), Data.HR.Raw{i}(:,4),...
                    'title', f_nam,'legend',false, 'fig_gen', false);
            else
                if p.Results.ignore_timestamps
                    plot(Data.HR.Raw{i}(:,3));
                    xlabel('Index number');
                else
                    plot(Data.HR.Raw{i}(:,1), Data.HR.Raw{i}(:,3));
                    xlabel('Timepoint (ms)');
                end
                title(f_nam);
                ylabel('RR-interval (ms)');
            end
        end
        
        % Plot HR Poincare
        figure;
        for i = 1:length(Data.HR.Raw)
            if p.Results.individual_plots == 0
                subplot(row, col, i);
            else
                figure;
            end
            
            [dump, f_nam] = fileparts(Data.HR.files{i});

            [SD1, SD2] = poincare_plot(Data.HR.Raw{i}(:,3));
            title([f_nam, strcat(' SD1: [',string(SD1),'] SD2: [',string(SD2),']')]);
        end
    end
    
    % Flag what data is present for plotting
    if isfield(Data, 'ECG')
        figure;
        % Plot ECG data
        for i = 1:length(Data.ECG.Raw)
            if p.Results.individual_plots == 0
                subplot(row, col, i);
            else
                figure;
            end
            
            [dump, f_nam] = fileparts(Data.ECG.files{i});

            if align && ~isempty(Data.ECG.Affect{i})
                Data.ECG.Raw{i} = affect_mark(Data.ECG.Raw{i}, Data.ECG.Affect{i}, aff_list);
                colored_lineplot(Data.ECG.Raw{i}(:,3), Data.ECG.Raw{i}(:,4),...
                    'title', f_nam, 'legend', false,'fig_gen', false);
            else
                plot(Data.ECG.Raw{i}(:,3));
                title(f_nam);
            end
        end
    end

end