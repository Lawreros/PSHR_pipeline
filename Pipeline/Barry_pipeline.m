function [] = Barry_pipeline(hr_files, ecg_files, varargin)
% Pipeline which analyzes the data collected by Barry and displays the
% RMSSD and pNN50 collected.
%
% Required Inputs:
%   hr_files: [1-by-n cell array] cell array containing a list of all of
%       the HR files to be included in the analysis
%   
%   ecg_files: [1-by-n cell array] cell array containing a list of all the
%       ECG files to be included in the analysis.
%
% Optional Inputs:
%   visualize: [bool] Whether to generate figures of  the data from the 
%   hr_files and ecg_files. Default is false.
%
%   individual_plots: [bool] Whether the generated figures are to be
%       seperate for each file, or put into a subplot. Default is false.

    p = inputParser;
    addRequired(p,'hr_files',@iscell);
    addRequired(p,'ecg_files',@iscell);
    addParameter(p,'visualize', true, @islogical);
    addParameter(p,'individual_plots',false, @islogical);
    parse(p,hr_files,ecg_files,varargin{:});

    % Load data from files
    Data = pshr_load('HR', hr_files, 'ECG', ecg_files, 'align', false);
    
    % Apply bandpass and acar filtering to the raw HR data in order to get rid of
    % outliers and compare model performance
    for i = 1:length(hr_files)
        Data.HR.PP.raw{i} = Data.HR.Raw{i}(:,[1:3]);
        Data.HR.PP.bandpass{i}(:,3) = bandpass(Data.HR.Raw{i}(:,3),300,1600,false);
        Data.HR.PP.bandpass{i}(:,1) = Data.HR.Raw{i}(:,1);
        
        Data.HR.PP.acar{i}(:,3) = acar(Data.HR.Raw{i}(:,3),5,false);
        Data.HR.PP.acar{i}(:,1) = Data.HR.Raw{i}(:,1);
    end

    % typ = list of fields in Data.HR.PP (in this case 'raw','bandpass' and
    % 'acar')
    typ = fieldnames(Data.HR.PP);
    
    % Iterate through all of the different preprocessing methods
    % to see how RMSSD and PNN50 are calculated
    for j = 1:length(hr_files)
        disp('=======================');
        for i=1:length(typ)
            ret = rmssd_calc(Data.HR.PP.(typ{i}){j}(:,3), false, false);
            disp(strcat(typ{i},' filtered ',Data.HR.files{j}));
            disp('RMSSD:' + string(ret));

            ret = pnnx_calc(Data.HR.PP.(typ{i}){j}(:,3), 50, false, false);
            disp('PNN50:' + string(ret));
            disp('Percent Removed:' + string(sum(isnan(Data.HR.PP.(typ{i}){j}(:,3))/size(Data.HR.PP.(typ{i}){j},1))));
        end
        disp('=======================');
    end
    
    % Visualize the collected data
    if p.Results.visualize
        % Depending on whether they want unique plots or a multiplot, generate
        % figure accordingly
        for i = 1:length(typ)
            if p.Results.individual_plots == false
                % Calculate dimensions of subplot
                row = floor(sqrt(length(Data.HR.PP.(typ{i}))));
                col = ceil(length(Data.HR.PP.(typ{i}))/row);
            end


            % Flag what data is present for plotting
            if isfield(Data, 'HR')
                % Plot HR data
                for j = 1:length(Data.HR.PP.(typ{i}))
                    if p.Results.individual_plots == 0
                        subplot(row, col, j);
                    else
                        figure;
                    end

                    [dump, f_nam] = fileparts(Data.HR.files{j});

                    plot(Data.HR.PP.(typ{i}){j}(:,1), Data.HR.PP.(typ{i}){j}(:,3));
                    xlabel('Timepoint (ms)');
                    title(f_nam);
                    ylabel('RR-interval (ms)');
                    
                end
                % Create title for entire subplot:
                sgtitle(strcat('Filtering type: ',typ{i}),'FontSize',30);

                % Plot HR Poincare
                figure;
                for j = 1:length(Data.HR.PP.(typ{i}))
                    if p.Results.individual_plots == 0
                        subplot(row, col, j);
                    else
                        figure;
                    end

                    [dump, f_nam] = fileparts(Data.HR.files{j});

                    [SD1, SD2] = poincare_plot(Data.HR.PP.(typ{i}){j}(:,3));
                    title([f_nam, strcat(' SD1: [',string(SD1),'] SD2: [',string(SD2),']')]);
                    % Create title for entire subplot:
                    sgtitle(strcat('Filtering type: ',typ{i}),'FontSize',30);
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

                    plot(Data.ECG.Raw{i}(:,3));
                    title(f_nam);
                end
            end
        end
    end
    
    % Test that PQRST works well on it
    if ~isempty(ecg_files)
        figure;
        ecg_PQRST(Data.ECG.Raw{1}(:,[1,3]), 'disp_plot', true);
    end
    
end