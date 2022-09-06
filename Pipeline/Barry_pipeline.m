function [] = Barry_pipeline(hr_files, ecg_files)
% Pipeline which analyzes the data collected by Barry

    % Load data from files
    Data = pshr_load('HR', hr_files, 'ECG', ecg_files, 'align', false);
    
    % Create structure of data under different preprocessing
    % methods
    for i = 1:length(hr_files)
        prep.raw{i} = Data.HR.Raw{i}(:,3);
        prep.acar{i} = acar(Data.HR.Raw{i}(:,3), 5, false);
        prep.bandpass{i} = bandpass(Data.HR.Raw{i}(:,3),300,1300,false);
        prep.kamath{i} = kamath(Data.HR.Raw{i}(:,3), false);
        prep.karlsson{i} = karlsson(Data.HR.Raw{i}(:,3), false);
        prep.malik{i} = malik(Data.HR.Raw{i}(:,3), false);
    end

    typ = fieldnames(prep);
    
    % Iterate through all of the different preprocessing methods
    % to see how RMSSD and PNN50 are calculated
    for j = 1:length(hr_files)
        disp('=======================');
        for i=1:length(typ)
            ret = rmssd_calc(prep.(typ{i}){j}, false, false);
            disp(strcat(typ{i},' filtered ',Data.HR.files{j}));
            disp('RMSSD:' + string(ret));

            ret = pnnx_calc(prep.(typ{i}){j}, 50, false, false);
            disp('PNN50:' + string(ret));
            disp('Percent Removed:' + string(sum(isnan(prep.(typ{i}){j}))/length(prep.(typ{i}){j})));
        end
        disp('=======================');
    end
    
    % Visualize the collected data
    visualization_pipeline('hr_files', hr_files, 'ecg_files', ecg_files, ...
        'ignore_timestamps', true, 'individual_plots', false);
    
    
    % Test that PQRST works well on it
    if ~isempty(ecg_files)
        figure;
        ecg_PQRST(Data.ECG.Raw{1}(:,[1,3]), 'disp_plot', true);
    end
    
end