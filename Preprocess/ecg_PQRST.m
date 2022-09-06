function [locations_matrix, time_matrix] = ecg_PQRST(mat, varargin)
% Function which takes an ECG signal and returns a table containing the
% index location of each of the 
% Required Inputs:
%   mat: [n-by-2 matrix] vector containing the ECG data you wish to
%       process, along with corresponding timestamps. This has ideally has 
%       already been preprocessed and can contain NaN values. 
%       Rows of the input matrix are in the form of [timestamp, ECG_value].
%       If you do not have/care about timestamps, have the first column be
%       all 1's (or any number).
%       

% Optional Parameters:
%       varargin: [cell array] Contains the various input arguments
%       (excluding mat) for parsing. Due to the large quantity of
%       parameters, this allows for default values to be used for every
%       argument that isn't passed to this function. To pass these input
%       argurment values into the function, input them as ('name of
%       parameter1', value, 'name of parameter2', value, etc.).
%       These parameters are used for the "findpeaks" functions 
%       (see mathworks documentation):
%
%           *_MinPeakProminence: [int]
%           *_MinPeakDistance: [int]
%           *_MinPeakWidth: [int]
%
%       where * represents either the P, QS, R, or T waves. For the default
%       values of these, see the code before "%% Start of analysis code".
%       

%   disp_plot: [bool] Whether to plot the data in mat with the
%       different waves indicated in the figure. Default is false.
%
%   bin_width: [int] The maximum amount of indices before and after
%       the R-wave to look for the other waves in the PQRST complex.
%       Default is 50.
%
%   sample_rate: [int] The sampling rate (Hz) for the ECG data being
%       analyzed. If provided, a matrix converting the information in 
%       locations_matrix will be made in terms of milliseconds. 
%       With the R-wave column converted into the RR-interval duration 
%       and the other waves converted into milliseconds before/after the R-wave. 
%       Default value is false, which does not run the process and sets
%       time_matrix to false.
%
%   min_waves: [int or 1-by-5 vector] Whether to cut all R-wave
%       locations and their associated P, Q, S, and T-wave rows from
%       the returned locations_matrix. Can either specify how many of
%       the wave components must be present to keep the R-wave
%       (inputing '3' will mean that rows containing R-waves with 3 or 
%       more wave components, including the R-wave, will be kept) or
%       the specific waves that need to be present (inputing a binary
%       vector for which waves have to be present with the template of:
%       [P,Q,R,S,T]. Thus requiring Q and T-waves be present is
%       [0,1,1,0,1]. Default value is false, which does not activate
%       this filtration process.
%           

% Returns:
%   locations_matrix: [m-by-6 matrix] matrix containing the indices of
%       the P, Q, S, and T-waves associated with each of the R-waves. Will
%       contain NaN values if no wave is found. The output format has the
%       columns: [time_stamp, P, Q, R, S, T]

%   time_matrix: [m-by-6 matrix] matrix contining the relative times
%       for the P, Q, S, and T-waves associated with each of the R-waves
%       (converted into RR-intervals). This is created if sample_rate is
%       provided a non-false integer, resulting in units of milliseconds.
%       The output format has the columns: 
%           [time_stamp, P-R, Q-R, RR-interval, S-R, T-R]
%       where *-R represents the difference in milliseconds between the *
%       and R wave. RR-interval is the amount of time that has passed since
%       the last R-wave.


    % Input parser so you don't have to define a bunch of parameters, only
    % the ones that you need.
    p = inputParser;
    p.KeepUnmatched=true;
    validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0); %specify type of input that is a positive scalar
    %addRequired(p,'mat',@ismatrix);
    %addRequired(p,'width',validScalarPosNum);
    addParameter(p,'R_MinPeakProminence',800,validScalarPosNum);
    addParameter(p,'R_MinPeakDistance',40,validScalarPosNum);

    addParameter(p,'QS_MinPeakProminence',100,validScalarPosNum);
    addParameter(p,'QS_MinPeakDistance',5,validScalarPosNum);

    addParameter(p,'T_MinPeakProminence',50,validScalarPosNum);
    addParameter(p,'T_MinPeakWidth',10,validScalarPosNum);
    addParameter(p,'T_MinPeakDistance',80,validScalarPosNum);

    addParameter(p,'P_MinPeakProminence',8,validScalarPosNum);
    addParameter(p,'P_MinPeakWidth',20,validScalarPosNum);
    addParameter(p,'P_MinPeakDistance',20,validScalarPosNum);

    %addOptional(p,'height',defaultHeight,validScalarPosNum);
    addParameter(p,'bin_width',50,validScalarPosNum);
    addParameter(p,'sample_rate',false,@isnumeric);
    addParameter(p,'disp_plot',false,@islogical);
    addParameter(p,'min_waves', false, @isvector);
    %addParameter(p,'shape',defaultShape, @(x) any(validatestring(x,expectedShapes)));
    %parse(p,mat,varargin{:});
    parse(p,varargin{:});
    
    
    %% Start of analysis code
    
    [pks_R, locs_R] = findpeaks(mat(:,2), 'MinPeakProminence', p.Results.R_MinPeakProminence, 'MinPeakDistance', p.Results.R_MinPeakDistance);

    % From there find the location of the Q and S waves (which should be very
    % close to the R-wave peak), this is done by multiplying section by -1 and
    % then finding the peaks on that, since it will ignore the RR-intervals
    % https://www.mathworks.com/matlabcentral/answers/267303-how-to-find-q-and-s-point-in-qrs-complex-of-ecg-signal


    %This finds Q and S fairly well, but picks up part of T and other noise,
    %add a step to this to only take the immediate surrounding 2 points to the
    %R-wave
    [pks_QS, locs_QS] = findpeaks(-mat(:,2), 'MinPeakProminence', p.Results.QS_MinPeakProminence, 'MinPeakDistance', p.Results.QS_MinPeakDistance);

    % Attmpt to find P and T waves, which are their own peaks which occur
    % before and after the Q and S waves, respectively

    [pks_T, locs_T] = findpeaks(mat(:,2),'MinPeakWidth',p.Results.T_MinPeakWidth, 'MinPeakProminence', p.Results.T_MinPeakProminence,'MinPeakDistance',p.Results.T_MinPeakDistance);


    [pks_P, locs_P] = findpeaks(mat(:,2),'MinPeakWidth', p.Results.P_MinPeakWidth,'MinPeakProminence', p.Results.P_MinPeakProminence,'MinPeakDistance', p.Results.P_MinPeakDistance);


    if p.Results.disp_plot
        plot([1:length(mat(:,2))],mat(:,2),locs_R,pks_R,'^b',locs_QS,-pks_QS,'vg', locs_T, pks_T,'^r',locs_P, pks_P, 'vr');
        legend('EKG', 'R', 'Q/S', 'T','P');
        title("Results from ecg_PQRST analysis");
        ylabel("Voltage (uV)");
    end

    %% Organize collections of PQRST, using the R-wave as a base:

    align_matrix = zeros(length(mat(:,2)),4);

    align_matrix(locs_P, 1) = 1; % Column of binary indicators for the P-wave locations
    align_matrix(locs_QS, 2) = 1; % Q and S wave locations, due to them being picked up by the same findpeaks pass
    align_matrix(locs_R,3) =1; % R-wave
    align_matrix(locs_T,4) = 1; % T-wave

    locations_matrix = NaN(length(locs_R), 6);
    locations_matrix(:,4) = locs_R; % Already know the answer for the R-waves

    max_dat = length(align_matrix(:,1));
    
    for i = 1:length(locs_R) % cycle through all of the R-waves and find the other waves

        % Go backward from the R-wave to find the P and Q wave locations
        p_check = 0;
        q_check = 0;
        for j = locations_matrix(i,4)-1 : -1 : locations_matrix(i,4) - p.Results.bin_width
            
            if j == 0
                break; %prevent from looking outside of bounds
            end
            
            if align_matrix(j,1)==1 && p_check==0
                locations_matrix(i,2) = j;
                p_check=1;
            elseif align_matrix(j,2)==1 && q_check==0
                locations_matrix(i,3) = j;
                q_check=1;
            elseif p_check==1 && q_check==1 %no need to keep checking if you've already filled the spots
                break
            end
        end

        % Go forward from the R-wave to find the S and T wave locations
        s_check = 0;
        t_check = 0;
        for j = locations_matrix(i,4)+1: locations_matrix(i,4) + p.Results.bin_width
            
            if j > max_dat
                break;
            end
            
            
            if align_matrix(j,2)==1 && s_check==0
                locations_matrix(i,5) = j;
                s_check = 1;
            elseif align_matrix(j,4)==1 && t_check==0
                locations_matrix(i,6) = j;
                t_check = 1;
            elseif s_check==1 && t_check==1
                break
            end
        end
    end

    %% Return only a list of RR-intervals which meet a specific specification

    %min_waves = [0,1,1,1,0];
    %min_waves = 3;
    min_waves = p.Results.min_waves;

    if islogical(min_waves)==0 % If they want to do the wave filtration
        cut = [];
        
        if isscalar(min_waves) %If they have specified a minimum amount of waves
            for i = 1:length(locations_matrix)
                % Find how many rows contain too many NaNs
                if sum(isnan(locations_matrix(i,2:end))) > (5 - min_waves) %if there are more NaNs than (max - scalar)
                    cut = [cut,i];
                end
            end
        else %If they have input the vector of which waves MUST be present
            
            % First just cut all of the columns not included in analysis
            dump = locations_matrix(:,2:end);
            dump(:,find(min_waves==0))=[];

            % Go through each row and catalog which RR-intervals fail
            for i = 1:length(dump)
                if sum(isnan(dump(i,:))) > 0
                    % A NaN has been found
                    cut = [cut,i];
                end
            end

        end

        % Cut the rows from locations_matrix which fail
        locations_matrix(cut,:) = [];
    end
    
    clear dump locs* pks*; %just some cleaning up of variables
    
    %% Get timestamps for the locations of each RR-interval
    
    for i=1:length(locations_matrix(:,4)) % Go through each R-wave index and find the closest corresponding timestamp in mat
        j = 0;
        
        while isnan(mat(locations_matrix(i,4)-j,1))
            j = j+1;
        end
        locations_matrix(i,1) = mat(locations_matrix(i,4)-j,1);
    end
    
    
    %% Return the PQRST matrix with everything converted to milliseconds
    
    % This option will convert the index locations into
    time_matrix = false;
    if p.Results.sample_rate %If they have specified wanting relative timepoints instead of indexes
        smp = 1000 / p.Results.sample_rate; %multiply by indexes to convert into milliseconds
        time_matrix = nan(size(locations_matrix));
        time_matrix(:,1) = locations_matrix(:,1);
        
        %Because the first entry is a little weird, just do it outside of
        %the for loop
        time_matrix(1,4) = locations_matrix(1,4)*smp;
        time_matrix(1,[2,3,5,6]) = (locations_matrix(1,[2,3,5,6])-locations_matrix(1,4))*smp;
        
        for i=2:length(locations_matrix(:,4)) %skip the first entry
            time_matrix(i,4) = (locations_matrix(i,4) - locations_matrix(i-1,4))*smp;
            time_matrix(i,[2,3,5,6]) = (locations_matrix(i,[2,3,5,6])-locations_matrix(i,4))*smp;
        end
    end
    
    
end