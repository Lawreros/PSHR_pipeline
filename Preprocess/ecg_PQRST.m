function [locations_matrix] = ecg_PQRST(mat, varargin)
% Function which takes an ECG signal and returns a table containing the
% index location of each of the 
%   Inputs:
%       mat: [n-by-1 vector] vector containing the ECG data you wish to
%       process. This has ideally already been preprocessed and is robust
%       to the presence of NaN values.
%       
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
%       There are also the input parameters:
%           disp_plot: [bool] Whether to plot the data in mat with the
%           different waves indicated in the figure. Default is false.
%
%           bin_width: [int] The maximum amount of indices before and after
%           the R-wave to look for the other waves in the PQRST complex.
%           Default is 50.
%
%           min_waves: [int or 1-by-5 vector] Whether to cut all R-wave
%           locations and their associated P, Q, S, and T-wave rows from
%           the returned locations_matrix. Can either specify how many of
%           the wave components must be present to keep the R-wave
%           (inputing '3' will mean that rows containing R-waves with 3 or 
%           more wave components, including the R-wave, will be kept) or
%           the specific waves that need to be present (inputing a binary
%           vector for which waves have to be present with the template of:
%           [P,Q,R,S,T]. Thus requiring Q and T-waves be present is
%           [0,1,1,0,1]. Default value is false, which does not activate
%           this filtration process.
%           

%   Returns:
%       locations_matrix: [m-by-5 matrix] matrix containing the indices of
%       the P, Q, S, and T-waves associated with each of the R-waves. Will
%       contain NaN values if no wave is found. The output format has the
%       columns: [P, Q, R, S, T]



    %% Add input parser due to the quantity of values that can be adjusted
    R_peakp = 800;
    R_peakd = 40;
    
    QS_peakp = 100;
    QS_peakd = 5;

    T_peakp = 50;
    T_peakw = 10;
    T_peakd = 80;

    P_peakp = 8;
    P_peakw = 20;
    P_peakd = 20;
    
    bin_width = 50;

    % Input parser so you don't have to define a bunch of parameters, only
    % the ones that you need.
    p = inputParser;
    validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0); %specify type of input that is a positive scalar
    %addRequired(p,'mat',@ismatrix);
    %addRequired(p,'width',validScalarPosNum);
    addParameter(p,'R_MinPeakProminence',R_peakp,validScalarPosNum);
    addParameter(p,'R_MinPeakDistance',R_peakd,validScalarPosNum);

    addParameter(p,'QS_MinPeakProminence',QS_peakp,validScalarPosNum);
    addParameter(p,'QS_MinPeakDistance',QS_peakd,validScalarPosNum);

    addParameter(p,'T_MinPeakProminence',T_peakp,validScalarPosNum);
    addParameter(p,'T_MinPeakWidth',T_peakw,validScalarPosNum);
    addParameter(p,'T_MinPeakDistance',T_peakd,validScalarPosNum);

    addParameter(p,'P_MinPeakProminence',P_peakp,validScalarPosNum);
    addParameter(p,'P_MinPeakWidth',P_peakw,validScalarPosNum);
    addParameter(p,'P_MinPeakDistance',P_peakd,validScalarPosNum);

    %addOptional(p,'height',defaultHeight,validScalarPosNum);
    addParameter(p,'bin_width',bin_width,validScalarPosNum);
    addParameter(p,'disp_plot',false,@islogical);
    addParameter(p,'min_waves', false, @isvector);
    %addParameter(p,'shape',defaultShape, @(x) any(validatestring(x,expectedShapes)));
    %parse(p,mat,varargin{:});
    parse(p,varargin{:});
    
    
    %% Start of analysis code
    
    [pks_R, locs_R] = findpeaks(mat, 'MinPeakProminence', p.Results.R_MinPeakProminence, 'MinPeakDistance', p.Results.R_MinPeakDistance);

    % From there find the location of the Q and S waves (which should be very
    % close to the R-wave peak), this is done by multiplying section by -1 and
    % then finding the peaks on that, since it will ignore the RR-intervals
    % https://www.mathworks.com/matlabcentral/answers/267303-how-to-find-q-and-s-point-in-qrs-complex-of-ecg-signal


    %This finds Q and S fairly well, but picks up part of T and other noise,
    %add a step to this to only take the immediate surrounding 2 points to the
    %R-wave
    [pks_QS, locs_QS] = findpeaks(-mat, 'MinPeakProminence', p.Results.QS_MinPeakProminence, 'MinPeakDistance', p.Results.QS_MinPeakDistance);

    % Attmpt to find P and T waves, which are their own peaks which occur
    % before and after the Q and S waves, respectively

    [pks_T, locs_T] = findpeaks(mat,'MinPeakWidth',p.Results.T_MinPeakWidth, 'MinPeakProminence', p.Results.T_MinPeakProminence,'MinPeakDistance',p.Results.T_MinPeakDistance);


    [pks_P, locs_P] = findpeaks(mat,'MinPeakWidth', p.Results.P_MinPeakWidth,'MinPeakProminence', p.Results.P_MinPeakProminence,'MinPeakDistance', p.Results.P_MinPeakDistance);


    if p.Results.disp_plot
        plot([1:length(mat)],mat,locs_R,pks_R,'^b',locs_QS,-pks_QS,'vg', locs_T, pks_T,'^r',locs_P, pks_P, 'vr');
        legend('EKG', 'R', 'Q/S', 'T','P');
        title("Results from ecg_PQRST analysis");
    end

    %% Organize collections of PQRST, using the R-wave as a base:

    align_matrix = zeros(length(mat),4);

    align_matrix(locs_P, 1) = 1; % Column of binary indicators for the P-wave locations
    align_matrix(locs_QS, 2) = 1; % Q and S wave locations, due to them being picked up by the same findpeaks pass
    align_matrix(locs_R,3) =1; % R-wave
    align_matrix(locs_T,4) = 1; % T-wave

    locations_matrix = NaN(length(locs_R), 5);
    locations_matrix(:,3) = locs_R; % Already know the answer for the R-waves

    for i = 1:length(locs_R) % cycle through all of the R-waves and find the other waves

        % Go backward from the R-wave to find the P and Q wave locations
        p_check = 0;
        q_check = 0;
        for j = locations_matrix(i,3)-1 : -1 : locations_matrix(i,3) - p.Results.bin_width
            if align_matrix(j,1)==1 && p_check==0
                locations_matrix(i,1) = j;
                p_check=1;
            elseif align_matrix(j,2)==1 && q_check==0
                locations_matrix(i,2) = j;
                q_check=1;
            elseif p_check==1 && q_check==1 %no need to keep checking if you've already filled the spots
                break
            end
        end

        % Go forward from the R-wave to find the S and T wave locations
        s_check = 0;
        t_check = 0;
        for j = locations_matrix(i,3)+1: locations_matrix(i,3) + p.Results.bin_width
            if align_matrix(j,2)==1 && s_check==0
                locations_matrix(i,4) = j;
                s_check = 1;
            elseif align_matrix(j,4)==1 && t_check==0
                locations_matrix(i,5) = j;
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
                if sum(isnan(locations_matrix(i,:))) > (5 - min_waves) %if there are more NaNs than (max - scalar)
                    cut = [cut,i];
                end
            end
        else %If they have input the vector of which waves MUST be present
            
            % First just cut all of the columns not included in analysis
            mat = locations_matrix;
            mat(:,find(min_waves==0))=[];

            % Go through each row and catalog which RR-intervals fail
            for i = 1:length(mat)
                if sum(isnan(mat(i,:))) > 0
                    % A NaN has been found
                    cut = [cut,i];
                end
            end

        end

        % Cut the rows from locations_matrix which fail
        locations_matrix(cut,:) = [];
    end
    
end