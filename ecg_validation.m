% File for exploration of collected ECG data from the Student
clear all;

ecg_file = '~/Downloads/OneDrive_1_6-6-2022/AI2_20191007_023246.mat';
%ecg_file = '~/Downloads/OneDrive_1_6-6-2022/AI_20191007_022054.mat';

load(ecg_file);

peak = 800;
dis = 40;

[pks, locs] = findpeaks(transpose(AI2_20191007_023246mffECG), 'MinPeakProminence', peak, 'MinPeakDistance', dis);
plot(AI2_20191007_023246mffECG)
% [pks, locs] = findpeaks(transpose(AI_20191007_022054mffECG), 'MinPeakProminence', peak, 'MinPeakDistance', dis);
% plot(AI_20191007_022054mffECG)


% Generate RR-intervals form the index locations of the peaks
for i =2:length(locs)
    locs(i,2) = (locs(i,1)-locs(i-1,1))*(1000/PNSSamplingRate);
end

% Calculate RMSSD
summation = 0;
for i = 2:length(locs)
    summation = summation+(locs(i,2)-locs(i-1,2))^2;
end
RMSSD = sqrt(1/(length(locs))*summation);

% Calculate pNN50
PNNX = pnnx_calc_2(locs(:,2),50,false,false);



function [ret] = pnnx_calc_2(mat,diff,bin,band)
    % Calculates the percentage of adjacent NN-intervals that differ from
    % each other by more than "diff" milliseconds
    %   Inputs:
    %       mat: A [n-by-1] vector which contains the data you want to
    %       calculate pNNX for
    %       diff: [int] The minimum difference in milliseconds between
    %       successive NN-intervalse that you want to count
    %
    %       bin: [1-by-2 cell array] Used for creating a vector of the pNNX
    %       results from a sliding bin of Y seconds or entries. This takes the
    %       format of {index, 'units'}, so if you want to have a bin of the
    %       last 5 seconds: {5, 'second'} or if you want the last 5 measurements: {5, 'measure'}
    %       If you don't want this, set bin to false.
    %
    %       band: [2 int vector] The range [start, end] of values you want
    %       to calculate the pnnX of. If false, then analyze the whole
    %       range
    
    if band
        r_1 = band(1);
        r_2 = band(2);
    else
        [r_2,c] = size(mat);
        r_1 = 1;
    end
    
    % If they've decided to use the bin values
    if iscell(bin)
        a = bin{1}; % value
        b = bin{2}; % units
        
        ret = zeros(r_2-r_1,1);
        
        if strcmp(b,'second')
            for i = r_1:r_2
                count = 0;
                j = 0;
                
                % Check loop backward until you have reached the 'b'
                % seconds in the past through summing
                while (sum(mat(i-j:i,1))/1000) <= a
                    j = j+1;
                    if j == i
                        break;
                    end
                end
                
                if j > 1 && j < i % If there is more than one entry
                    for k = (i-j+2):i
                        if abs(mat(k,1) - mat(k-1,1))>= diff
                            count = count+1;
                        end
                    end
                    ret(i-r_1+1,1) = count/(j-1);
                else
                    ret(i-r_1+1,1) = 0;
                end
        
            end
        else % Looking at the past 'b' entries for the calculation
            for i = (a+1):(r_2-r_1+1)
                count = 0;
                for j = 0:(a-1)
                    if abs(mat(i-j,1) - mat(i-j-1,1))>= diff
                        count = count+1;
                    end
                end
                ret(i,1) = count/a;
            end
        end
        
    else
        % If they just want a percentage for a matrix
        count = 0;
        for i = (r_1+1):r_2
            if abs(mat(i,1)-mat(i-1,1))>= diff
                count = count+1;
            end
        end
        ret = count/(r_2-r_1);
    end
end