function ecg_rr_alignment(rr, ecg)
% Function which attempts to align the ecg data with the RR-interval data
% Inputs:
%   rr: [matrix]
%   ecg: [matrix]


% Get the rr-interval estimates from the ECG data
peak = 800;
dist = 40;
freq = 130;
locs = ecg_rr_conversion(ecg, peak, dist, freq);



% Convert RR and ECG vectors into very long strings that are comma
% seperated. Also divide by 10 in order to ignore very minor differences
rr_string = "";
for i = 1:length(rr)
    rr_string = strcat(rr_string, string(rr(i)));
end


ecg_string = "";
for i=1:length(locs)
    ecg_string = strcat(ecg_string,string(locs(i)));
end


%editDistance calculates the maximum number of insertions/deletions in
%order to convert the first string into the second string
% TODO: Look into the cost function for editDistance, want to make
% deletions more common than insertions
d = editDistance(ecg,rr);


end

function []=ignore_this()
    %Placeholder function to store previous work
    format long g;

    % PhaseSpace file path
    PS_path = {'/data/PHYSMON/EGI/'};
    % List of PhaseSpace files
    PS_File_list = {'Ross_190707.txt'};%'Ross_190602_1.txt', 'Ross_190602_3.txt', 'Ross_190602_4.txt'};


    % EGI file path
    EGI_path = {'/data/PHYSMON/EGI/'};
    % List of EGI files
    EGI_File_list = {'Ross_190707_20190707_070500.mat'};%'Ross_190602_1_20190602_111912.mat','Ross_190602_3_20190602_121609.mat','Ross_190602_4_20190602_123756.mat'};
    % Start times
    EGI_Start = [25499337.183];%40750890.381,44167873.944,45475559.147];


    %% Load PolarStrap Data
    plot_prep = [0,0,0];

    for qq=1:length(PS_File_list)
        clear RR_raw;
        File_name = strcat(PS_path{1},PS_File_list{qq});

        data(qq).file_names = File_name;

        % Read information from MRI_session
        q = 1;
        RR_raw = {};
        for z = 1:length(File_name)
            fid = fopen(File_name);

            line = fgetl(fid);

            % Read in information
            i = 1;
            while line ~= -1
                RR_raw(i,:) = textscan(line, '%f:%f:%f %s %d %d %d', 'delimiter', '\t');
                line = fgetl(fid);
                i=i+1;
            end

        end
        % Eliminate "Heart Rate"
        RR_raw(:,8) = RR_raw(:,5);
        RR_raw(:,4:5) = [];

        % Normalize times
        for i = 1:length(RR_raw)
                RR_raw{i,3} = RR_raw{i,1}*3600+RR_raw{i,2}*60+RR_raw{i,3};
        end
        RR_raw(:,1:2) = [];

        for i = 1:length(RR_raw)
            RR_raw{i,1} = RR_raw{i,1}*1000;

            if isempty(RR_raw{i,3})
                RR_raw{i,3} = 0;
            end
            if isempty(RR_raw{i,2})
                RR_raw{i,2} = 0;
            end

        end


        %% Process data


        % Organize for Plotting
        for i = 2:(length(RR_raw)-1)

            if RR_raw{i,3} ~= 0
                plot_prep(end+1,1) = RR_raw{i,1};
                plot_prep(end,2) = RR_raw{i,2};
                if isempty(RR_raw{i,4})
                else
                plot_prep(end,3) = RR_raw{i,4};
                end
                plot_prep(end+1,1) = NaN; %((RR_raw{i+1,1}-RR_raw{i,1})/2)+RR_raw{i,1};
                plot_prep(end,2) = RR_raw{i,3};
            else
                plot_prep(end+1,1) = RR_raw{i,1};
                plot_prep(end,2) = RR_raw{i,2};
                if isempty(RR_raw{i,4})
                else
                plot_prep(end,3) = RR_raw{i,4};
                end
            end

        end

        % Remove blanks and over 1.6s RR-interval
        for i=1:length(plot_prep)
            if plot_prep(i,2) == 0
                plot_prep(i,2) = NaN;
            elseif plot_prep(i,2) > 1600 || plot_prep(i,2) < 100
                plot_prep(i,2) = NaN;
            end
        end

        data(qq).PS = plot_prep(2:end,:);
        plot_prep = [0,0];

    end

    %% Load in EGI data
    for qq = 1:length(EGI_File_list)
        temp = load(strcat(EGI_path{1},EGI_File_list{qq}), '*ECG');
        a = temp.(strcat(EGI_File_list{qq}(1:end-4),'mffECG'));

        plot(a);

        clear temp;
        [pks_a,locs_a]=findpeaks(a,'MinPeakProminence',1500);
        locs_a = locs_a';

        locs_a = EGI_Start(qq) + locs_a*4;   % Convert time from 250 Hz reading to 1000 Hz reading
        for i = 2:length(locs_a)
            locs_a(i,2) = locs_a(i,1)-locs_a(i-1,1);
        end

        data(qq).EGI = locs_a(2:end,:);
        data(qq).align(:,3:4) = locs_a(2:end,:);

    end

    clear locs_a pks_a a;

    %% Levenshtein Alignment (adapt PS to fit EGI?)

    for qq = 1:length(EGI_File_list)

        lev = zeros(length(data(qq).EGI),length(data(qq).PS));
        lev(1,1) = 0;

        for i = 2:length(data(qq).EGI)
            lev(i,1) = i-1;
        end
        for j = 2:length(data(qq).PS)
            lev(1,j) = j-1;
        end


        for j = 2:length(data(qq).PS)
            for i = 2:length(data(qq).EGI)
                if data(qq).EGI(i,2) == data(qq).PS(j,2)
                    substitutionCost = 0;
                else
                    substitutionCost = abs(data(qq).EGI(i,2) - data(qq).PS(j,2))/10;
                end
                lev(i,j) = min([lev(i-1,j)+1, lev(i,j-1)+1, lev(i-1, j-1)+substitutionCost]);
            end
        end

        [Y,I] = min(lev(end,:));

        [r,c] = size(lev);
        c = I;

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

        % Flip move
        move = flip(move,2);

        clear i j c r I line pks_a q q_c q_r Y substitutionCost;

        % Realign
        j = 1;
        plot_prep_2=[];

        for i = 1:length(move)

            if move(1,i) == 0
                plot_prep_2(end+1,1:2) = data(qq).PS(j,1:2);
                j = j+1;
            elseif move(1,i) == 1
                j = j+1;
            elseif move(1,i) == 2
                plot_prep_2(end+1,1:2) = NaN;
            end

        end
        data(qq).align(:,1:2) = plot_prep_2;


        for zz = 1:length(data(qq).align)
            data(qq).align(zz,5) = data(qq).align(zz,1) - data(qq).align(zz,3);

        end


    end

    clear plot_prep_2 i j z zz qq plot_prep move lev;


    %% Compare

    % plot(locs_a(2,:), 'LineWidth',2);
    % hold on;
    % plot(plot_prep_2);
    % 
    % 
    % plot_prep_2(2,1:length(locs_a)) = locs_a(2,:);
    % 
    % for j = 1:length(plot_prep_2)
    %     plot_prep_2(3,j) = abs(plot_prep_2(1,j) - plot_prep_2(2,j));
    %     
    % end
    % disp(nanmean(plot_prep_2(3,2:end)));
    % disp(nanstd(plot_prep_2(3,2:end)));
    % 
    % 

end
