% Coding_comp
%   This script is meant to compare the results between coding done by
%   different individuals, saving the results to an excel file

clear all;

% Location of directory containing the files to be compared
inp_dir = 'Reliability_Training/';

aff_1 = 'Chat&Chew_2022-03-18_1255_V01_ANSWER_KEY.csv';
aff_2 = {'Chat&Chew_2022-03-18_1255_V01_Liu.csv',...
    'Chat&Chew_2022-03-18_1255_V01_Baldie.csv',...
    'Chat&Chew_2022_03-18_1255_V01_Zhao.csv'};

% % aff_1 the "ground truth"
% aff_1 = 'SalivaPlay_Chat & Chew_2022-06-03_1255_V01_Kessler.csv';

% aff_2 the file(s) to compared to aff_1
% aff_2 = {'SalivaPlay_Chat & Chew_2022-06-03_1255_V01_Montanez.csv',...
%         'SalivaPlay_Chat & Chew_2022-06-03_1255_V01_Kessler.csv'};

% Where to save the results of the coding comparisons. If there is a 
out_file = 'test_comp.xlsx';


% What column to use when aligning the files (generally the Time_sec column
% for timepoint alignment)

align_col = 1;


% Which columns to compare across, including the timestamp column (which is
% assumed to be the first entry in comp_cols in both files. The columns to compare across
% are generally the 'Affect1', 'Affect2', 'Affect3' columns).

comp_cols = [2,3,4];


for z = 1:length(aff_2)


    tab_1 = readtable(strcat(inp_dir,aff_1), 'Format', 'auto');
    tab_2 = readtable(strcat(inp_dir,aff_2{z}), 'Format', 'auto');

    % Change column names so you don't have duplicate issues
%     tab_2.Properties.VariableNames{2} = 'Affect1_2';
%     tab_2.Properties.VariableNames{3} = 'Affect2_2';
%     tab_2.Properties.VariableNames{4} = 'Affect3_2';
%     tab_2.Properties.VariableNames{17} = 'on_camera_2'; %TODO: This is for the new coding format, will be wrong for 03-18 coding
%     tab_2.Properties.VariableNames{18} = 'problem_yn_2'; % ""

    
    comp_tab = tab_1(:,[align_col,comp_cols]);
    comp_tab = outerjoin(comp_tab, tab_2(:,[align_col,comp_cols]),'Key',tab_1.Properties.VariableNames{align_col});
    comp_tab = convertvars(comp_tab, @isnumeric, @nanblank);
    
    
    % Standardize everything:
    for i = 1:length(comp_tab.Properties.VariableNames)
        %t.Var1(strcmp(t.Var1,'N/A')) = {''};
        var_nam = comp_tab.Properties.VariableNames{i};
%         sum(strcmp(comp_tab.(var_nam),' '))
        if sum(strcmp(comp_tab.(var_nam),' ')) > 0
            comp_tab.(var_nam)(strcmp(comp_tab.(var_nam),' ')) = {''};
        end
    end
    
    %Make table for comparison
%     comp_tab = tab_1(:,[1:4,17,18]);%25,26]);
%     comp_tab = [comp_tab, tab_2(:,[2:4,17,18])];%25,26])];

    %free up some memory
    clear tab_1 tab_2;


    %Begin comparison iterations
    [r,c] = size(comp_tab);

    comp_index = {'Disagreement',aff_2{z};aff_1,0};
    comp_agree = {'Affect'};

    for i = 1:r
%         comp_1 = table2cell(comp_tab(i,2:4));
        comp_1 = table2cell(comp_tab(i,2:length(comp_cols)+1));
        %comp_1 = comp_1(~isempty(comp_1));
%         comp_2 = table2cell(comp_tab(i,7:9));
        comp_2 = table2cell(comp_tab(i,length(comp_cols)+3:end));
        %comp_2 = comp_2(~isempty(comp_2));

        qq = comp_1(~ismember(comp_1, comp_2));
        zz = comp_2(~ismember(comp_2, comp_1));

        %Create list of agreements
        pp = unique([comp_1(ismember(comp_1, comp_2)),comp_2(ismember(comp_2, comp_1))]);

        for q = 1:length(pp)
            if ismember(pp{q},comp_agree(:,1))
                j = find(ismember(comp_agree(:,1),pp{q}));
                comp_agree{j,2} = comp_agree{j,2}+1;
            else
                comp_agree{end+1,1} = pp{q};
                comp_agree{end,2} = 1;
            end

        end

        if ~isempty(qq) || ~isempty(zz) %If there is a disagreement
            disp(strcat('Disagreement found at time:', string(comp_tab{i,1})));

            if isempty(qq)
                b = strjoin(zz,'&');
                a = 'n/a';
            elseif isempty(zz)
                b = 'n/a';
                a = strjoin(qq, '&');
            else % Both groups have a unique affect
                b = strjoin(zz,'&');
                a = strjoin(qq,'&');
            end

            %key = strcat(a, '_vs_', b);

            %Add to the array

            if ismember(a,comp_index(:,1)) && ismember(b,comp_index(1,:))
                j = find(ismember(comp_index(:,1),a));
                k = find(ismember(comp_index(1,:),b));
                comp_index{j,k} = [comp_index{j,k},comp_tab{i,1}];

            elseif ismember(a,comp_index(:,1)) && ~ismember(b,comp_index(1,:)) %only a already exists
                j = find(ismember(comp_index(:,1),a));
                comp_index{1,end+1} = b; %add the b category
                comp_index{j,end} = comp_tab{i,1};

            elseif ismember(b,comp_index(1,:)) %only b already exists
                k = find(ismember(comp_index(1,:),b));
                comp_index{end+1,1} = a;
                comp_index{end,k} = comp_tab{i,1};

            else %neither exist
                comp_index{end+1,1} = a;
                comp_index{1,end+1} = b;
                comp_index{end,end} = comp_tab{i,1};

            end

        end

    end

    % Convert the comp_index to cells with strings so they can be saved as
    % excell files

    bb = comp_index;
    cc = comp_index;
    [r,c] = size(comp_index);
    for i = 3:r
        for j = 3:c
            bb{i,j} = strjoin(string(bb{i,j}),', ');
            cc{i,j} = length(comp_index{i,j});
        end
    end
    

    % Add the seconds of disagreement to comp_agree for easy exporting
    aff_list = unique([transpose(comp_index(3:end,1)), comp_index(1,3:end)]);

    for i = 1:length(aff_list)
        c = find(ismember(comp_index(1,:),aff_list{i}));
        r = find(ismember(comp_index(:,1),aff_list{i}));

        count = 0;

        if r
            count = count + sum(cellfun('length',comp_index(r,3:end)));
        end

        if c
            count = count + sum(cellfun('length',comp_index(3:end, c)));
        end

        if isempty(r)==0 && isempty(c)==0
            count = count - cellfun('length',comp_index(r,c));
        end

        % Append to comp_agree or make new row
        r = find(ismember(comp_agree(:,1), aff_list{i}));
        if r
            comp_agree{r,3} = count;
        else
            comp_agree{end+1,1} = aff_list{i};
            comp_agree{end,2} = 0;
            comp_agree{end,3} = count;
        end
    end

    % Calculate agreement rate for total problematic behavior
%     comp_agree{end+1,1} = "Existence of Problem Behavior";
%     comp_agree{end,3} = sum(abs(comp_tab.problem_yn - comp_tab.problem_yn_2));
%     comp_agree{end,2} = height(comp_tab) - comp_agree{end,3};

    % Calculate rates of agreement
    comp_agree{1,4} = "Agreement Rate";
    comp_agree{1,2} = "Seconds of Agreement";
    comp_agree{1,3} = "Seconds of Disagreement";
    for i = 2:length(comp_agree)
        
        if isempty(comp_agree{i,2})
            comp_agree{i,2} = 0;
        end
        
        if isempty(comp_agree{i,3})
            comp_agree{i,3} = 0;
        end
        
        comp_agree{i,4} = comp_agree{i,2} / (comp_agree{i,2}+comp_agree{i,3});
    end
    
    
    % Save information into excell file
    if ischar(out_file)
        writecell(cc, out_file, 'Sheet', z, 'Range', 'A2:Z15')
        writecell(bb, out_file, 'Sheet', z, 'Range', 'A16:Z29')
        writecell(comp_agree, out_file, 'Sheet',z,'Range','A30:Z50');
    end
        
    clear comp_agree comp_tab comp_tab bb;   
end
    
%Plot agreement and disagreement histogram







function output = nanblank(values)
    mask = isnan(values);
    if nnz(mask)
      output = string(values);
      output(mask) = '';
      output = char(output);
    else
        output = values;
    end
end