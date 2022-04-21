%Coding comparison
%This script is meant to compare the results between coding done by
%different people
clear all;

aff_1 = 'ignore_samples/2022-03-25_Montanez.csv';
aff_2 = 'ignore_samples/2022-03-25_Kessler.csv';


tab_1 = readtable(aff_1, 'Format', 'auto');
tab_2 = readtable(aff_2, 'Format', 'auto');

% Change column names so you don't have duplicate issues
tab_2.Properties.VariableNames{2} = 'Affect1_2';
tab_2.Properties.VariableNames{3} = 'Affect2_2';
tab_2.Properties.VariableNames{4} = 'Affect3_2';
tab_2.Properties.VariableNames{25} = 'on_camera_2';
tab_2.Properties.VariableNames{26} = 'problem_yn_2';

%Make table for comparison
comp_tab = tab_1(:,[1:4,25,26]);
comp_tab = [comp_tab, tab_2(:,[2:4,25,26])];

%free up some memory
clear tab_1 tab_2;


%Begin comparison iterations
[r,c] = size(comp_tab);

comp_index = {'Disagreement',aff_2;aff_1,0};
comp_agree = {'Affect'};

for i = 1:r
    comp_1 = table2cell(comp_tab(i,2:4));
    comp_1 = comp_1(~isempty(comp_1));
    comp_2 = table2cell(comp_tab(i,7:9));
    comp_2 = comp_2(~isempty(comp_2));
    
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

%Pad out and organize comp_index so that it has the same number of affects
%for both reviewers

unique([transpose(comp_index(3:end,1)), comp_index(1,3:end)]);


%Plot agreement and disagreement histogram


disp('done');