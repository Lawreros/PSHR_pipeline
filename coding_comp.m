%Coding comparison
%This script is meant to compare the results between coding done by
%different people
clear all;

aff_1 = "ignore_samples/2022-03-25_Montanez.csv";
aff_2 = "ignore_samples/2022-03-25_Kessler.csv";


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

comp_results = {'Disagreement','nothing';'nothing',0};

for i = 1:r
    comp_1 = table2cell(comp_tab(i,2:4));
    comp_2 = table2cell(comp_tab(i,7:9));
    
    qq = comp_1(~ismember(comp_1, comp_2));
    zz = comp_2(~ismember(comp_2, comp_1));
    
    if ~isempty(qq) || ~isempty(zz) %If there is a disagreement
        disp('Disagreement found');
        
        if isempty(qq)
            b = strjoin(zz,'&');
            a = 'nothing';
        elseif isempty(zz)
            b = 'nothing';
            a = strjoin(qq, '&');
        else % Both groups have a unique affect
            b = strjoin(zz,'&');
            a = strjoin(qq,'&');
        end
        
        %key = strcat(a, '_vs_', b);
        
        
%         % Add to the array
%         if ismember(key, comp_results(:,1))
%             %find where the row for this category of disagreement is
%             j = find(ismember(comp_results(:,1),key));
%             %append time of disagreement
%             comp_results{j,2} = [comp_results{j,2},comp_tab{i,1}];
%         else
%             comp_results{end+1,1} = key;
%             comp_results{end,2} = comp_tab{i,1};
%         end
        
        %Add to the array
        
        if ismember(a,comp_results(:,1)) && ismember(b,comp_results(1,:))
            j = find(ismember(comp_results(:,1),a));
            k = find(ismember(comp_results(1,:),b));
            comp_results{j,k} = [comp_results{j,k},comp_tab{i,1}];
            
        elseif ismember(a,comp_results(:,1)) && ~ismember(b,comp_results(1,:)) %only a already exists
            j = find(ismember(comp_results(:,1),a));
            comp_results{1,end+1} = b; %add the b category
            comp_results{j,end} = comp_tab{i,1};
        
        elseif ismember(b,comp_results(1,:)) %only b already exists
            k = find(ismember(comp_results(1,:),b));
            comp_results{end+1,1} = a;
            comp_results{end,k} = comp_tab{i,1};
            
        else %neither exist
            comp_results{end+1,1} = a;
            comp_results{1,end+1} = b;
            comp_results{end,end} = comp_tab{i,1};
            
        end
        


    end
    
end

disp('done');