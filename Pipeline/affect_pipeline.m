function [] = affect_pipeline(aff_files)
% Pipeline which analyzes the prevalence and clustering of different
% affects

    aff_list = {'SIB','ISB','innappropriate face related behavior','polar strap adjustment/removal'...
        'repetitive behaviors','inappropriate movement','crying', 'pulling at pants'};

    %aff_list = false;
    
    Data = pshr_load('Affect',aff_files,'align',false,'verbose',false, 'meta', 10);

    dist_vec = [];
    for i =1:length(Data.Affect.Times)
        %Create a vector of zeros with affects on aff_list marked with
        %nonzero values
        aff_vec = affect_mark([], Data.Affect.Times{i},aff_list,false);
        
        %Calculate distances amongst and across affects
        markers = unique(aff_vec);
        
        for j =2:length(markers)
            inst = find(aff_vec==markers(j));
            inst = inst(2:end)-inst(1:end-1);
            dist_vec = [dist_vec; inst(inst>1)];
        end
        
    end
    histogram(dist_vec(dist_vec<180),'Binwidth', 2, 'Binlimits',[2,180]);
    title('Time between problematic behaviors [2, 180]');
    xlabel('Time (sec)');
    
    cap = 10;
    % Test for creating meta chunks:
    % Make a new affect called something like 'meta' and append it onto the
    % existing affect table. From there everything can run as normal...
%     inst = find(aff_vec==markers(j));
%     vec = inst(2:end)-inst(1:end-1);
%     idx = find(vec < cap & vec > 1);
%     starts =[];
%     ends =[];
%     for i=1:length(idx)
%         starts = [starts,inst(idx(i))];
%         ends = [ends,inst(idx(i)+1)];
%     end
%     Data.Affect.Times{1}{end+1,1} = 'meta_chunk';
%     Data.Affect.Times{1}{end,2} = starts;
%     Data.Affect.Times{1}{end,3} = ends;
    
end