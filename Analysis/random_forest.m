function [maccu, mprecision, mrecall, c_matrix, feature_importance] = random_forest(mdata, is_affect, tree_num, hold_out, feature_names)
    %creates a random forest model from the given data & number of trees
    %returns the percent accuracy of the model (maccu) & a feature
    %importance plot
    %   inputs:
    %       mdata: A [m-by-n] matrix, where n is the number of variables
    %       present in the data (n>=3).
    %       In addition, the columns must be in the following order:
    %           time  ...... 
    %       (i.e. time must be 1st col & everything else follows)
    %
    %       is_affect: A [n-by-1] vector of 0s and 1s.
    %       1 = problematic behavior
    %       0 = nonproblematic behavior
    %
    %       tree_num: an integer that represents the number of trees the
    %       model needs to contain when training. If false, then it is set
    %       to default (100 trees)
    %
    %       hold_out: a decimal between (0,1) that represent the percentage
    %       of data that's held out for testing
    %
    %       feature_names: a list containing names of all variables in order of columns.
    %       if list is nonempty: feature importance plot is needed, graphed
    %       in the order of variable names in the list given.
    %       if list is empty: feature importance plot is not needed.
    %
    %   outputs:
    %       maccu: model accuracy, in percentage
    % 
    %       mprecision: precision of the classification model
    %
    %       mrecall: recall of the classification model
    %
    %       c_matrix: confusion matrix in the form 
    %           [True Positives False Negatives, False Positives True Negatives]
    %
    %       OPTIONAL: a feature importance plot


    %determine num of trees
    if tree_num
        numTrees=tree_num;
    else
        numTrees=100;
    end

    %convert matrix to table
    hr_data=array2table(mdata);

    p_data=array2table(is_affect);
    hr_data=[hr_data,p_data];


    %fill in NaNs in time & 0 in bpm
    rows=size(hr_data,1);
    for row = 1:rows
        A=hr_data{row,:};
        TF=isnan(A);
        if TF(1)==true
            hr_data(row,1)=hr_data(row-1,1);
            hr_data(row,2)=hr_data(row-1,2);
        end
    end


    %training the model
    cv=cvpartition(size(hr_data,1),'HoldOut', hold_out);
    index=cv.test;
    dTrain=hr_data(~index,2:end);
    dTest=hr_data(index,2:end);
    testing=dTest(1:end,1:end-1);
    model=fitensemble(dTrain,'is_affect','Bag',numTrees,'Tree','Type','classification');
    prediction=predict(model,testing);

    %accuracy of prediction (this part I looked up)
    maccu=(sum(prediction==table2array(dTest(:,end)))/size(dTest,1))*100;

    %precision, recall, confusion matrix
    dTest_array = table2array(dTest);
    TP=0;
    TN=0;
    FP=0;
    FN=0;
    for row=1:size(dTest,1)
        if(prediction(row,end)==1 && dTest_array(row,end)==1)
            TP = TP+1;
        elseif(prediction(row,end)==1 && dTest_array(row,end)==0)
            FP = FP+1;
        elseif(prediction(row,end)==0 && dTest_array(row,end)==1)
            FN = FN+1;
        elseif(prediction(row,end)==0 && dTest_array(row,end)==0)
            TN = TN+1;
        end
    end

    mprecision = TP / (TP+FP);
    mrecall = TP / (TP + FN);
    c_matrix = [TP FN, FP TN];

    %prompt if feature importance plot is needed
    if ~isempty(feature_names)
        
        si = size(hr_data,2)-1;  %since hr_data has is_affect, need to exclude
        hr_data.Properties.VariableNames(1) = {'time'};
        hr_data.Properties.VariableNames(2:si) = feature_names;

        %feature importance
        feature_importance=oobPermutedPredictorImportance(model);
        figure
        bar(feature_importance)
        title('Feature Importance Estimates')
        xlabel('Predictor variable')
        ylabel('Importance')
        h=gca;
        h.XTickLabel=hr_data.Properties.VariableNames(2:end-1);
        h.XTickLabelRotation=45;
        h.TickLabelInterpreter='none';
    end

end
