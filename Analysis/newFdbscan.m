function idx = newFdbscan(M1,labs,M2,epsilon,minpoints,plotting)
%Creates 3D plots of every combination of M1 features and returns id values
%of clusters
    %   Inputs:
    %       M1: [n-by-m matrix] a matrix where each column is a different feature and 
    %           each row is a different datapoint. Must have at least 2
    %           columns
    %       labs: [1-by-m cell array] of string labels for each column in M1
    %       M2: [n-by-1 vector] of 0's or 1's where 1 is problematic behavior
    %       epsilon: distance DBscan used to determine which points belong
    %       in the same cluster. The larger the data set, most likely the
    %       larger the epsilon. Set default to 20 and play around with it. 
    %       minpoints: Number of minimum points to make a cluster. The larger the data set, most likely the
    %       larger the number of minpoints. Set default to 20 and play around with it. 
    %       Plotting: boolean value, if true plots will generate

    %       
    %   Outputs:
    %       idx: id values for the each point, assigning each point to a cluster
    %       If plotting is set to 1 (true), function displays 3d plots of
    %       every combination of features from M1 with M2 as the z axis

    for z = 2:size(M1,2)
       
        idx(:,z) = dbscan(M1,epsilon,minpoints); 
                                               
        %PLOTTING  | | |
        %          V V V   

        % create n, value used for number of loop iterations 
        uni = unique(idx(:,z));  %returns array with no reptitions 
        numUni = numel(uni);  %counts number of elements 
        n=numUni-1;             %set n to total number of unique elements in idx (account for outliers by subtracting 1)      
    end     

    if plotting
        % colors
        colorbank = {};
        connectedColorBank = {};
        for i = 1:n
            colorbank{i} = {rand,rand,abs(rand-0.5)}; %abs(rand-0.5) is to avoid getting 1,1,1 which would be white
                                           
            connectedColorBank{i} = string("[" + colorbank{i}{1} + "," + colorbank{i}{2} + "," + colorbank{i}{3} + "]");
        end

        
        %loop generating clusters without outliers

        %(n choose k) = n!/(k!*(n-k)!)  = number of plots returned
        %where n = size(M1,2) and k = 2

        numPlots = factorial(size(M1,2))/(2*factorial(size(M1,2)-2)); 

        pos = 1; %pos of subplot, will increase by one every iteration of inner loop
     
        for j = 1:size(M1,2)-1
            for k = j+1:size(M1,2)
                subplot(2,ceil(numPlots/2),pos)
                pos = pos + 1;

                for i = 1:n
                    %plot3(arr1(idx1==-1,1),arr1(idx1==-1,2),arr1(idx1==-1,3),'.','color','b')  %all points not assigned a cluster 
                    plot3(M1(idx(:,j+1)==i,j),M1(idx(:,j+1)==i,k),M2(idx(:,j+1)==i,1),'.','color',connectedColorBank{i},'MarkerSize',12) % all points in diff clusters
                    hold on
                end

                xlabel(labs{j})
                ylabel(labs{k})
                zlabel("affect")
                title(labs{j} + " vs " + labs{k} + " with affect");
                legend;
                grid on
                grid minor
                hold off;


            end 
        end
    end
    idx(:,1) = [];  %deletes first column of idx    
end
