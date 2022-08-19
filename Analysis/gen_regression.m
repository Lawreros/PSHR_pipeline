function [mdl, pihat] = gen_regression(mat, target, regression_type, varargin)
% Function which taken in input matrix of features and returns the results
% of a variety of regressions. MATLAB's regression functions exclude any
% datapoints which contain NaN values in their regression.

% Inputs:
%   mat: [n-by-m matrix] the matrix you wish to use for regression. Each
%   row is assumed to be a different feature vector.
%
%   target: [n-by-1 matrix] vector containing the target classification/value
%   for each datapoint
%
%   ordinal: [bool] whether the data you are providing has ordinal target
%   values (i.e. it is categorical with some order involved with the
%   categories). Default is false, meaning 'nominal' will be used, denoting no
%   inherent ordering to categories.
%
%   test: [1-by-n vector] vector containing the ratio of each category
%   you want to have in the regression you are preforming. This will remove
%   any data points containing NaNs and create a subset of datapoints that
%   have the ratio you described. For example [0.2, 0.3, 0.5] will result
%   in the regression being performed on a subset of datapoints that are
%   20% category 0, 30% category 1, and 50% category 2. All datapoints excluded
%   for the sake of this ratio will be used for testing.
%   Default is [0.5, 0.5], if [] then all datapoints will be used.
%
%   verbose: [bool] whether to have the results of the various tests
%   printed to the Command Window.
%
%   regression_type: [string] the type of regression you wish to perform.
%   Current options are:
%       'linear': a linear regression
%       'log': a multinomial logistic regression


    % Parse any additional settings:
    p = inputParser;
    validSplit = @(x) isnumeric(x) && ismatrix(x);
    addParameter(p,'ordinal',false, @islogical);
    addParameter(p,'test',[0.5, 0.5], validSplit);
    addParameter(p,'verbose',true, @islogical);
    
    parse(p,varargin{:});
    
    
    % Get the correct ratio of categories: (for now just assume 50/50 for
    % problematic and nonproblematic)
    if p.Results.test
        target(any(isnan(mat),2),:) = [];
        mat(any(isnan(mat),2),:) = [];
    
        np_mat = mat(target==0,:);
        p_mat = mat(target==1,:);
        
        % Calculate based off of smallest group (problematic)
        samp = randsample(length(np_mat), round(length(p_mat)*1));
        not_samp = [1:length(np_mat)];
        not_samp(samp) = [];
        
        new_mat = [p_mat; np_mat(samp,:)];
        % Create testing data
        test_dat = np_mat(not_samp,:);
        new_target = [ones(length(p_mat),1);zeros(length(samp),1)];
        
        if p.Results.verbose
            disp(strcat('Subset of :',string(length(samp)),' from :',string(length(np_mat)),'for non-problematic'));
            disp(strcat('Problematic :',string(length(samp))));
        end
        
        mat = new_mat;
        target = new_target;
        clear new_mat new_target
    end
    
    switch regression_type
        case 'linear'
            mdl = fitlm(mat, target);
            if p.Results.verbose
                disp(mdl);
            end
        case 'log'
            %https://www.mathworks.com/help/stats/mnrfit.html
            if p.Results.ordinal
                [B, dev, stats] = mnrfit(mat, categorical(target),'model','ordinal');
                mdl = [stats.beta, stats.se, stats.t, stats.p];
                if p.Results.verbose
                    fprintf('    Beta\t SE \t tStats\t  pValue\n');
                    disp(mdl);
                end
            else
                [B, dev, stats] = mnrfit(mat, categorical(target), 'model','nominal');
                mdl = [stats.beta, stats.se, stats.t, stats.p];
                if p.Results.verbose
                    fprintf('    Beta\t SE \t tStats\t  pValue\n');
                    disp(mdl);
                end
                
                % Test the generated model:
                pihat = mnrval(B, test_dat);
                
                % Convert probabilities into binary values
                pihat = round(pihat);
                
                disp(strcat('Accuracy :', string(nansum(pihat(:,1))),'/',string(length(pihat)),...
                    '=',string(nansum(pihat(:,1))/length(pihat))));
                
            end
        otherwise
            disp('Unknown regression type, terminating function')
            return
    end
    

end