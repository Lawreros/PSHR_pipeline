function [mdl, pihats, AUC] = gen_regression(mat, target, regression_type, varargin)
% Function which taken in input matrix of features and returns the results
% of a variety of regressions. MATLAB's regression functions exclude any
% datapoints which contain NaN values in their regression.

% Required Inputs:
%   mat: [n-by-m matrix] the matrix you wish to use for regression. Each
%       row is assumed to be a different feature vector.
%
%   target: [n-by-1 matrix] vector containing the target classification/value
%       for each datapoint

% Optional Parameters:
%   ordinal: [bool] whether the data you are providing has ordinal target
%       values (i.e. it is categorical with some order involved with the
%       categories). Default is false, meaning 'nominal' will be used, denoting no
%       inherent ordering to categories.
%
%   test: [1-by-n vector] vector containing the ratio of each category
%       you want to have in the regression you are preforming. This will remove
%       any data points containing NaNs and create a subset of datapoints that
%       have the ratio you described. For example [0.2, 0.3, 0.5] will result
%       in the regression being performed on a subset of datapoints that are
%       20% category 0, 30% category 1, and 50% category 2. All datapoints excluded
%       for the sake of this ratio will be used for testing.
%       Default is [0.5, 0.5], if [] then all datapoints will be used.
%
%   verbose: [bool] whether to have the results of the various tests
%       printed to the Command Window.
%
%   regression_type: [string] the type of regression you wish to perform.
%       Current options are:
%           'linear': a linear regression
%           'log': a multinomial logistic regression

% Returns:
%   mdl: [p-by-4 matrix] matrix containing regression coefficients along
%       with the SE, tStats, and pValue for each coefficient

%   pihat: [q-by-n matrix] vector containing the predictive results of
%       running the regression model on the remaining datapoints not used for
%       training. Each column represents the probability that the datapoint
%       belongs to that category.


    % Parse any additional settings:
    p = inputParser;
    validSplit = @(x) isnumeric(x) && ismatrix(x);
    addParameter(p,'ordinal',false, @islogical);
    addParameter(p,'test',[0.5, 0.5], validSplit);
    addParameter(p,'verbose',true, @islogical);
    
    parse(p,varargin{:});
    
    [train, test, unused] = train_test_split(mat, target, [0.5 0.5], 'split', 0.8);
    
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
                [B,dev,stats] = mnrfit([train.cat_0;train.cat_1],categorical([zeros(length(train.cat_0(:,1)),1);ones(length(train.cat_1(:,1)),1)]), 'model', 'nominal');
                LL = stats.beta - 1.96.*stats.se;
                UL = stats.beta + 1.96.*stats.se;
                mdl = [stats.beta, stats.se, stats.t, stats.p, LL, UL];
                if p.Results.verbose
                    fprintf('    Beta\t SE \t tStats\t  pValue\n');
                    disp(mdl);
                    disp('Upper bounds');
                    disp(LL);
                    disp('Lower bounds');
                    disp(UL);
                end
                
                % Test the generated model:
                [pihat, dlow, dhi] = mnrval(B, [test.cat_0;test.cat_1],stats);
                cat_0_len = size(test.cat_0, 1);
                cat_1_len = size(test.cat_1, 1);
                pihats = [];
                
                
                for q = 1:3
                    if q == 1 % Worst case scenario
                        pdump(1:cat_0_len,1) = pihat(1:cat_0_len,1) - dlow(1:cat_0_len,1);
                        
                        pdump(cat_0_len+1:length(pihat),1) = pihat(cat_0_len+1:end,1) + dhi(cat_0_len+1:end,1);
                        
                    elseif q == 2 % Nothing
                        pdump = round(pihat(:,1));
                
                    else % Best case scenario
                        pdump(1:cat_0_len,1) = pihat(1:cat_0_len,1) + dhi(1:cat_0_len,1);
                        
                        pdump(cat_0_len+1:length(pihat),1) = pihat(cat_0_len+1:end,1) - dlow(cat_0_len+1:end,1);
                    end
                    
                    acc = [sum(round(pdump(1:cat_0_len,1)) == 1)/cat_0_len, sum(round(pdump(cat_0_len+1:end,1)) == 0)/cat_1_len];
                    
                    if p.Results.verbose
                        disp('Accuracy for non-target: '+string(acc(1)));
                        disp('Accuracy for target: '+string(acc(2)));
                    end
                    pdump = acc;
                
                    [X,Y,T,AUC] = perfcurve([zeros(cat_0_len,1);ones((cat_1_len),1)],pihat(:,2),1);
                    
                    pihats = [pihats; pdump];
                    
                end
                
            end
        otherwise
            disp('Unknown regression type, terminating function')
            return
    end
    
end