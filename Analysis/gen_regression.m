function [mdl, pihats, t_stats, AUC] = gen_regression(mat, target, regression_type, varargin)
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
    
    
    % Get the correct ratio of categories: (for now just assume 50/50 for
    % problematic and nonproblematic)
    if p.Results.test
        
        [train, test, unused] = train_test_split(mat, target, [0.5 0.5], 'split', 0.8);
        
        % Run t-test on data
        [r,f] = size(train.cat_0);
        t_stats=[];
        for q = 1:f
            [h, p_, ci, stats] = ttest2([train.cat_0(:,q);test.cat_0(:,q);unused.cat_0(:,q)],[train.cat_1(:,q);test.cat_1(:,q);unused.cat_1(:,q)], 'Vartype', 'unequal');
            t_stats(1,q) = p_;
            disp(strcat('two-sample t-test result for feature ', string(q),' = ', string(p_)));
            disp(strcat('MEANS: NON-PROB : ', string(mean([train.cat_0(:,q);test.cat_0(:,q);unused.cat_0(:,q)])),' PROB : ', string(mean([train.cat_1(:,q);test.cat_1(:,q);unused.cat_1(:,q)]))));
            disp('95% estimated difference:');
            disp(ci);
            disp('Esitmated standard deviation');
            disp(stats.sd);
            %a = length([train.cat_0(:,q);test.cat_0(:,q);unused.cat_0(:,q)]);
            %b = length([train.cat_1(:,q);test.cat_1(:,q);unused.cat_1(:,q)]);
            %boxchart([repmat(0,a,1);repmat(1,b,1)],[train.cat_0(:,q);test.cat_0(:,q);unused.cat_0(:,q);train.cat_1(:,q);test.cat_1(:,q);unused.cat_1(:,q)]);
            
            [h,p_] = ttest([train.cat_0(:,q);test.cat_0(:,q)],[train.cat_1(:,q);test.cat_1(:,q)]);
            t_stats(1,f+q) = p_;
            disp(strcat('paired t-test result for feature: ',string(q),' =', string(p_)));
        end        
        
        %%%%
%         target(any(isnan(mat),2),:) = [];
%         mat(any(isnan(mat),2),:) = [];
%     
%         np_mat = mat(target==0,:);
%         p_mat = mat(target==1,:);
%         
%         % Run t-test on data
%         [r, f] = size(np_mat);
%         for q = 1:f
%             [h,p_] = ttest2(np_mat(:,q), p_mat(:,q));
%             disp(strcat('two-sample t-test result for feature: ',string(q),' =', string(p_)));
%         end
        
        
        % Calculate based off of smallest group (problematic)
%         samp = randsample(length(np_mat), round(length(p_mat)*1));
%         not_samp = [1:length(np_mat)];
%         not_samp(samp) = [];
%         
%         for q = 1:f
%             [h,p_] = ttest(np_mat(samp,q), p_mat(:,q));
%             disp(strcat('paired t-test result for feature: ',string(q),' =', string(p_)));
%         end
        
%         new_mat = [p_mat; np_mat(samp,:)];
%         % Create testing data
%         test_dat = np_mat(not_samp,:);
%         new_target = [ones(length(p_mat),1);zeros(length(samp),1)];
%         
%         if p.Results.verbose
%             disp(strcat('Subset of :',string(length(samp)),' from :',string(length(np_mat)),'for non-problematic'));
%             disp(strcat('Problematic :',string(length(samp))));
%         end
%         
%         mat = new_mat;
%         target = new_target;
%         clear new_mat new_target
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
%                 [B, dev, stats] = mnrfit(mat, categorical(target), 'model','nominal');
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
%                 pihat = mnrval(B, test_dat);
                [pihat, dlow, dhi] = mnrval(B, [test.cat_0;test.cat_1],stats);
                pihats = [];
                
                % Convert probabilities into binary values and account for
                % both lower and higher bounds
                for q = 1:3
                    if q == 1
                        pdump = round(pihat-dlow);
                        % if both probs are > 50% then set both to 0
                        pdump(sum(pdump,2)>1,:)=0;
                    elseif q ==2
                        pdump = round(pihat);
                    else
                        pdump = round(pihat+dhi);
                        % if both probs are > 50% then set both to 0
                        pdump(sum(pdump,2)>1,:)=0;
                    end

    %                 disp(strcat('Accuracy :', string(nansum(pihat(:,1))),'/',string(length(pihat)),...
    %                     '=',string(nansum(pihat(:,1))/length(pihat
                    acc = pdump - [2*ones(length(test.cat_0(:,1)),1);-2*ones(length(test.cat_1(:,1)),1)];

                    disp('Accuracy for non-prob: '+string(length(find(acc(:,1)==-1))/length(test.cat_0(:,1))));
                    disp('Accuracy for prob: '+string(length(find(acc(:,2)==3))/length(test.cat_1(:,1))));
                    pdump = [length(find(acc(:,1)==-1))/length(test.cat_0(:,1)), length(find(acc(:,2)==3))/length(test.cat_1(:,1))];
                
                    [X,Y,T,AUC] = perfcurve([zeros(length(test.cat_0(:,1)),1);ones(length(test.cat_1(:,1)),1)],pihat(:,2),1);
                    
                    if AUC > 0.6
                        disp("hey!");
                    end
                    
                    pihats = [pihats; pdump];
                end
                
            end
        otherwise
            disp('Unknown regression type, terminating function')
            return
    end
    

end