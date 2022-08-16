function [mdl] = gen_regression(mat, target, regression_type, varargin)
% Function which taken in input matrix of features and returns the results
% of a variety of regressions.

% Inputs:
%   mat: [n-by-m matrix] the matrix you wish to use for regression. Each
%   row is assumed to be a different feature vector.
%
%   target: [n-by-1 matrix] vector containing the target classification/value
%   for each datapoint
%
%   ordinal: [bool] whether the data you are providing has ordinal target
%   values (i.e. it is categorical with some order involved with the
%   categories). Default is true.
%
%   regression_type: [string] the type of regression you wish to perform.
%   Current options are:
%       'linear': a linear regression
%       'log': a multinomial logistic regression


    % Parse any additional settings:
    p = inputParser;
    addRequired(p,'ordinal',true, @islogical);
    
    parse(p,varargin{:});


    switch regression_type
        case 'linear'
            mdl = fitlm(mat, target);
        case 'log'
            %https://www.mathworks.com/help/stats/mnrfit.html
            if p.Results.ordinal
                mdl = mnrfit(mat, target,'model','ordinal');
            else
                mdl = mnrfit(mat, target);
            end
        otherwise
            disp('Unknown regression type, terminating function')
            return
    end


end