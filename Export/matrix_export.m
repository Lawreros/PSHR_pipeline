function [] = matrix_export(mat,fil_name, varargin)
%Simple function which exports the contents of a given matrix to a .csv
%file

% Required Inputs:
%   mat: [n-by-m matrix] The matrix you wish to save to a .csv file

%   fil_name: [string] What you want the created csv file to be called. This will be
%   saved in your current directory unless the full file path is specified 
%   in file name.

% Optional Parameters:
%   fil_type: [string] What file type you want to save the matrix as.
%   Either 'csv' for a csv file,'txt' for a tab delimited text file, or
%   'mat' for a matlab file of the matrix. Default is 'csv'.
%
%
%   var_names: [m-by-1 cell array] A cell array containing the title you
%   wish to give to each column in the csv file (requires fil_type = 'csv').
%   If false, then generic titles are generated. Default is false.


    p = inputParser;
    addParameter(p, 'var_names',false, @iscell);
    addParameter(p, 'fil_type', 'csv', @ischar);
    parse(p,varargin{:});

    if ~isempty(p.Results.var_names)
        mat = array2table(mat,'VariableNames',p.Results.var_names);
        writetable(mat,strcat(fil_name,'.csv'),'Delimiter',',');   
    elseif p.Results.fil_type == 'txt'
        writematrix(mat,strcat(fil_name,'.txt'),'Delimiter','tab');
    elseif p.Results.fil_type == 'mat'
        save(strcat(fil_name,'.mat', mat));
    else
        writematrix(mat,strcat(fil_name,'.',p.Results.fil_type));

    end
end
