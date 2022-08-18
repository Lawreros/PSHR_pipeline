# How to Contribute to the PSHR Pipeline
This document will outline the philosophy of the PSHR_pipeline and the general guidelines behind its construction
in order to minimize the amount of friction involved in implementing code from collaborators.

## The `Data` Structure
The PSHR_pipeline is built around the `Data` structure, where all loaded data and where all loaded data and analysis results are stored. 
The thinking behind this is to have a consolidated `DataFrame` (as seen in the R programing language) which different analysis functions
can operate on and store their results in.

- Each loaded `HR`, `ECG`, or `Affect` file will be stored in either `Data.HR.Raw`, `Data.ECG.Raw`, `Data.Affect.Raw` 

### Data Guidelines:
Here are some general rules that should be applied to 

- If timestamps are included in a given matrix, they should be prioritized to the earliest column possible


### `Data` example:
An example of this would be loading in the Data 

- see README.md for more information

```
% In the example below, a function is called which averages the RR-interval values from the first entry of Data.HR.Raw.
% The result is then stored in first entry of cell array Data.HR.Average.

Data = pshr_load('HR', {'./sample/HR_A.txt', './sample/HR_B.txt', NaN}, ...
		'ECG', {'./sample/ECG_A.txt','./sample/ECG_B.txt', './sample/ECG_C.txt'} ,
		'Affect', {'./sample/A_coding.csv', './sample/B_coding.csv', NaN});


Data.HR.Average{1} = calc_average(Data.HR.Raw{1}(:,3));

```

## Function assumptions

- Should always assume that any analysis of timeseries data will be on data in the form of a matrix or vector, where each 
- Do not make any changes to `Data.HR.Raw` or `Data.ECG.Raw`, as they contain the unaltered data from the text files.
- (?) Preprocessing should be put into `Data.HR.PP` or `Data.ECG.PP` as the "most recent/best data to use"
- The input arguments for any given function should be the minimum required to preform the analysis and should ideally not rely on any past analysis/functions. Asking for `Data` structure itself should be done a infrequently as possible
- Example of function (see below)
- Variable names should be descriptive without being too verbose in order to improve readiblity
- Input arguments should be organized from left to right as `(dataset_1, dataset_2, ..., parameter_1, parameter_2,...)`, where the data/variable being analyzed/acted upon is given priority in the listing.
- A docstring explaining the inputs and outputs for each function is REQUIRED, with additional commenting of the code so that others can understand what is occurring
- The acceptable data types for each input, as well as a short description for what each parameter influences

### Modular function example:

```
sample_matrix_1 = [1, 0;
				 2, 1;
				 6, 1;
				 5, 1;
				 2, 0;
				 6, 0;
				 8, 1]

sample_matrix_2 = [7, 1, 0;
				 5, 2, 1;
				 8, 6, 1;
				 2, 5, 1;
				 6, 2, 0;
				 9, 6, 0;
				 10, 8, 1]


% Note that the function only asks for one vector, this is to make it robust
% against any addition/removal of columns unrelated to the column in question.
% For example, this funcition can be called on both sample_matrix_1(:,1) and
% sample_matrix_2(:,2), returning the same result.

% Where possible you want to make functions modular, i.e. self-contained and 
% able to be called easily in any order


output = bandpass(sample_matrix_1(:,1), 2, 5, false);
output_2 = bandpass(sample_matrix_2(:,2), 2, 5, false);


function [ret] = bandpass(mat, l_band, u_band, rang)
    % Applies a bandpass filtering to vector provided. Any vector
    % value outside of the range specified by the lower and upper bounds is
    % replaced with a NaN
    %   Inputs:
    %       mat: [n-by-1 vector] vector of values which are being bandpass
    %       filtered
    %       source: [string], Which matrix from the structure you want to
    %       use
    %       l_band: [int], the lower bounding value for the bandpass filter
    %       u_band: [int], the upper bounding value for the bandpass filter
    %       rang: [2 int vector] The range [start, end] of values you want
    %       to use the bandpass on. If false, then analyze the whole
    %       range
    
    %   Returns:
    %       ret: [m-by-1 vector] vector of values that have been bandpass
    %       filtered. This should be the length defined by rang

    if rang
        r_1 = rang(1);
        r_2 = rang(2);
    else
        [r_2,c] = size(mat);
        r_1 = 1;
    end

    %Create copy of matrix to edit
    ret = mat(r_1:r_2,1);
    
    for i = 1:length(ret)
        if (ret(i,1)>= u_band) || (ret(i,1) <= l_band)
            ret(i,1) = NaN;
        end
    end 
end


end
```
