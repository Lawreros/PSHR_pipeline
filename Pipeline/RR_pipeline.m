function [outputArg1,outputArg2] = RR_pipeline(inputArg1,inputArg2)
%RR_pipeline
%   This function runs all of the analysis processes requested by the user
outputArg1 = inputArg1;
outputArg2 = inputArg2;

%eval allows for you to run a string like a code statement
eval('function_name(arguments)','O');

%Pass in dataframe and loop through all the steps of the
%preprocessing/analysis list. Create new entries in the structure based off
%of the steps taken.

%pass in structure and cell array of processing steps. This allows for each
%function to be made in column 1 with columns 2+ being different parameter
%values

%Only the preprocessing has order matter, so the rest of the analysis
%should be able to be modular


end
