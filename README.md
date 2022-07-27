# PSHR_pipeline
A small repository for analysis of the RR-interval and ECG data collected by the PSHR app found here:
https://github.com/Lawreros/PSHR_v2

The purpose of this pipeline is to provide a easy-to-use framework and function library for RR-interval and ECG data organization and analysis.
Future development for the repository will include a GUI for analysis and visualization of RR-interval and ECG data by non-coders.

## Repository Contents:
Summary of Repository contents:
```
/linear.m : a linear version of the MATLAB script for testing new functions/development.
    This file is recommended as the base for constructing an analysis pipeline.
/gui.m : the main file for the analysis GUI. Is currently in very eary development.
/coding_comp.m : a script to compare across the coding of a video by multiple people.
    This is used to calculate agreement between coders.

/sample
    A directory containing sample HR, ECG, and Affect coded data for testing purposes. All
    three files are from the same recording session.
    /A_coding.csv : sample coding file of the video coresponding to the recording session
        where the ECG and RR-interval data was collected.
    /ECG_A.txt : sample ECG data
    /HR_A.txt : sample HR data

/Analysis
    The directory where all analysis functions are contained. Any functions which result in
    metrics/measurements should be found here. See documentation inside of function for more
    detailed information.
    /pnnx_calc.m : function to calculate the pNN50 metric (other numbers aside from 50 can be specified)
    /rmssd_calc.m : calculates root mean standard deviation of successive differences between RR-intervals
    /sdnn_calc.m : calculates standard deviation of RR-intervals
    /sdsd_calc.m : calculates standard devaition of successive differences between RR-intervals

/Export
    The directory where all functions focused on either displaying or exporting data in different forms
    are located. See documentaiton inside of function for more detailed information.
    /poincare_plot.m : function which plots the Poincare plot for a given vector. Also exports the SD1
        and SD2 measurements for the plot.

/Import
    The directory where all functions related to loading data from files into usable formats for analysis
    are located. See documentation inside of function for more detailed information.
    /load_affect.m : 
    /pshr_load_data.m : 

/Pipeline
    The directory where all "approved" pipelines are contained. Each function here should be a self-contained
    analysis pipeline which performs a "complete" analysis of some kind.
    /RR_pipeline.m : [Not functional]

/Preprocess
    The direcotry where all functions related to the preprocessing of data before being analyzed are contained.
    See documentaiton inside of function for more detailed information.
    /acar.m : 
    /affect_mark.m : 
    /bandpass.m : 
    /ecg_PQRST.m : takes an ECG signal and locates the P, Q, R, S, and T-waves in the signal
    /ecg_rr_alignment.m : attempts to align an RR-interval vector with the RR-intervals estimated from the ECG
        data. This serves to align the RR and ECG data.
    /ecg_rr_conversion.m : takes an ECG signal and returns the location of the R-waves and the duration of the
        RR-intervals.
    /kamath.m : 
    /karlsson.m : 
    /malik.m : 

```

## Tutorial (subject to change as pipeline is developed)
1) Download the repository and open the `linear.m` file in your MATLAB editor
2) You should be able to run the program as long as you are running it from the directory `linear.m` is in
3) The resulting variable `Data` is a `struct` which contains the relevant RR-interval, ECG, and Affect information


## Data structure organization:
```
Data.HR
        .Raw : cell array containing the raw HR data collected from the provided file and vectorized
        .path : path to the directory containing the HR file(s) used
        .files : name(s) of the HR files used
        .Affect : cell array containing the start and stop indexes in .Raw for each of the affects observed

    .ECG
        .Raw : cell array containing the raw ECG data collected from the provided file and vectorized
        .path : path to the directory containing the ECG file(s) used
        .files : name(s) of the ECG files used
        .Affect : cell array containing the start and stop indexes in .Raw for each of the affects observed

    .Affect
        .Raw : cell array containing the raw ECG data collected from the provided file and vectorized
        .path : path to the directory containing the ECG file(s) used
        .files : name(s) of the ECG files used
        .Times : cell array containing the start and stop times (in seconds) for each affect noted in .Raw

```
