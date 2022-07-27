# PSHR_pipeline
A small repository for analysis of the RR-interval and ECG data collected by the PSHR app found here:
https://github.com/Lawreros/PSHR_v2

The purpose of this pipeline is to provide a easy-to-use framework and function library for RR-interval and ECG data organization and analysis.
Future development for the repository will include a GUI for analysis and visualization of RR-interval and ECG data by non-coders.

## Repository Contents:
Summary of Repository contents:
```
/linear.m : a linear version of the MATLAB script for testing new functions/development (debugging GUI's in MATLAB is a pain)
/gui.m : 
/coding_comp.m :

/sample
    directory containing sample HR, ECG, and Affect coded data
    /A_coding.csv
    /ECG_A.txt
    /HR_A.txt

/Analysis
    /pnnx_calc.m
    /rmssd_calc.m
    /sdnn_calc.m
    /sdsd_calc.m
/Export
    /poincare_plot.m
/Import
    /load_affect.m
    /pshr_load_data.m
/Pipeline
    /RR_pipeline.m
/Preprocess
    /acar.m
    /affect_mark.m
    /bandpass.m
    /ecg_PQRST.m
    /ecg_rr_alignment.m
    /ecg_rr_conversion.m
    /kamath.m
    /karlsson.m
    /malik.m

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
