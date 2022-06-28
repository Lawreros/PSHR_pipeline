# PSHR_pipeline
A small repository for analysis of the RR-interval and ECG data collected by the PSHR app found here:
https://github.com/Lawreros/PSHR_v2

The purpose of this pipeline is to provide a easy to use GUI for analysis and visualization of RR-interval and ECG data. It also serves to output analysis data into formats that allow for easy analysis by other software.

## Repository Contents:
Relevant files for casual users:
- `sample` : folder containing sample HR, ECG, and Affect coded data
- `linear.m` : a linear version of the MATLAB script found in `main.m`, used for new function development and debugging
- `main.m` : the MATLAB script which calls the main GUI

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
