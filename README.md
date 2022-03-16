# PSHR_pipeline
Repository for analysis of the RR-interval and ECG data collected by the PSHR app found here:
https://github.com/Lawreros/PSHR_v2

The purpose of this pipeline is to provide a GUI for analysis and visualization of RR-interval and ECG data. It also serves to output analysis data into formats that allow for easy analysis by other software.

## Repository Contents:
Relevant files for casual users:
- `sample` : folder containing sample HR and ECG data
- `main.m` : the MATLAB script which calls the main GUI. This is the function that you want to run in you MATLAB editor.


## Tutorial (subject to change as pipeline is developed)
1) Download the repository and open the `main.m` file in your MATLAB editor
2) Run the program, a GUI should appear
3) In the menu at the top of the window, click on `Import`. You will be given two options as to what type of data you would like to load. Selecting one will make another window appear asking you to select the files. You can select more than one file at a time and the click `Open` to load the files. They should now appear on the GUI. You can use the `Import` button to load files of a different type, or a new set of files.
4) Highlighting a file (or multiple files) that you have loaded in the GUI and pressing the `View` button will create a figure that displays a basic line plot of the data. This allows for catching any errors that may appear.

