# About This Folder
This Repository contains MATLAB and Python code files used in the work accompanying the thesis entitled **"Metabolic Engineering of Bacteria to Convert CO<sub>2</sub> to Products of Value"**. The structure and use of the files found here will be explaine briefly in this document.

## A General Remark
Many of the coding files here contain references to file locations and file names. These are purely left in for demonstration purposes, if these files were to be used by the Reader, they would need to point the scripts to wherever they might have equivalent files in their working environment.

## MATLAB_Data_Analysis
This folder contains general files for data analysis using MATLAB. All files are at least compatible with version 2020b or newer, older versions have not been checked. The purpose of each file will be briefly pointed out here.
  1. [ALE_growth_rates_meta_analysis_BioLector.m](MATLAB_Data_Analysis/ALE_growth_rates_meta_analysis_BioLector.m):
       - Generates a scatter plot displaying maximum growth rates for a number of analysed growth curves.
       - Requires analysed growth curves and estimated growth rates for each of the displayed samples.
  2. [BB_Batch_Analysis_w_online_data.m](MATLAB_Data_Analysis/BB_Batch_Analysis_w_online_data.m):
       - Generates plots using previously analysed offline biomass data and online-measured data from the ROSITA software from Bionet to display batch fermentation profiles (run in F0 Baby Bionet units).
       - Accepts data for multiple fermentations and plots them individually.
       - Required the respective data files.
  3. [BB_Conti_Analysis.m](MATLAB_Data_Analysis/BB_Conti_Analysis.m)
       - Generates plots using previously analysed offline biomass data and online-measured data from the ROSITA software from Bionet to display continuous fermentation profiles (run in F0 Baby Bionet units).
       - Required the respective data files.
  4. [BB_Conti_Analysis_Compare_multiple_runs.m](MATLAB_Data_Analysis/BB_Conti_Analysis_Compare_multiple_runs.m):
       - Compares multiple continuous fermentation profiles from F0 Baby Bionet units.
       - Required the respective data files.
  5. [BioLector_Analysis.m](MATLAB_Data_Analysis/BioLector_Analysis.m):
       - Analyses growth curves obtained with a BioLector I device.
       - Requires data files from the machine as well as a file matching well IDs with sample names.
       - Will automatically compute mean and standard deviation across replicates.
       - Converts measured signal into CDW [g/L] values based on conversion factors.
       - Calculates growth rates for each growth curve.
       - Plots biomass and growth rates in a tiledlayout plot.
  6. [HPLC_Analysis_Shimadzu.m](MATLAB_Data_Analysis/HPLC_Analysis_Shimadzu.m):
       - Analyses HPLC peak tables collected from a SIL-20ACHT HPLC machine (Shimadzu, Kyoto, Japan).
       - Requires HPLC data files as well as an excel file containing the order of samples that were run.
       - Will plot standard curves and save curve parameters into a matlab variable that can be loaded later.
       - Will calculate concentrations of unknown samples based on previously measured and analysed standards.
  7. [batch_triplicate_w_HPLC.m](MATLAB_Data_Analysis/batch_triplicate_w_HPLC.m):
       - 
