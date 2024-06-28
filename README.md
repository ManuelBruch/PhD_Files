# About This Folder
This repository contains MATLAB and Python code files used in the work accompanying the thesis entitled **"Metabolic Engineering of Bacteria to Convert CO<sub>2</sub> to Products of Value"**. The structure and use of the files found here will be explaine briefly in this document.

## A General Remark
Many of the coding files here contain references to file locations and file names. These are purely left in for demonstration purposes, if these files were to be used by the Reader, they would need to point the scripts to wherever they might have equivalent files in their working environment.
<br>
<br>

## [MATLAB_Data_Analysis](MATLAB_Data_Analysis)
This folder contains general files for data analysis using MATLAB. All files are at least compatible with version R2021a or newer, older versions have not been checked. The purpose of each file will be briefly pointed out here.
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
       - Plots growth curve triplicates together with measured HPLC data from batch fermentations.
       - Requires respective data files.
  8. [conversion_factor_OD_to_CDW.m](MATLAB_Data_Analysis/conversion_factor_OD_to_CDW.m):
       - Calculates a conversion factor from OD<sub>600</sub> measurements to CDW in [g/L].
  9. [flask_triplicate_w_HPLC.m](MATLAB_Data_Analysis/flask_triplicate_w_HPLC.m):
       - Plots growth curve triplicates together with measured HPLC data from shake flask cultivations.
       - Requires respective data files.

### The [functions](MATLAB_Data_Analysis/functions) Folder
This folder contains some general helper functions written to streamline the code in the other MATLAB files.
  1. [HPLC_Sample_Analysis_Shimadzu.m](MATLAB_Data_Analysis/functions/HPLC_Sample_Analysis_Shimadzu.m):
       - Used in the [HPLC_Analysis_Shimadzu.m](MATLAB_Data_Analysis/HPLC_Analysis_Shimadzu.m) file to analyse unknown samples and create a time series of values.
  2. [HPLC_Sample_Analysis_Shimadzu_noT.m](MATLAB_Data_Analysis/functions/HPLC_Sample_Analysis_Shimadzu_noT.m):
       - Used in the [HPLC_Analysis_Shimadzu.m](MATLAB_Data_Analysis/HPLC_Analysis_Shimadzu.m) file to analyse unknown samples when they are not coming from a time series (such as a growth curve).
  3. [HPLC_Standard_Analysis_Shimadzu.m](MATLAB_Data_Analysis/functions/HPLC_Standard_Analysis_Shimadzu.m):
       - Used in the [HPLC_Analysis_Shimadzu.m](MATLAB_Data_Analysis/HPLC_Analysis_Shimadzu.m) file to analyse standard curves of known concentrations.
  4. [myplot.m](MATLAB_Data_Analysis/functions/myplot.m):
       - Older function to standardise plotting across various files.
       - Used in some of the older scripts
       - Requires input parameters to be added in a specific order.
  5. [myNewPlot.m](MATLAB_Data_Analysis/functions/myNewPlot.m):
       - Newer function for plotting standardisation in various files.
       - Allows for input parameters to be handled as name/value pair arguments and therefore does not require a specific order.
       - Used in newer files.
  6. [slidingWindow.m](MATLAB_Data_Analysis/functions/slidingWindow.m):
       - Sliding Window algorithm to calculate linear regression lines across the given window.
       - Returns parameters of the linear regression, such as slope, y-intercepts, error estimates and x-value ranges.
       - Used exclusively in this work for the calculation of growth rates.
<br>

## [Proteomics_Data_Analysis](Proteomics_Data_Analysis)
This folder contains files associated with the proteomics analysis conducted in this work.
### [Perseus_Analysis](Proteomics_Data_Analysis/Perseus_Analysis)
The file doing in which the initial proteome analysis and data transformation was performed is found in [Perseus_proteomics_analysis.sps](Proteomics_Data_Analysis/Perseus_Analysis/Perseus_proteomics_analysis.sps). A snapshot of the workflow is shown in [Perseus_workflow.png](Proteomics_Data_Analysis/Perseus_Analysis/Perseus_workflow.png).
### [PRISM_Analysis](Proteomics_Data_Analysis/PRISM_Analysis)
This folder contains the [GraphPad PRISM file](Proteomics_Data_Analysis/PRISM_Analysis/20230725statistical_analysis_proteins_in_both_strains_post_imputation.prism), used to conduct statistical analysis of the proteomics data. Additionally, folders with extracted tables for the analysis results for proteins detected [only in _C. necator_ H16](Proteomics_Data_Analysis/PRISM_Analysis/H16_only), [only in _C. necator_ ALE26](Proteomics_Data_Analysis/PRISM_Analysis/ALE26_only) or in [both strains](Proteomics_Data_Analysis/PRISM_Analysis/both) are provided.
### [Python_code](Proteomics_Data_Analysis/Python_code)
This folder contains the python code used in analysing and graphing the proteomics results. The purpose of each file will briefly be laid out here. <br>
Required Python packages: _biopython_, _gurobi_ (with an active license), _Ipython_ when working with Spyder, _matplotlib_, _numpy_, _pandas_, _scipy_, _seaborn_
  
  1. [differential_expression_analysis.py](Proteomics_Data_Analysis/Python_code/differential_expression_analysis.py):
       - Reads in data after exporting them from Perseus.
       - uses the Protein IDs assigned in the MaxQuant search to search NCBI for associated entries and fetch genome locations of the corresponding genes.
       - writes these locations (chromosome 1, chromosome 2 or pHG1) into the original dataframe and saves it.
  2. [proteomics_PCA.py](Proteomics_Data_Analysis/Python_code/proteomics_PCA.py):
       - Runs a principal component analysis (PCA) on the data extracted from Perseus.
  3. [format_PRISM_results.py](Proteomics_Data_Analysis/Python_code/format_PRISM_results.py):
       - Reformats the results from the GraphPad PRISM Sidak analysis into a simpler table and extracts -log(adjusted-p-values) and log(fold changes).
  4. [build_pool_table_plot_dataframe.py](Proteomics_Data_Analysis/Python_code/build_pool_table_plot_dataframe.py):
       - Reads in the reformatted PRISM results and based on the Uniprot protein IDs fetches KEGG IDs and pathways associated with each protein based on KEGG and a list manually downloaded from Uniprot.
       - Creates dataframes for each of the original three tables (detected in: _C. necator_ H16 only, _C. necator_ ALE26 only or both) and one merged dataframe of the three individual ones.
  5. [group_pathways_in_meta_poolTables.py](Proteomics_Data_Analysis/Python_code/group_pathways_in_meta_poolTables.py):
       - Defines some overarching pathways based on the KEGG BRITE hierarchy to add to the dataframes created in [build_pool_table_plot_dataframe.py](Proteomics_Data_Analysis/Python_code/build_pool_table_plot_dataframe.py).
       - Requires the BRITE hierarchy as an additional input (can be downloaded from the KEGG website, strain identifier: reh).
  6. [proteomics_poolTablePlot_plotting.py](Proteomics_Data_Analysis/Python_code/proteomics_poolTablePlot_plotting.py):
       - Creates plots displaying the overarching pathway for each protein (y-axis) against it's log(fold change) (x-axis) for the comparison between two growth conditions in a scatter plot.
       - colour-codes each marker based on the associated adjusted p-value.
       - plots proteins detected in _C. necator_ H16 only, _C. necator_ ALE26 only or both strains on seperate plots.
  7. [proteomics_poolTablePlot_plotting_accounting_for_nan.py](Proteomics_Data_Analysis/Python_code/proteomics_poolTablePlot_plotting_accounting_for_nan.py):
       - creates a similar plot to the one in [proteomics_poolTablePlot_plotting.py](Proteomics_Data_Analysis/Python_code/proteomics_poolTablePlot_plotting.py).
       - combines all three comparisons by plotting proteins detected in only one strain at arbitrarily large or small fold changes to seperate them from those found in both strains.
  8. [proteomics_heat_maps_for_pathways.py](Proteomics_Data_Analysis/Python_code/proteomics_heat_maps_for_pathways.py):
       - Creates heatmaps for the protein expression in each condition after imputation and averaging the expression data.
       - Creates individual heat maps for each of the overarching pathways determined before.
