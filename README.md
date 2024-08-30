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
  10. [DASbox_conti_analysis.m](MATLAB_Data_Analysis/DASbox_conti_analysis.m):
       - Generates plots using previously analysed offline biomass data and online-measured data from the DASware Control software from Eppendorf to display continuous fermentation profiles (run in DASbox Mini Fermenter units).
       - Required the respective data files.

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
  6. [PHB_GC_evaluation.m](MATLAB_Data_Analysis/functions/PHB_GC_evaluation.m):
       - Function to evaulate GC data based on user-provided table.
       - Uses formula developed in this lab to quantify PHA monomer mass based on the peak areas of a sample peak and an internal standard (methyl benzoate).
  7. [slidingWindow.m](MATLAB_Data_Analysis/functions/slidingWindow.m):
       - Sliding Window algorithm to calculate linear regression lines across the given window.
       - Returns parameters of the linear regression, such as slope, y-intercepts, error estimates and x-value ranges.
       - Used exclusively in this work for the calculation of growth rates.
<br>

## [Proteomics_Data_Analysis](Proteomics_Data_Analysis)
This folder contains files associated with the proteomics analysis conducted in this work.
The excel file [poolTable_imports_of_python_dataframes.xlsx](Proteomics_Data_Analysis/poolTable_imports_of_python_dataframes.xlsx) contains tables with the statistical results for differential protein expression between all tested conditions. It also features the pathway associations of each of the identified proteins. This table was the basis from which the pool-table plots were created.
### [Perseus_Analysis](Proteomics_Data_Analysis/Perseus_Analysis)
The file doing in which the initial proteome analysis and data transformation was performed is found in [Perseus_proteomics_analysis.sps](Proteomics_Data_Analysis/Perseus_Analysis/Perseus_proteomics_analysis.sps). A snapshot of the workflow is shown in [Perseus_workflow.png](Proteomics_Data_Analysis/Perseus_Analysis/Perseus_workflow.png).
### [PRISM_Analysis](Proteomics_Data_Analysis/PRISM_Analysis)
This folder contains the [GraphPad PRISM file](Proteomics_Data_Analysis/PRISM_Analysis/20230725statistical_analysis_proteins_in_both_strains_post_imputation.prism), used to conduct statistical analysis of the proteomics data. Additionally, folders with extracted tables for the analysis results for proteins detected [only in _C. necator_ H16](Proteomics_Data_Analysis/PRISM_Analysis/H16_only), [only in _C. necator_ ALE26](Proteomics_Data_Analysis/PRISM_Analysis/ALE26_only) or in [both strains](Proteomics_Data_Analysis/PRISM_Analysis/both) are provided.
### [Python_code](Proteomics_Data_Analysis/Python_code)
This folder contains the python code used in analysing and graphing the proteomics results. The purpose of each file will briefly be laid out here. <br>
Required Python packages: _biopython_, _gurobi_ (with an active license), _matplotlib_, _numpy_, _pandas_, _scipy_, _seaborn_
  
  1. [differential_expression_analysis.py](Proteomics_Data_Analysis/Python_code/differential_expression_analysis.py):
       - Reads in data after exporting them from Perseus.
       - Uses the Protein IDs assigned in the MaxQuant search to search NCBI for associated entries and fetch genome locations of the corresponding genes.
       - Writes these locations (chromosome 1, chromosome 2 or pHG1) into the original dataframe and saves it.
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
       - Colour-codes each marker based on the associated adjusted p-value.
       - Plots proteins detected in _C. necator_ H16 only, _C. necator_ ALE26 only or both strains on seperate plots.
  7. [proteomics_poolTablePlot_plotting_accounting_for_nan.py](Proteomics_Data_Analysis/Python_code/proteomics_poolTablePlot_plotting_accounting_for_nan.py):
       - Creates a similar plot to the one in [proteomics_poolTablePlot_plotting.py](Proteomics_Data_Analysis/Python_code/proteomics_poolTablePlot_plotting.py).
       - Combines all three comparisons by plotting proteins detected in only one strain at arbitrarily large or small fold changes to seperate them from those found in both strains.
  8. [proteomics_heat_maps_for_pathways.py](Proteomics_Data_Analysis/Python_code/proteomics_heat_maps_for_pathways.py):
       - Creates heatmaps for the protein expression in each condition after imputation and averaging the expression data.
       - Creates individual heat maps for each of the overarching pathways determined before.
  9. [convert_proteomics_to_mmol_per_gCDW.py](Proteomics_Data_Analysis/Python_code/convert_proteomics_to_mmol_per_gCDW.py):
       - Converts imputed LFQ data into protein abundancies normalised by the cell dry weight concentration in each original sample.
       - Uses hard-coded sample quantities like volume and CDW.
       - Calculates mean and standard deviation for each condition and saves the resulting dataframe.
  10. [proteome_quantity_per_reaction_in_C_necator_central_carbon_metabolism.py](Proteomics_Data_Analysis/Python_code/proteome_quantity_per_reaction_in_C_necator_central_carbon_metabolism.py):
       - Uses gene-reaction relations from the [R-Notebook around the Cupriavidus proteome](https://github.com/m-jahn/R-notebook-ralstonia-proteome/blob/main/data/input/model_reactions.csv) to match protein abundancies with their respective reactions.
       - Plots bar charts for each reaction with all involved protein abundancies (mean and standard deviation) for each tested condition.
### [Additional_Figures](Proteomics_Data_Analysis/Additional_Figures)
This folder contains additional figures to the thesis that were either too many or too large to fit them in the printed document.
  - [reaction_resolved_protein_expression](Proteomics_Data_Analysis/Additional_Figures/reaction_resolved_protein_expression):
      - Contains Figures showing protein abundancies associated with the reactions of the central carbon metabolism in _C. necator_
      - [C_necator_protein_expression_central_carbon_metabolism.png](Proteomics_Data_Analysis/Additional_Figures/reaction_resolved_protein_expression/C_necator_protein_expression_central_carbon_metabolism.png):
        - Shows all proteins detected in the proteomics analysis in conjunction with the reaction they are associated with.
        - Reactions having no associated proteins detected are only represented by their reaction ID.
        - This figure is referenced in the thesis in chapter 3, section 3.4.2.
      - [plots](Proteomics_Data_Analysis/Additional_Figures/reaction_resolved_protein_expression/plots):
        - Contains the individual plots used in [C_necator_protein_expression_central_carbon_metabolism.png](Proteomics_Data_Analysis/Additional_Figures/reaction_resolved_protein_expression/C_necator_protein_expression_central_carbon_metabolism.png).
  - [heatmaps](Proteomics_Data_Analysis/Additional_Figures/heatmaps):
      - Contains heatmaps depicting the average protein expression in each condition for each protein in each meta-pathway in _C. necator_ in heatmaps.
<br>

## [Genomic_Analysis_Graphics](Genomic_Analysis_Graphics)
This folder contains additional graphics supplementary to the genomic analyses described in the thesis.
### [First_Attempt](Genomic_Analysis_Graphics/First_Attempt)
This folder contains results for the first attempt of the analysis with faulty data provided by the sequencing service.
#### [Bandage_graphs](Genomic_Analysis_Graphics/First_Attempt/Bandage_graphs)
This folder contains Bandage graphs created using the Bandage Image tool in Galaxy.
  - [C_necator_ALE31_Bandage_graph.svg](Genomic_Analysis_Graphics/First_Attempt/Bandage_graphs/C_necator_ALE31_Bandage_graph.svg):
      - Bandage graph created in the Galaxy project to display interactions amongst scaffolds in the assembly from the sample to sequence _C. necator_ ALE31.
  - [C_necator_ALE42_Bandage_graph.svg](Genomic_Analysis_Graphics/First_Attempt/Bandage_graphs/C_necator_ALE42_Bandage_graph.svg):
      - Bandage graph created in the Galaxy project to display interactions amongst scaffolds in the assembly from the sample to sequence _C. necator_ ALE42.
#### [Circos_plots](Genomic_Analysis_Graphics/First_Attempt/Circos_plots)
This folder contains circos plots created using the Quast tool in Galaxy.
  - [circos_plot_H16.png](Genomic_Analysis_Graphics/First_Attempt/Circos_plots/circos_plot_H16.png):
      - Circos plot for the alignment of assembled scaffolds for _C. necator_ H16 sequencing reads with the [reference genome](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA531660).
  - [circos_plot_ALE26.png](Genomic_Analysis_Graphics/First_Attempt/Circos_plots/circos_plot_ALE26.png):
      - Circos plot for the alignment of assembled scaffolds for _C. necator_ ALE26 sequencing reads with the [reference genome](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA531660).
  - [circos_plot_ALE42.png](Genomic_Analysis_Graphics/First_Attempt/Circos_plots/circos_plot_ALE42.png):
      - Circos plot for the alignment of assembled scaffolds for _C. necator_ ALE42 sequencing reads with the [reference genome](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA531660).
   
### [Second_Attempt_with_corrected_sequences](Genomic_Analysis_Graphics/Second_Attempt_with_corrected_sequences)
This folder contains additional graphs from the corrected sequence analysis
#### [Bandage_graphs](Genomic_Analysis_Graphics/Second_Attempt_with_corrected_sequences)
This folder contains Bandage graphs created using the Bandage Image tool in Galaxy.
  - [C_necator_ALE25_Bandage_graph.svg](Genomic_Analysis_Graphics/Second_Attempt_with_corrected_sequences/Bandage_Graphs/C_necator_ALE25_Bandage_graph.svg):
      - Bandage graph created in the Galaxy project to display interactions amongst scaffolds in the assembly from the sample to sequence _C. necator_ ALE25.
  - [C_necator_ALE31_Bandage_graph.svg](Genomic_Analysis_Graphics/Second_Attempt_with_corrected_sequences/Bandage_Graphs/C_necator_ALE31_Bandage_graph.svg):
      - Bandage graph created in the Galaxy project to display interactions amongst scaffolds in the assembly from the sample to sequence _C. necator_ ALE31.
  - [C_necator_ALE42_Bandage_graph.svg](Genomic_Analysis_Graphics/Second_Attempt_with_corrected_sequences/Bandage_Graphs/C_necator_ALE42_Bandage_graph.svg):
      - Bandage graph created in the Galaxy project to display interactions amongst scaffolds in the assembly from the sample to sequence _C. necator_ ALE42.
  - [C_necator_ALE25_Bandage_graph_with_filtered_scafolds.svg](Genomic_Analysis_Graphics/Second_Attempt_with_corrected_sequences/Bandage_Graphs/C_necator_ALE25_Bandage_graph_with_filtered_scafolds.svg):
      - Bandage graph created in the Galaxy project to display interactions amongst scaffolds in the assembly from the sample to sequence _C. necator_ ALE25.
      - This version was created after filtering out scaffolds with a low coverage.
  - [C_necator_ALE31_Bandage_graph_with_filtered_scafolds.svg](Genomic_Analysis_Graphics/Second_Attempt_with_corrected_sequences/Bandage_Graphs/C_necator_ALE31_Bandage_graph_with_filtered_scafolds.svg):
      - Bandage graph created in the Galaxy project to display interactions amongst scaffolds in the assembly from the sample to sequence _C. necator_ ALE31.
      - This version was created after filtering out scaffolds with a low coverage.
  - [C_necator_ALE42_Bandage_graph_with_filtered_scafolds.svg](Genomic_Analysis_Graphics/Second_Attempt_with_corrected_sequences/Bandage_Graphs/C_necator_ALE42_Bandage_graph_with_filtered_scafolds.svg):
      - Bandage graph created in the Galaxy project to display interactions amongst scaffolds in the assembly from the sample to sequence _C. necator_ ALE42.
      - This version was created after filtering out scaffolds with a low coverage.
#### [Circos_plots](Genomic_Analysis_Graphics/Second_Attempt_with_corrected_sequences/Circos_plots)
This folder contains circos plots created using the Quast tool in Galaxy.
  - [circos_plot_ALE25.png](Genomic_Analysis_Graphics/Second_Attempt_with_corrected_sequences/Circos_plots/circos_plot_ALE25.png):
      - Circos plot for the alignment of assembled scaffolds for _C. necator_ ALE25 sequencing reads with the [reference genome](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA531660).
  - [circos_plot_ALE31.png](Genomic_Analysis_Graphics/Second_Attempt_with_corrected_sequences/Circos_plots/circos_plot_ALE31.png):
      - Circos plot for the alignment of assembled scaffolds for _C. necator_ ALE31 sequencing reads with the [reference genome](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA531660).
  - [circos_plot_ALE42.png](Genomic_Analysis_Graphics/Second_Attempt_with_corrected_sequences/Circos_plots/circos_plot_ALE42.png):
      - Circos plot for the alignment of assembled scaffolds for _C. necator_ ALE42 sequencing reads with the [reference genome](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA531660).
  - [circos_plot_ALE25_with_filtered_scaffolds.png](Genomic_Analysis_Graphics/Second_Attempt_with_corrected_sequences/Circos_plots/circos_plot_ALE25_with_filtered_scaffolds.png):
      - Circos plot for the alignment of assembled scaffolds for _C. necator_ ALE25 sequencing reads with the [reference genome](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA531660).
      - This version was created after filtering out scaffolds with a low coverage.
  - [circos_plot_ALE31_with_filtered_scaffolds.png](Genomic_Analysis_Graphics/Second_Attempt_with_corrected_sequences/Circos_plots/circos_plot_ALE31_with_filtered_scaffolds.png):
      - Circos plot for the alignment of assembled scaffolds for _C. necator_ ALE31 sequencing reads with the [reference genome](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA531660).
      - This version was created after filtering out scaffolds with a low coverage.
  - [circos_plot_ALE42_with_filtered_scaffolds.png](Genomic_Analysis_Graphics/Second_Attempt_with_corrected_sequences/Circos_plots/circos_plot_ALE42_with_filtered_scaffolds.png):
      - Circos plot for the alignment of assembled scaffolds for _C. necator_ ALE42 sequencing reads with the [reference genome](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA531660).
      - This version was created after filtering out scaffolds with a low coverage.
<br>

## [Resource_Balance_Analysis](Resource_Balance_Analysis)
This folder contains files for RBA on _C. necator_. The model used in these simulations is the one forked from [Michael Jahn's repository](https://github.com/m-jahn/Bacterial-RBA-models/tree/master/Ralstonia-eutropha-H16) as part of [Jahn et al., 2021](https://doi.org/10.7554/eLife.69019).
Required Python packages: _biopython_, _gurobi_ (with an active license), _matplotlib_, _numpy_, _pandas_, _scipy_, _seaborn_, _RBApy_ (requires python 3.7 or older)
  - [JMMmedium.tsv](Resource_Balance_Analysis/JMMmedium.tsv):
      - A tab-seperated file containing the concentrations of metabolites in the growth medium.
  - [original_medium.tsv](Resource_Balance_Analysis/original_medium.tsv):
      - The growth medium composition used in the original model by Michael Jahn.
  - [solve_model_V2.py](Resource_Balance_Analysis/solve_model_V2.py):
      - This is a mostly unaltered copy of the [solve_model.py](https://github.com/m-jahn/Bacterial-RBA-models/blob/master/Ralstonia-eutropha-H16/solve_model.py) function by Michael Jahn, with some minor alterations for file paths.
      - Copied the write_proteins() method from the RBApy module into this script as a standalone function to alter the output to a full output, rather than just the proteins that carry a non-zero amount of flux to get a full list of the proteins/cellular machinery in the model.
  - [Proteomics_RBA_comparison.py](Resource_Balance_Analysis/Proteomics_RBA_comparison.py):
      - Script to compare predicted proteins from RBA analysis and measured data from proteomics.
      - Plots correlation graphs.
      - Requires outputs from [solve_model_V2.py](Resource_Balance_Analysis/solve_model_V2.py) and [convert_proteomics_to_mmol_per_gCDW.py](Proteomics_Data_Analysis/Python_code/convert_proteomics_to_mmol_per_gCDW.py), as well as UniProt search results from manually looking up UniProt ID and Gene Names found in the RBA model and a file containing the gene-reaction.relation rules.
### [model](Resource_Balance_Analysis/model)
This folder contains the RBA model files created by Michael Jahn.
### [simulation](Resource_Balance_Analysis/simulation)
This folder is the output folder for simulations run on the model. It also contains the [substrate_input.csv](Resource_Balance_Analysis/simulation/substrate_input.csv) file used to hand growth medium parameters to the model.
### [RBVA](Resource_Balance_Analysis/RBVA)
This folder contains scripts created to run resource balance variability analysis (RBVA) using the [_RBAtools_](https://rba.inrae.fr/rbatools.html) python module ([Bodeit et al., 2023](https://doi.org/10.1093/bioadv/vbad056)).
  - [RBVA_MJahn_model.ipynb](Resource_Balance_Analysis/RBVA/RBVA_MJahn_model.ipynb):
      - Jupyter notebook to conduct RBVA.
      - Reads the same files as [solve_model_V2.py](Resource_Balance_Analysis/solve_model_V2.py) and allows for manual selection of the condition to test.
      - Allows for either using the automatically calculated maximum growth rate or for a user-defined input.
      - Performs FVA on the reactions present in the model and records predicted protein concentrations.
  - [RBVA_proteomics_comparison.py](Resource_Balance_Analysis/RBVA/RBVA_proteomics_comparison.py):
      - Similar to [Proteomics_RBA_comparison.py](Resource_Balance_Analysis/Proteomics_RBA_comparison.py) but considers the RBVA results.
      - Also computes some descriptive statistics of the RBVA data, such as mean, median, min and max concentration for each protein.
  - [RBVA_data_variability.py](Resource_Balance_Analysis/RBVA/RBVA_data_variability.py):
      - Plots range for each predicted protein that carries flux in any of the RBVA conditions as well as mean measured protein concentrations for each strain for a given growth condition.
