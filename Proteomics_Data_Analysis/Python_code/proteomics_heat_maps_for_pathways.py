# -*- coding: utf-8 -*-
"""
Created on Wed Nov 15 14:13:40 2023

@author: Manuel Bruch
"""

#%% clean up
from IPython import get_ipython
get_ipython().magic('reset -sf') # clears all variables

#%% import statements
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


#%%
def findOccurrences(s, ch):
    return [i for i, letter in enumerate(s) if letter == ch]


#%% load in data
curFolder = os.getcwd()
backSlash = findOccurrences(curFolder, os.sep)
baseFolder = curFolder[:backSlash[-1]]
baseDataFolder = os.path.join(baseFolder, 'data', 'proteomics', '20230620',
                              'data_analysis')

# load in table from pool table plots with meta pathways
fullFileName = os.path.join(baseDataFolder, 'poolTablePlots',
                            'poolTable_df_merged_metaPWs.csv')
pwDF = pd.read_csv(fullFileName, sep="\t")
# only keep columns with protein info
cols2keep = ['UniProtIDs', #'EntryName', 'ProteinName', 'geneLocations',
             # 'UniProt_Pathways', 'KEGG_IDs', 'KEGG_Pathways', 
             # 'KEGG_Pathway_Modules',
             'KEGG_and_UniProt_Meta_Pathways']
pwDF = pwDF[cols2keep]

# load in table from RBA comparison for mean expression level
fullFileName = os.path.join(baseDataFolder, 'model_comparisons',
                            'mean_allImputedData_mmol_p_g_CDW_not_BCA.csv')
expressionDF = pd.read_csv(fullFileName, sep="\t")

#%% join the dataframes
expressionDF.rename(columns = {'Protein IDs': 'UniProtIDs'}, inplace = True)

allData = pwDF.merge(expressionDF, how = 'left', on = 'UniProtIDs')
# remove standard deviation columns
cols2keep = list(allData.columns)
cols2keep = [i for i in cols2keep if '_std_' not in i]
allData = allData[cols2keep]
backUpData = allData.copy()

#%% normalise the data in each row by the maximum value
dataCols = [i for i in cols2keep if '_mean_' in i]
allData[dataCols] = allData[dataCols].div(allData[dataCols].max(axis = 1), axis = 0)

allData.set_index('UniProtIDs', inplace = True)
# rename data containing columns
for col in list(allData.columns):
    if '_mean_' in col:
        allData.rename(columns = {col: col[:-15]}, inplace = True)

# change column order
allData = allData[['KEGG_and_UniProt_Meta_Pathways',
                  'H16_80_full', 'H16_80_lim', 'H16_150_lim', 
                  'ALE26_80_full', 'ALE26_80_lim', 'ALE26_150_lim']]

#%% have a dataframe with the not-normalised data
backUpData.set_index('UniProtIDs', inplace = True)
# rename data containing columns
for col in list(backUpData.columns):
    if '_mean_' in col:
        backUpData.rename(columns = {col: col[:-15]}, inplace = True)

# change column order
backUpData = backUpData[['KEGG_and_UniProt_Meta_Pathways',
                         'H16_80_full', 'H16_80_lim', 'H16_150_lim', 
                         'ALE26_80_full', 'ALE26_80_lim', 'ALE26_150_lim']]

backUpData_sum = backUpData.sum(numeric_only = True)

#%% get individual pathways for dictionary
tableHeaders = list(allData.columns)
metaPWHeader = [i for i in tableHeaders if 'Meta_Pathways' in i][0]
pathways = list(allData[metaPWHeader])
pathways_split = [item.split('//') for item in pathways]
pathways_flat = [item for l in pathways_split for item in l]
pathways_set = list(set(pathways_flat))

#%% create dictionary of dataframes for each pathway
pWsepDFs = {}
# remove meta pathways column
cols2keep = [i for i in tableHeaders if 'Meta_Pathways' not in i]

for pW in pathways_set:
    pWsepDFs[pW] = allData[allData[metaPWHeader].str.contains(pW)]
    pWsepDFs[pW] = pWsepDFs[pW][cols2keep]#.transpose()

#%% plot heat maps
saveFilePath = os.path.join(baseDataFolder, 'pathway_heatmaps')
for pW in pWsepDFs:
    rowNums = len(pWsepDFs[pW].index)
    figsizeY = rowNums*(0.2) + 0.2
    rowCols = len(pWsepDFs[pW].columns)
    figsizeX = rowNums*(0.15) + 0.2
    # plt.figure(figsize = (12,figsizeY))
    fig, axs = plt.subplots(1, 2,
                             gridspec_kw={#'height_ratios': [50, 1], 
                                          'width_ratios': [50, 1], 
                                          'wspace': 0.1}, 
                             sharex='col',
                             # sharey='row',
                             figsize = (12 ,figsizeY))
    curPlot = sns.heatmap(pWsepDFs[pW], cmap="crest",# yticklabels = True,
                          ax=axs[0], cbar_ax=axs[1])
    plt.axes(axs[0])
    plt.yticks(rotation = 0) 
    curPlot.set(ylabel = None)
    plt.title(pW)
    
    #%%% save the figure
    figName = pW.replace(':','')
    figName = figName.replace(' ','_')
    
    # figureSavePath = os.path.join(saveFilePath, 
    #                               figName + '.svg')
    # curPlot.figure.savefig(figureSavePath, format='svg')
    # figureSavePath = os.path.join(saveFilePath, 
    #                               figName + '.png')
    # curPlot.figure.savefig(figureSavePath, bbox_inches='tight', dpi = 500)