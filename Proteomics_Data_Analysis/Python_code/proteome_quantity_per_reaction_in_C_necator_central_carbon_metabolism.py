# -*- coding: utf-8 -*-
"""
Created on Mon Jul  1 10:45:42 2024

@author: Manuel Bruch
"""

#%% clean up
from IPython import get_ipython
get_ipython().magic('reset -sf') # clears all variables

#%% import statements
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

#%% set parameters


#%% function definitions
def findOccurrences(s, ch):
    return [i for i, letter in enumerate(s) if letter == ch]


def flatten(l):
    return [item for sublist in l for item in sublist]


#%% load in data
curFolder = os.getcwd()
backSlash = findOccurrences(curFolder, os.sep)
baseFolder = curFolder[:backSlash[-2]]
fileLocation = os.path.join(baseFolder, 'Cupriavidus-necator-modelling',
                            'data', 'proteomics', '20230620', 'data_analysis')
subFolder = ['model_comparisons',
             'central_carbon_metabolism_protein_quantities']
figSaveFolder = os.path.join(fileLocation, subFolder[1], 'plots')
protFileName = 'mean_allImputedData_mmol_p_g_CDW_BCA.csv'
# contains mean and stddev in mmol(protein)/g(CDW) for each detected protein
# for each condition after imputation
protData = pd.read_csv(os.path.join(fileLocation, subFolder[0], protFileName),
                       sep="\t")

"""
# the following are likely unneccessary:
gprFile = os.path.join(baseFolder, 'R-notebook-ralstonia-proteome', 'data',
#                        'input', 'model_reactions.csv')
# # original reaction set of Michael's R-notebook
gprDf = pd.read_csv(gprFile, index_col = 0)

gprUniProtIDFile = os.path.join(fileLocation, subFolder[0],
#                                 'GPR_gene_IDs_idmapping_2023_09_04.tsv')
# # maps each gene with the corresponding UniProt ID
gprUniProtIDs = pd.read_csv(gprUniProtIDFile, sep='\t')
"""

gprDfUniProtSaveFile = os.path.join(fileLocation, subFolder[0],
                                    'GPR_UniProt_IDs.csv')
# like Michel's file, but with the gene IDs replaced with UniProt protein IDs
gprDfUniProt = pd.read_csv(gprDfUniProtSaveFile, sep = '\t')

# load in file relating Michal's reaction IDs to the ones used in the graphic
rnNameFileName = 'CCM_reactions.csv'
rnNames = pd.read_csv(os.path.join(fileLocation, subFolder[1],
                                    rnNameFileName))

#%% get rows of gpr df relevant to this script
# clean up rnNames dataframe by removing nan-containing rows
rnNames = rnNames.dropna()
# converting Michael's reaction IDs into list and selecting df rows
modelRnIDs = list(rnNames[list(rnNames.columns)[1]])
mask = gprDfUniProt['reaction_id'].isin(modelRnIDs)
CCMgprs = gprDfUniProt[mask]
cols2keep = ['reaction_id', 'reaction_name', 'genes']
CCMgprsRed = CCMgprs[cols2keep]

CCMgprsRed_myIDs = pd.merge(CCMgprsRed, rnNames, left_on = 'reaction_id',
                            right_on = 'ID_Michael')

#%% build a dictionary
# using my reaction IDs as keys and the corresponding proteins and rows with
# quantities as values
dataDict = dict()
for index, row in CCMgprsRed_myIDs.iterrows():
    dataDict[row['ID_Image']] = dict()
    # get all the involved proteins as a list
    if isinstance(row['genes'], str):
        curProts = row['genes'].split(', ')
    elif np.isnan(row['genes']):
        curProts = []
    dataDict[row['ID_Image']]['proteins'] = curProts
    
    # use list as mask to select rows from protData DF
    mask = protData['Protein IDs'].isin(curProts)
    dataDict[row['ID_Image']]['expression'] = protData[mask]
    # print('# detected proteins in ' + row['ID_Image'] + ': ' + str(len(protData[mask].index)))

#%% reshape dataframes for plotting
# possibly try this: https://stackoverflow.com/questions/75240070/how-can-custom-errorbars-be-aligned-on-grouped-bars

# rename conditions with roman numerals
condDictMean = {'ALE26_80_lim_mean_mmol_gCDW':  'V',
                'H16_150_lim_mean_mmol_gCDW':   'III',
                'ALE26_80_full_mean_mmol_gCDW': 'IV',
                'H16_80_full_mean_mmol_gCDW':   'I',
                'ALE26_150_lim_mean_mmol_gCDW': 'VI',
                'H16_80_lim_mean_mmol_gCDW':    'II'}

condDictStd = {'ALE26_80_lim_std_mmol_gCDW':  'V',
                'H16_150_lim_std_mmol_gCDW':   'III',
                'ALE26_80_full_std_mmol_gCDW': 'IV',
                'H16_80_full_std_mmol_gCDW':   'I',
                'ALE26_150_lim_std_mmol_gCDW': 'VI',
                'H16_80_lim_std_mmol_gCDW':    'II'}

for key in dataDict:
    exprDF = dataDict[key]['expression']
    protNum = 0
    if not len(exprDF) == 0:
        for index, row in exprDF.iterrows():
            meanCols = [i for i in list(row.index) if 'mean' in i]
            meanDF = row[meanCols]
            meanDF = meanDF.to_frame().reset_index()
            oldNames = list(meanDF.columns)
            meanDF.rename(columns = {oldNames[0]: 'cond', oldNames[1]: row['Protein IDs']},
                          inplace = True)
            meanDF.replace({'cond': condDictMean}, inplace = True)
            # meanDF['Protein'] = [row['Protein IDs']] * len(meanDF)
            
            stdCols = [i for i in list(row.index) if 'std' in i]
            stdDF = row[stdCols]
            stdDF = stdDF.to_frame().reset_index()
            oldNames = list(stdDF.columns)
            stdDF.rename(columns = {oldNames[0]: 'cond', oldNames[1]: row['Protein IDs'] + '_std'},
                         inplace = True)
            stdDF.replace({'cond': condDictStd}, inplace = True)
            # stdDF['Protein'] = [row['Protein IDs']] * len(stdDF)
            rowDF = pd.merge(meanDF, stdDF, how = 'inner')
            if protNum == 0:
                rxnDF = rowDF
            else:
                rxnDF = pd.merge(rxnDF, rowDF, how = 'inner')
            protNum += 1
        dataDict[key]['orgData'] = rxnDF.sort_values(by = 'cond').reset_index(drop = True)
    else:
        dataDict[key]['orgData'] = 'Proteins not detected'
    
#%% plot data

font1 = {'family': 'arial', 'color': 'black', 'size': 10} # for legends, axes
font2 = {'family': 'arial', 'color': 'black', 'weight': 'bold', 'size': 14} # for titles
for key in dataDict:
    if not isinstance(dataDict[key]['orgData'], str):
        data = dataDict[key]['orgData']
        data = data.set_index('cond') # automatically sets "cond" as x-tick labels
        
        cols = list(data.columns)
        meanCols = [i for i in cols if not 'std' in i and not 'cond' in i]
        stdCols = [i for i in cols if 'std' in i or 'cond' in i]
        
        yerr = data[stdCols]  # need to rename them to match data columns
        stColsRenamed = [sub.replace('_std', '') for sub in stdCols]
        yerr.columns = stColsRenamed
        
        ax = data.plot(kind = 'bar', y = meanCols, yerr = yerr,
                       rot = 0, fontsize = 10, figsize = (5, 4))
        ax.legend(bbox_to_anchor=(1.0, 1.0), fontsize = 10)
        ax.set_title(key, fontdict = font2)
        ax.set_ylim(bottom = 0)
        ax.set_xlabel('')
        ax.set_ylabel('c(Protein) [mmol/g(CDW)]', fontdict = font1)
        
        # save the figure
        figureSavePath = os.path.join(figSaveFolder, key + '.svg')
        ax.figure.savefig(figureSavePath, format = 'svg')
        figureSavePath = os.path.join(figSaveFolder, key + '.png')
        ax.figure.savefig(figureSavePath, bbox_inches='tight', dpi = 300)