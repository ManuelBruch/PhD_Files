# -*- coding: utf-8 -*-
"""
Created on Tue Aug 29 12:41:14 2023

@author: Manuel Bruch
relative filepaths corrected for all OS
"""

#%% clean up
from IPython import get_ipython
get_ipython().magic('reset -sf') # clears all variables

#%% import statements
import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

#%% function definitions
def findOccurrences(s, ch):
    return [i for i, letter in enumerate(s) if letter == ch]


#%% read data frame with pathway and expression data
curFolder = os.getcwd()
backSlash = findOccurrences(curFolder, os.sep)
baseFolder = curFolder[:backSlash[-1]]

csvFileLocation = os.path.join(baseFolder, 'data', 'proteomics', 
                            '20230620', 'data_analysis', 'poolTablePlots')

dataFile = 'poolTable_df_merged_metaPWs.csv'

dataTable = pd.read_csv(os.path.join(csvFileLocation, dataFile), sep="\t")
dataTable = dataTable.sort_values('UniProtIDs', ascending = True)
dataTableOriginal = dataTable.copy() # to check if the artificial assigning of values works as intended

#%% get individual pathways that are now allocated
tableHeaders = list(dataTable.columns)
metaPWHeader = [i for i in tableHeaders if 'Meta_Pathways' in i][0]
pathways = list(dataTable[metaPWHeader]) # 2885 items, 1709 of which are 'not annotated' --> 1176 annotated. only 354 in UniProt, 1150 in KEGG
pathways_split = [item.split('//') for item in pathways]
pathways_flat = [item for l in pathways_split for item in l] # 3695 items
pathways_set = list(set(pathways_flat)) # 29 items

#%% replace nan values in log(fold change) column depending on if and in which strain that protein would have been detected
# load in the protein count dataframe
protCountLocation = os.path.join(baseFolder, 'data', 'proteomics', 
                                 '20230620', 'data_analysis', 'proteomics_condition_count.csv')
protCountDf = pd.read_csv(protCountLocation, sep="\t")
countHeaders = list(protCountDf.columns)
countHeaders = [i for i in countHeaders if 'count' in i]
protCountDf = protCountDf[['proteinIDs'] + countHeaders]
# change the headers to match with dataTable
colNameChange = {'proteinIDs': 'UniProtIDs',
                 'H16 full 80 count': 'H16_80_full',
                 'H16 CN50 80 count': 'H16_80_lim',
                 'H16 CN50 150 count': 'H16_150_lim',
                 'ALE26 full 80 count': 'ALE26_80_full',
                 'ALE26 CN50 80 count': 'ALE26_80_lim',
                 'ALE26 CN50 150 count': 'ALE26_150_lim'}
protCountDf = protCountDf.rename(columns = colNameChange)
protCountDf = protCountDf.sort_values('UniProtIDs', ascending = True)

# get groupong variables for the conditions
conditionHeaders = [i for i in tableHeaders if 'log(' in i]
allConditions = [i.replace(' log(fold change)','') for i in conditionHeaders]
allConditions = [i.replace(' -log(adjusted p)','') for i in allConditions]
allConditions = list(set(allConditions))

for condition in allConditions:
    if 'H16' in condition and 'ALE26' in condition:
        # print(condition)
        foldChange = list(dataTable[condition +  ' log(fold change)'])
        pValue = list(dataTable[condition +  ' -log(adjusted p)'])
        # nanCount = np.count_nonzero(np.isnan(foldChange))
        comparedCond = condition.split('-')
        # get the condition dependent counts relevant in the comparison
        for cond in comparedCond:
            if 'H16' in cond:
                H16count = list(protCountDf[cond])
            else:
                ALE26count = list(protCountDf[cond])
        increaseH16 = 0
        increaseALE26 = 0
        for i in range(len(foldChange)):
            if np.isnan(foldChange[i]):
                if H16count[i] > 1:
                    if condition.startswith('H16'):
                        foldChange[i] = 15 + increaseH16
                        pValue[i] = 16
                    else:
                        foldChange[i] = -15 - increaseH16
                        pValue[i] = 16
                    increaseH16 += 0.01
                elif ALE26count[i] > 1:
                    if condition.startswith('ALE26'):
                        foldChange[i] = 15 + increaseALE26
                        pValue[i] = 26
                    else:
                        foldChange[i] = -15 - increaseALE26
                        pValue[i] = 26
                    increaseALE26 += 0.01
        dataTable[condition +  ' log(fold change)'] = foldChange
        dataTable[condition +  ' -log(adjusted p)'] = pValue

#%% split rows with multiple pathways in the Meta_Pathways column
# convert dataframe into dictionary of rows
dataDict = dataTable.to_dict(orient='index')
originalIndexes = list(dataDict.keys())
# iterate over dictionary and populate a second one with the duplicates (can't expand the dictionary one is iterating over)
nextKey = originalIndexes[-1] + 1
dataDict2 = {}
for key in dataDict:
    if '//' in dataDict[key][metaPWHeader]:
        dataDict[key][metaPWHeader] = dataDict[key][metaPWHeader].split('//')
        while len(dataDict[key][metaPWHeader]) > 1:
            poppedPw = dataDict[key][metaPWHeader].pop(-1)
            dataDict2[nextKey] = dataDict[key].copy()
            dataDict2[nextKey][metaPWHeader] = poppedPw
            nextKey += 1
        dataDict[key][metaPWHeader] = dataDict[key][metaPWHeader][0]

merged_dict = dataDict | dataDict2 # merge two dictionaries

# create a new dataframe from the dictionary
appendedDF = pd.DataFrame(merged_dict).transpose()
appendedDF = appendedDF.sort_values(['UniProtIDs', metaPWHeader],
              ascending = [True, True])
appendedDF.reset_index(inplace = True, drop = True)
appendedDF = appendedDF.infer_objects() # automatically convert float columns to float (and would convert others to better objects as well)

#%% alter pathway names for plotting
pathwayNameReallocation = {'Folding, sorting and degradation':                      'Folding, sorting, degradation', 
                           'Unclassified: metabolism':                              'Unclassified',
                           'Membrane transport':                                    'Membrane transport',
                           'Translation':                                           'Translation',
                           'Protein families: signaling and cellular processes':    'Signaling & cellular procs.',
                           'Information processing in viruses':                     'Viral information processing',
                           'not annotated':                                         'Not annotated',
                           'Biosynthesis of other secondary metabolites':           'Secondary metabolite synth.',
                           'Metabolism of cofactors and vitamins':                  'Cofactor/Vitamin metab.',
                           'Cell motility':                                         'Cell motility',
                           'Degradation of aromatic compounds':                     'Aromatics degradation',
                           'Xenobiotics biodegradation and metabolism':             'Xenobiotics metab.',
                           'Protein families: genetic information processing':      'Genetic information proc.',
                           'Amino acid metabolism':                                 'Amino acid metab.',
                           'Drug resistance: antimicrobial':                        'Antimicrobial drug resistance',
                           'Carbohydrate metabolism':                               'Carbohydrate metab.',
                           'Biosynthesis of cofactors':                             'Biosynthesis of cofactors',
                           'Carbon metabolism':                                     'Carbon metab.',
                           'Energy metabolism':                                     'Energy metab.',
                           'Cell growth and death':                                 'Cell growth and death',
                           'Metabolism of other amino acids':                       'Other amino acids',
                           'Transcription':                                         'Transcription',
                           'Nucleotide metabolism':                                 'Nucleotide metab.',
                           'Glycan biosynthesis and metabolism':                    'Glycan synthesis & metab.',
                           'Fatty acid metabolism':                                 'Fatty acid metab.',
                           'Cellular community - prokaryotes':                      'Quorum sensing',
                           'Lipid metabolism':                                      'Lipid metab.',
                           'Replication and repair':                                'Replication and DNA repair',
                           'Metabolism of terpenoids and polyketides':              'Terpenoid/polyketid metab.'}

newPWList = list(appendedDF[metaPWHeader])
newPWList = [pathwayNameReallocation[item] for item in newPWList]
appendedDF['Meta pathways'] = newPWList
tableHeaders = list(appendedDF.columns)

#%% split the appended dataframe into dataframes per condition in a dictionary
headersToKeep = [i for i in tableHeaders if not 'log(' in i]

conditionDFDict = {}
pCutoff = -np.log(0.05)

for condition in allConditions:
    curKey = condition.replace('-',' vs ')
    conditionDFDict[curKey] = appendedDF[headersToKeep].copy()
    statCols = [i for i in tableHeaders if condition in i]
    colNames = [i.replace(condition + ' ','') for i in statCols]
    for i in range(len(colNames)):
        conditionDFDict[curKey][colNames[i]] = appendedDF[statCols[i]]
    # add column to distinguish statistically significant and not significant values
    pValueCol = [i for i in colNames if 'adjusted p' in i][0]
    pValues = list(conditionDFDict[curKey][pValueCol])
    pCutoffList = []
    for val in pValues:
        if val == 16:
            pCutoffList.append('H16 only')
        elif val == 26:
            pCutoffList.append('ALE26 only')
        elif val > pCutoff:
            pCutoffList.append('significant')
        else:
            pCutoffList.append('not significant')
    conditionDFDict[curKey]['significance'] = pCutoffList

#%% plot some stuff
# get the global min and max values
conditionDF = appendedDF[conditionHeaders]
globalMin = min(list(conditionDF.min()))
# need to remove artificial p-Values from list
maxList = list(conditionDF.max())
maxList = [i for i in maxList if i != 26]
globalMax = max(maxList)

categoryOrder = list(pathwayNameReallocation.values())
categoryOrder.sort()
hueOrder = ['not significant', 'significant', 'H16 only', 'ALE26 only']

for curKey in conditionDFDict:
    curPlot = sns.catplot(data = conditionDFDict[curKey], x = "log(fold change)",
                          y = 'Meta pathways', hue = "significance",
                          kind = "strip", size = 3, order = categoryOrder, 
                          height = 5, hue_order = hueOrder, 
                          aspect = 1.5,
                          legend = False)#.set(title = curKey + '\n')
    curPlot.set(xlim=(globalMin, globalMax))
    
    # curPlot._legend.prop.set_family(['arial'])
    # curPlot._legend.prop.set_size(9)
    # plt.legend(title='Significance', fontsize = 9)
    font1 = {'family':'arial', 'color':'black', 'size':10}
    plt.xlabel('log (fold change)', fontdict = font1);
    plt.ylabel('', fontdict = font1);
    # plt.title(curKey + '\n', fontsize = 9)
    plt.tick_params(axis = 'both', which = 'major', labelsize = 9)
    
    # figureSavePath = os.path.join(csvFileLocation, 
    #                               'figures_accounting_for_nan', 'publication',
    #                               curKey + '.svg')
    # # curPlot.figure.savefig(figureSavePath, format='svg')
    # figureSavePath = os.path.join(csvFileLocation, 
    #                               'figures_accounting_for_nan', 'publication',
    #                               curKey + '.png')
    # # curPlot.figure.savefig(figureSavePath, bbox_inches='tight', dpi = 500)

#%% plot data for the ALE publication
# select conditions to plot:
# H16_80_full vs ALE26_80_full
# H16_80_lim vs ALE26_80_lim
# H16_150_lim vs ALE26_150_lim
desiredKeys = ['H16_80_full vs ALE26_80_full',
               'H16_80_lim vs ALE26_80_lim',
               'H16_150_lim vs ALE26_150_lim']
plotDict = {k: v for (k,v) in conditionDFDict.items() if k in desiredKeys}
plotTitles = ['80 mM FA, C/N = 5.3\n',
              '80 mM FA, C/N = 50\n',
              '150 mM FA, C/N = 50\n'];

# follow this website to get plots: https://datavizpyr.com/seaborn-join-two-plots-with-shared-y-axis/
f, axs = plt.subplots(1, len(plotDict), sharey = True, figsize = (9, 6)) # for A4 page: figsize = (6, 5)
for i in range(len(plotDict)):
    
    # save data to file for sharing with Victor and Marci
    PHDfolder = curFolder[:backSlash[2]]
    dataFileLocation = os.path.join(PHDfolder, 'Publications', 'ALE_paper',
                                    'Manuel_data_for_origin', 'figure6_a')
    resDataFileName = desiredKeys[i] + '.csv'
    plotDict[desiredKeys[i]].to_csv(os.path.join(dataFileLocation, resDataFileName), index = False)
    
    # plot data
    curPlot = sns.stripplot(data = plotDict[desiredKeys[i]],
                            x = "log(fold change)", y = 'Meta pathways',
                            hue = "significance", size = 3,
                            palette = 'colorblind',
                            order = categoryOrder, hue_order = hueOrder,
                            legend = False, ax = axs[i])
    curPlot.set(xlim=(globalMin, globalMax))
    if i == 0:
        sns.despine(ax = axs[i])
    else:
        sns.despine(ax = axs[i], left = True)
        axs[i].tick_params(left=False)
    
    font1 = {'family': 'arial', 'color': 'black', 'size': 10} # for A4 page: 'size':9
    font2 = {'family': 'arial', 'color': 'black', 'weight': 'bold', 'size': 12} # for A4 page: 'size':9
    if i == 1:
        axs[i].set_xlabel('log (fold change)', fontdict = font1);
    else:
        axs[i].set_xlabel('', fontdict = font1);
    axs[i].set_ylabel('', fontdict = font1);
    axs[i].set_title(plotTitles[i], fontdict = font2)
    axs[i].tick_params(axis = 'both', which = 'major', labelsize = 9) # for A4 page: labelsize = 9

# plt.ylabel('', fontdict = font1);
f.tight_layout()

# save the figure
# figureSavePath = os.path.join(csvFileLocation, 
#                               'figures_accounting_for_nan', 'poster',
#                               'combined_poolTables_ppt.svg')
# curPlot.figure.savefig(figureSavePath, format='svg')
# figureSavePath = os.path.join(csvFileLocation, 
#                               'figures_accounting_for_nan', 'poster',
#                               'combined_poolTables_ppt.png')
# curPlot.figure.savefig(figureSavePath, bbox_inches='tight', dpi = 500)