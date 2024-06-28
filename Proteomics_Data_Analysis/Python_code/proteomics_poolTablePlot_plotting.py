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

#%% get individual pathways that are now allocated
tableHeaders = list(dataTable.columns)
metaPWHeader = [i for i in tableHeaders if 'Meta_Pathways' in i][0]
pathways = list(dataTable[metaPWHeader]) # 2885 items, 1709 of which are 'not annotated' --> 1176 annotated. only 354 in UniProt, 1150 in KEGG
pathways_split = [item.split('//') for item in pathways]
pathways_flat = [item for l in pathways_split for item in l] # 3695 items
pathways_set = list(set(pathways_flat)) # 29 items

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
                           'Protein families: signaling and cellular processes':    'Signaling and cellular processes',
                           'Information processing in viruses':                     'Viral information processing',
                           'not annotated':                                         'Not annotated',
                           'Biosynthesis of other secondary metabolites':           'Secondary metabolite synthesis',
                           'Metabolism of cofactors and vitamins':                  'Cofactor/Vitamin metab.',
                           'Cell motility':                                         'Cell motility',
                           'Degradation of aromatic compounds':                     'Aromatics degradation',
                           'Xenobiotics biodegradation and metabolism':             'Xenobiotics metab.',
                           'Protein families: genetic information processing':      'Genetic information processing',
                           'Amino acid metabolism':                                 'Amino acid metab.',
                           'Drug resistance: antimicrobial':                        'antimicrobial drug resistance',
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
conditionHeaders = [i for i in tableHeaders if 'log(' in i]
allConditions = [i.replace(' log(fold change)','') for i in conditionHeaders]
allConditions = [i.replace(' -log(adjusted p)','') for i in allConditions]
allConditions = list(set(allConditions))

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
        if val > pCutoff:
            pCutoffList.append('significant')
        else:
            pCutoffList.append('not significant')
    conditionDFDict[curKey]['significance'] = pCutoffList

#%% plot some stuff
# get the global min and max values
conditionDF = appendedDF[conditionHeaders]
globalMin = min(list(conditionDF.min()))
globalMax = max(list(conditionDF.max()))

categoryOrder = list(pathwayNameReallocation.values())
categoryOrder.sort()
for curKey in conditionDFDict:
    curPlot = sns.catplot(data = conditionDFDict[curKey], x = "log(fold change)",
                          y = 'Meta pathways', hue = "significance",
                          kind = "strip", size = 3, order = categoryOrder, 
                          height = 6, aspect = 1.5).set(title = curKey + '\n')
    curPlot.set(xlim=(globalMin, globalMax))
    figureSavePath = os.path.join(csvFileLocation, 'figures', curKey + '.svg')
    # curPlot.figure.savefig(figureSavePath, format='svg')
    figureSavePath = os.path.join(csvFileLocation, 'figures', curKey + '.png')
    # curPlot.figure.savefig(figureSavePath, bbox_inches='tight', dpi = 500)