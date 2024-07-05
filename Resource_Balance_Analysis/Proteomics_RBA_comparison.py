# -*- coding: utf-8 -*-
"""
Created on Mon Sep  4 10:37:44 2023

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
import matplotlib.pyplot as plt
import seaborn as sns
import scipy as sp

#%% set parameters
printComparisons = False
plotAll = False
logPlot = False

#%% function definitions
def findOccurrences(s, ch):
    return [i for i, letter in enumerate(s) if letter == ch]


def flatten(l):
    return [item for sublist in l for item in sublist]


#%% load in data
curFolder = os.getcwd()
backSlash = findOccurrences(curFolder, os.sep)
baseFolder = curFolder[:backSlash[-1]]

fileLocation = os.path.join(baseFolder, 'data', 'proteomics',
                            '20230620', 'data_analysis', 'model_comparisons')

# uniprot data from RBA output
fileNames = ['UniProt_GeneName_output_idmapping_2023_09_04.tsv',
             'UniProt_UniProt_ID_output_idmapping_2023_09_04.tsv']
searchDBs = ['UniProt_ID', 'GeneName']
uniprotData = {}
allUniProtIds = []

for db in searchDBs:
    for fN in fileNames:
        if db in fN:
            curFile = fN
            continue
    fullFileName = os.path.join(fileLocation, curFile)
    uniprotData[db] = pd.read_csv(fullFileName, sep="\t")
    allUniProtIds.append(list(uniprotData[db]['Entry']))

# unpack the two listsin the UniProtIds variable
allUniProtIds = flatten(allUniProtIds) # 1445 items
# check and remove duplicates
allUniProtIds = list(set(allUniProtIds)) # 1443 items


# load in proteomics data
protFileName = 'mean_allImputedData_mmol_p_g_CDW_BCA.csv'
protData = pd.read_csv(os.path.join(fileLocation, protFileName), sep="\t")
protDataIDs = list(protData['Protein IDs'])
# make sure to unpack proteins with multiple IDs
protDataIDs = [item.split(';') for item in protDataIDs]
protDataIDs = flatten(protDataIDs) # 2940 items
# check for duplicates via set
protDataIDs = list(set(protDataIDs)) # 2940 items

#%% compare which proteins are only in the proteomics data, only in the model or in both
protInBoth = [ID for ID in protDataIDs if ID in allUniProtIds] # 900 items
proteomicsOnly = [ID for ID in protDataIDs if ID not in allUniProtIds] # 2040 items
modelOnly = [ID for ID in allUniProtIds if ID not in protDataIDs] # 543 items
if printComparisons:
    print('\nProtein comparison between model and proteomics data')
    print('complete model')
    print('proteins in both: {}'.format(len(protInBoth)))
    print('proteins in proteomics only: {}'.format(len(proteomicsOnly)))
    print('proteins in model only: {}'.format(len(modelOnly)))

#%% compare proteomics data with actual flux carrying machinery and enzymes
subfolders = ['RBA_N_full_simulation', 'RBA_N_lim_simulation', 'RBA_N_lim_moreC_simulation']
modelResults = {}
# load in data
for folder in subfolders:
    modelResults[folder] = {}
    extPath = os.path.join(fileLocation, folder)
    allFiles = os.listdir(extPath)
    for file in allFiles:
        dictKey = file[:-4]
        if not 'mapping' in file:
            modelResults[folder][dictKey] = pd.read_csv(os.path.join(extPath, file), sep="\t", 
                                                        header = None, names = ['ID', 'value'])
        else:
            modelResults[folder][dictKey] = pd.read_csv(os.path.join(extPath, file), sep="\t")
    # map identifiers in the simulation results to the actual UniProt IDs
    modelResults[folder]['IDmapping'] = {}
    # get the mapping pairs from the uniprot search files from the 'From' and 'Entry' columns
    for key in modelResults[folder]:
        if 'UniProt_ID_mapping' in key:
            # convert two columns to dictionary
            mappingDict = modelResults[folder][key].set_index('From')['Entry'].to_dict()
            # merge the dictionary to the existing one
            modelResults[folder]['IDmapping'] = modelResults[folder]['IDmapping'] | mappingDict
    # add column to modelling results with the proper UniProt IDs
    # might not catch if a gene name is associated with more than one UniProt ID
    for key in modelResults[folder]:
        if 'mapping' not in key:
            origIDs = list(modelResults[folder][key]['ID'])
            newIDs = []
            for ID in origIDs:
                try:
                    newIDs.append(modelResults[folder]['IDmapping'][ID])
                except:
                    newIDs.append('ignore')
            modelResults[folder][key]['UniProtIDs'] = newIDs

# check which proteins are in the proteomics data
for folder in subfolders:
    for key in modelResults[folder]:
        if 'mapping' not in key:
            modelledIDs = list(modelResults[folder][key]['UniProtIDs'])
            modelledIDs = [i for i in modelledIDs if 'ignore' not in i]
            
            protInBoth = [ID for ID in protDataIDs if ID in modelledIDs]            # full:  377; lim  382
            proteomicsOnly = [ID for ID in protDataIDs if ID not in modelledIDs]    # full: 2563; lim 2558
            modelOnly = [ID for ID in modelledIDs if ID not in protDataIDs]         # full:   57; lim   60
            if printComparisons:
                print('\nProtein comparison between model and proteomics data')
                print(folder)
                print('proteins in both: {}'.format(len(protInBoth)))
                print('proteins in proteomics only: {}'.format(len(proteomicsOnly)))
                print('proteins in model only: {}'.format(len(modelOnly)))

#%% integrate isoenzyme information into proteomics data
# sum up isoenzyme concentration from MS data
# get those from Michael's R-notebook github
# import gpr information
if 'not_BCA' in protFileName:
    protDataGRRfileName = 'mean_not_BCA_imputed_data_with_GRR_sums.csv'
    BCAFNextension = 'not_BCA_'
else:
    protDataGRRfileName = 'mean_BCA_imputed_data_with_GRR_sums.csv'
    
    BCAFNextension = 'BCA_'
GRRprotDatafile = os.path.join(fileLocation, protDataGRRfileName)

if not os.path.exists(GRRprotDatafile):
    gprDfUniProtSaveFile = os.path.join(fileLocation, 'GPR_UniProt_IDs.csv')
    if not os.path.exists(gprDfUniProtSaveFile):
        gprFile = os.path.join(baseFolder, 'Modelling',
                               'R-notebook-ralstonia-proteome', 'data', 'input',
                               'model_reactions.csv')
        gprDf = pd.read_csv(gprFile, index_col = 0)
        geneNamesSaveFile = os.path.join(fileLocation, 'GPR_geneNames.csv')
        if not os.path.exists(geneNamesSaveFile):
            # get list of gene names to search on uniprot for identifiers
            geneNames = list(gprDf['genes'])
            geneNames = [i.split(', ') for i in geneNames if str(i) != 'nan']
            geneNames = list(set(flatten(geneNames))) # 1344 items
            geneNames = [i + '\n' for i in geneNames]
            with open(geneNamesSaveFile, 'w') as gprFile:
                gprFile.writelines(geneNames)
        
        gprUniProtIDFile = os.path.join(fileLocation, 'GPR_gene_IDs_idmapping_2023_09_04.tsv')
        gprUniProtIDs = pd.read_csv(gprUniProtIDFile, sep='\t')
        gprMappingDict = gprUniProtIDs.set_index('From')['Entry'].to_dict()
        gprDfUniProt = gprDf.copy()
        for key in gprMappingDict:
            gprDfUniProt = gprDfUniProt.replace(key, gprMappingDict[key], regex=True)
        gprDfUniProt.to_csv(gprDfUniProtSaveFile, sep = '\t', index = False)
    else:
        gprDfUniProt = pd.read_csv(gprDfUniProtSaveFile, sep = '\t')
    
    # get list of all gene-reaction-relations with more than one gene (protein) associated with it
    allGprs = list(gprDfUniProt['genes'])
    allGRRs = [gpr.split(', ') for gpr in allGprs if str(gpr) != 'nan'] # 1051 items
    allGRRs = [gpr for gpr in allGRRs if len(gpr) > 1] # 484 items
    
    # add all rows in the original dataframe that contain IDs belonging to the same
    # group in the allGRRs and remove the original rows
    originalProtData = protData.copy()
    allRowIndices = []
    for grr in allGRRs:
        rowIndex = []
        for protein in grr:
            try:
                # should error if the protein ID is not found (i.e. the protein was not detected)
                protIndex = protData.loc[protData['Protein IDs'].str.contains(protein)].index[0]
                rowIndex.append(protIndex)
            except:
                continue
        if len(rowIndex) > 1:
            # build the sum of all the rows
            # build an ID string
            grrDf = protData.loc[rowIndex]
            grrProtIDs = list(grrDf['Protein IDs'])
            combinedIDs = ';'.join(grrProtIDs)
            IDseries = pd.Series(data = {'Protein IDs': combinedIDs})
            # sum all the numeric values
            rowSum = grrDf.sum(numeric_only = True)
            # join the two series together
            newRow = pd.concat([IDseries, rowSum])
            newRow = newRow.replace(0, np.NaN)
            newRow = pd.DataFrame(newRow).transpose()
            newRow = newRow.infer_objects()
            # append the dataframe
            protData = pd.concat([protData, newRow], ignore_index = True)
            # protData = protData.drop(rowIndex)
            allRowIndices.append(rowIndex)
    
    for rowIndices in allRowIndices:
        for index in rowIndices:
            try:
                protData = protData.drop(index)
            except:
                # assume the row has been dropped already
                continue
    
    protData.to_csv(GRRprotDatafile, sep='\t', index = False)
else:
    protData = pd.read_csv(GRRprotDatafile, sep='\t')

#%% compare/correllate proteomics data with modelling results
# focus on result files with only active enzymes
# find groups of columns and split proteomics dataframe into condition specific dataframes
allColumns = list(protData.columns)
allConditions = [i for i in allColumns if 'Protein' not in i]
trimmedConditions = []
for condition in allConditions:
    undPositions = findOccurrences(condition, '_')
    trimmedConditions.append(condition[:undPositions[2]])
uniqueConditions = list(set(trimmedConditions))

exprByCond = {}
for cond in uniqueConditions:
    selectedCols = [i for i in allColumns if cond in i or 'Protein' in i]
    exprByCond[cond] = protData[selectedCols]
    exprByCond[cond] = exprByCond[cond].rename(columns = {'Protein IDs':'UniProtIDs'})
    # expand the dataframe by a column with the modelling results
    if 'full' in cond:
        modelKey = 'RBA_N_full_simulation'
    else:
        if '150' in cond:
            modelKey = 'RBA_N_lim_moreC_simulation'
        else:
            modelKey = 'RBA_N_lim_simulation'
    for key in modelResults[modelKey]:
        if 'proteins' in key:
            subKey = key
    modelData = modelResults[modelKey][subKey]
    reducedModelData = modelData[modelData['UniProtIDs'].str.contains('ignore') == False]
    reducedModelData = reducedModelData[['UniProtIDs', 'value']]
    reducedModelData.rename(columns = {'value':'predicted_value'}, inplace = True)
    mappingDict = reducedModelData.set_index('UniProtIDs')['predicted_value'].to_dict()
    proteomicsIDs = list(exprByCond[cond]['UniProtIDs'])
    # create list of nans to write into from the mapping dict
    predictionResults = np.empty((len(proteomicsIDs)))
    predictionResults[:] = np.nan
    predictionResults = list(predictionResults)
    mappingDictOrig = mappingDict.copy()
    for key in mappingDictOrig:
        indices = [i for i, s in enumerate(proteomicsIDs) if key in s]
        for idx in indices:
            # make sure that if multiple IDs are mapped all of the results are captured
            # --> now the mean for that is calculated whilst ignoring nan (np.nanmean())
            if isinstance(predictionResults[idx], list):
                predictionResults[idx] = predictionResults[idx].append(mappingDict[key])
            else:
                predictionResults[idx] = [mappingDict[key]]
        if len(indices) != 0:
            mappingDict.pop(key)
    for i in range(len(predictionResults)):
        if isinstance(predictionResults[i], list):
            predictionResults[i] = np.nanmean(predictionResults[i])
    exprByCond[cond]['predicted_value'] = predictionResults
    remainingPredictions = pd.DataFrame.from_dict(mappingDict, orient='index',
                                                  columns = ['predicted_value'])
    remainingPredictions = remainingPredictions.reset_index()
    remainingPredictions = remainingPredictions.rename(columns = {'index':'UniProtIDs'})
    exprByCond[cond] = pd.concat([exprByCond[cond], remainingPredictions], join = 'outer')

#%% plot the results
# 1) only plot data that are present in the proteomics and the predicted set
# 2) plot all data, replacing nan with 0 in the condition dependent dataframes

# get the order of plots
allConditions = list(exprByCond.keys())
allConditions = [i.replace('80','080') for i in allConditions]
allConditions.sort()
allConditions = [i.replace('080','80') for i in allConditions]
# create the plotting axes
nRows = 2
nCols = 3
plt.rcParams["font.family"] = "Arial"
fig, axes = plt.subplots(nRows, nCols, sharex=True, sharey=True, figsize=(7,5))

axisNumber = 0
titleDict = {'ALE26_80_full':   'ALE26, 80 C-mM, C/N = 5.3',
             'ALE26_80_lim':    'ALE26, 80 C-mM, C/N = 50',
             'ALE26_150_lim':   'ALE26, 150 C-mM, C/N = 50',
             'H16_80_full':     'H16, 80 C-mM, C/N = 5.3',
             'H16_80_lim':      'H16, 80 C-mM, C/N = 50',
             'H16_150_lim':     'H16, 150 C-mM, C/N = 50'}

# plot per condition
for cond in allConditions:
    plottingDf = exprByCond[cond]
    # multiply all values by 1000 to get it to umol/L
    plottingDf[plottingDf.select_dtypes(include=['number']).columns] *= 1000
    # rename columns for plotting
    pltColsOrig = list(plottingDf.columns)
    renameDict = {}
    for col in pltColsOrig:
        if '_mean_' in col:
            renameDict[col] = 'proteomics c_mean [umol/g(CDW)]'
        elif '_std_' in col:
            renameDict[col] = 'proteomics c_std [umol/g(CDW)]'
        elif 'predicted' in col:
            renameDict[col] = 'RBA c_predicted [umol/g(CDW)]'
    
    plottingDf = plottingDf.rename(columns = renameDict)
    if logPlot:
        plottingDf['log(c_measured [umol/g(CDW)])'] = np.log(plottingDf['proteomics c_mean [umol/g(CDW)]'])
        plottingDf['log(c_measured (std) [umol/g(CDW)])'] = np.log(plottingDf['proteomics c_std [umol/g(CDW)]'])
        plottingDf['log(c_predicted [umol/g(CDW)])'] = np.log(plottingDf['RBA c_predicted [umol/g(CDW)]'])
        xDataName = 'log(c_predicted [umol/g(CDW)])'
        yDataName = 'log(c_measured [umol/g(CDW)])'
        logFNextension = 'log_'
        rSquaredY = 0.05
    else:
        plottingDf['c_measured [umol/g(CDW)]'] = plottingDf['proteomics c_mean [umol/g(CDW)]']
        plottingDf['c_measured (std) [umol/g(CDW)]'] = plottingDf['proteomics c_std [umol/g(CDW)]']
        plottingDf['c_predicted [umol/g(CDW)]'] = plottingDf['RBA c_predicted [umol/g(CDW)]']
        xDataName = 'c_predicted [umol/g(CDW)]'
        yDataName = 'c_measured [umol/g(CDW)]'
        logFNextension = ''
        rSquaredY = 0.8
    
    if plotAll:
        #replace the nan values
        plottingDf = plottingDf.fillna(0)
        nanFNextension = 'withNan_'
    else:
        nanFNextension = 'withoutNan_'
    
    # save data to file for sharing with Victor and Marci
    # PHDfolder = curFolder[:backSlash[2]]
    # dataFileLocation = os.path.join(PHDfolder, 'Publications', 'ALE_paper',
    #                                 'Manuel_data_for_origin', 'figure6_b')
    # resDataFileName = cond + '_measured_v_predicted_data.csv'
    # plottingDf.to_csv(os.path.join(dataFileLocation, resDataFileName), index = False)
    
    # plot data
    plt.sca(fig.axes[axisNumber])
    rp = sns.regplot(data = plottingDf,
                     x = xDataName,
                     y = yDataName,
                     scatter_kws = {"color": "black", "alpha": 0.5},
                     line_kws = {"color": "red"})
    rp.set(ylim=(0, 0.5))
    rp.set(xlim=(0, 0.5))
    x_lim = rp.get_xlim()
    y_lim = rp.get_ylim()
    plt.xticks(np.arange(min(x_lim), max(x_lim)+0.1, 0.25))
    plt.yticks(np.arange(min(y_lim), max(y_lim)+0.1, 0.25))
    font2 = {'family': 'arial', 'color': 'black', 'size': 10, 'weight': 'bold'} # 
    fig.axes[axisNumber].set_title(titleDict[cond], fontdict = font2)
    fig.axes[axisNumber].set(xlabel = None)
    fig.axes[axisNumber].set(ylabel = None)
    fig.axes[axisNumber].ticklabel_format(axis = 'x', style = 'plain')#, scilimits = (0, 0))
    fig.axes[axisNumber].ticklabel_format(axis = 'y', style = 'plain')#, scilimits = (0, 0))
    plt.tick_params(axis = 'both', which = 'major', labelsize = 9)
    # compute correlation statistics of potentially transformed data
    statsDf = plottingDf.dropna()
    curAx = plt.gca()
    r, p = sp.stats.pearsonr(statsDf[xDataName], statsDf[yDataName])
    # display R^2 and p value of correllation
    # double {{}} to escape brackets used in the .format() argument
    curAx.text(.1, rSquaredY, '$\mathregular{{R^{{2}}}}$={:.2f}, p={:.2g}'.format(r, p),
               transform = curAx.transAxes, family = 'arial', color = 'black', size = 9)
    
    axisNumber += 1
    
if logPlot:
    xDataName = 'log($\mathregular{c_{predicted}}$ ($\mu$mol/g(CDW)))'
    yDataName = 'log($\mathregular{c_{measured}}$ ($\mu$mol/g(CDW)))'
else:
    xDataName = '$\mathregular{c_{predicted}}$ ($\mu$mol/g(CDW))'
    yDataName = '$\mathregular{c_{measured}}$ ($\mu$mol/g(CDW))'
fig.text(0.5, 0.02, xDataName, ha ='center', family = 'arial', color = 'black', size = 10)
fig.text(0.04, 0.5, yDataName, va ='center', rotation='vertical', family = 'arial', color = 'black', size = 10)

# export the plot
plotFileName = BCAFNextension + logFNextension + nanFNextension + 'RBA_proteomics_comparison_thesis'
fullPlotSavePath = os.path.join(fileLocation, 'figures', plotFileName + '.svg')
fig.figure.savefig(fullPlotSavePath, format='svg')

fullPlotSavePath = os.path.join(fileLocation, 'figures', plotFileName + '.png')
fig.figure.savefig(fullPlotSavePath, bbox_inches='tight', dpi = 1000)