# -*- coding: utf-8 -*-
"""
Created on Wed Nov  1 12:06:45 2023

@author: Manuel Bruch
"""

#%% clean up
from IPython import get_ipython
get_ipython().magic('reset -sf') # clears all variables

#%% import statements
import os
import pickle
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy as sp


#%% set parameters
plotAll = False
logPlot = False
meanOrMedian = 'mean'
muMax = '' # '' or '_muMax'
setMu = 'mu_0_09_' # '' or 'mu_0_09_'
setDate = '240115' # '240109' or '240115'

#%% function definitions
def findOccurrences(s, ch):
    return [i for i, letter in enumerate(s) if letter == ch]


#%% load in data
# load in dictionaries with RBVA data
curFolder = os.getcwd()
backSlash = findOccurrences(curFolder, '\\')
RBVAbaseFolder = curFolder[:backSlash[-1]]
fileList = [setDate + '_RBVA_' + setMu + '80mMFA_fullN'  + muMax + '_ProteinData.csv',
            setDate + '_RBVA_' + setMu + '80mMFA_limN'   + muMax + '_ProteinData.csv',
            setDate + '_RBVA_' + setMu + '150mMFA_limN'  + muMax + '_ProteinData.csv']
allRBVAdata = {}

for file in fileList:
    fullPath = os.path.join(RBVAbaseFolder, 'RBVA_MJahn_model',
                            'manual_RBVA_results', file)
    underScore = findOccurrences(file, '_')
    # get the RBVA conditions for the dict keys
    dictKey = file[underScore[1] + 1 : underScore[-1]]
    allRBVAdata[dictKey] = {}
    statFile = dictKey + 'meanMaxMinData.csv'
    statFileFullPath = os.path.join(RBVAbaseFolder, 'RBVA_MJahn_model',
                                    'manual_RBVA_results', statFile)
    if not os.path.exists(statFileFullPath):
        allRBVAdata[dictKey]['fullData'] = pd.read_csv(fullPath, index_col = 0)
        allRBVAdata[dictKey]['meanMaxMinData'] = pd.DataFrame()
        allRBVAdata[dictKey]['meanMaxMinData']['mean'] = allRBVAdata[dictKey]['fullData'].mean(axis = 1)
        allRBVAdata[dictKey]['meanMaxMinData']['std'] = allRBVAdata[dictKey]['fullData'].std(axis = 1)
        allRBVAdata[dictKey]['meanMaxMinData']['median'] = allRBVAdata[dictKey]['fullData'].median(axis = 1)
        allRBVAdata[dictKey]['meanMaxMinData']['min'] = allRBVAdata[dictKey]['fullData'].min(axis = 1)
        allRBVAdata[dictKey]['meanMaxMinData']['max'] = allRBVAdata[dictKey]['fullData'].max(axis = 1)
        allRBVAdata[dictKey]['meanMaxMinData'].to_csv(statFileFullPath, 
                                                      sep='\t', index = True)
    else:
        allRBVAdata[dictKey]['meanMaxMinData'] = pd.read_csv(statFileFullPath,
                                                             sep='\t',
                                                             index_col = 0)


# load in a dictionary that maps RBA reaction names with gene identifiers
fullPath  = os.path.join(RBVAbaseFolder, 'C_necator_MJahn_RBA_model_proteins.pkl')
with open(fullPath, 'rb') as f:
    RBAproteins = pickle.load(f)
    
# get all the keys to search in UniProt
searchListPath = os.path.join(curFolder, 'RBAproteinIDs.txt')
if not os.path.isfile(searchListPath):
    RBAproteinIDs = list(RBAproteins.keys())
    with open(searchListPath, mode = 'x', encoding = 'utf-8') as f:
        f.write('\n'.join(RBAproteinIDs))

# load in data from uniprot search
uniprotIDfiles = ['idmapping_2023_11_01.tsv', 'idmapping_2023_11_01_run02.tsv']
allRBVAidMaps = {}
for f in uniprotIDfiles:
    uniprotIDs = pd.read_csv(f, sep="\t")
    allRBVAidMaps = allRBVAidMaps | uniprotIDs.set_index('From')['Entry'].to_dict()


# load in proteomics data
baseFolder = curFolder[:backSlash[-3]]
fileLocation = os.path.join(baseFolder, 'data', 'proteomics',
                            '20230620', 'data_analysis', 'model_comparisons')
protDataGRRfileName = 'mean_not_BCA_imputed_data_with_GRR_sums.csv'
GRRprotDatafile = os.path.join(fileLocation, protDataGRRfileName)
protData = pd.read_csv(GRRprotDatafile, sep='\t')


#%% add column with uniprot IDs to RBVA results
RBVAresults = {}
for key in allRBVAdata:
    RBVAresults[key] = allRBVAdata[key]['meanMaxMinData']
    RBVAresults[key]['UniProtIDs'] = RBVAresults[key].index.map(allRBVAidMaps)

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
    if 'full' in cond:
        modelKey = setMu + '80mMFA_fullN' + muMax
    else:
        if '150' in cond:
            modelKey = setMu + '150mMFA_limN' + muMax
        else:
            modelKey = setMu + '80mMFA_limN' + muMax
    modelData = RBVAresults[modelKey]
    reducedModelData = modelData[modelData['UniProtIDs'].str.contains('nan') == False]
    reducedModelData = reducedModelData[['UniProtIDs', meanOrMedian]]
    reducedModelData.rename(columns = {meanOrMedian:'predicted_value'}, inplace = True)
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
fig, axes = plt.subplots(nRows, nCols, sharex=True, sharey=True, figsize=(7,4))

axisNumber = 0
titleDict = {'ALE26_80_full':   'ALE26, 80 C-mM',
             'ALE26_80_lim':    'ALE26, 80 C-mM, N-lim',
             'ALE26_150_lim':   'ALE26, 150 C-mM, N-lim',
             'H16_80_full':     'H16, 80 C-mM',
             'H16_80_lim':      'H16, 80 C-mM, N-lim',
             'H16_150_lim':     'H16, 150 C-mM, N-lim'}

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
    
    # plot data
    plt.sca(fig.axes[axisNumber])
    rp = sns.regplot(data = plottingDf,
                     x = xDataName,
                     y = yDataName,
                     scatter_kws = {"color": "black", "alpha": 0.5},
                     line_kws = {"color": "red"})
    fig.axes[axisNumber].set_title(titleDict[cond], fontsize = 9)
    fig.axes[axisNumber].set(xlabel = None)
    fig.axes[axisNumber].set(ylabel = None)
    fig.axes[axisNumber].ticklabel_format(axis = 'x', style = 'sci', scilimits = (0, 0))
    fig.axes[axisNumber].ticklabel_format(axis = 'y', style = 'sci', scilimits = (0, 0))
    plt.tick_params(axis = 'both', which = 'major', labelsize = 9)
    # compute correlation statistics of potentially transformed data
    statsDf = plottingDf.dropna()
    curAx = plt.gca()
    r, p = sp.stats.pearsonr(statsDf[xDataName], statsDf[yDataName])
    # display R^2 and p value of correllation
    # double {{}} to escape brackets used in the .format() argument
    curAx.text(.2, rSquaredY, '$\mathregular{{R^{{2}}}}$={:.2f}, p={:.2g}'.format(r, p),
               transform = curAx.transAxes, family = 'arial', color = 'black', size = 9)
    
    axisNumber += 1
fig.text(0.5, 0.004, xDataName, ha ='center', family = 'arial', color = 'black', size = 9)
fig.text(0.05, 0.5, yDataName, va ='center', rotation='vertical', family = 'arial', color = 'black', size = 9)

# export the plot
if len(muMax) > 0:
    plotFileName = logFNextension + nanFNextension + 'RBA_proteomics_comparison_word_rbatools' + 'mu_0_09'
else:
    plotFileName = logFNextension + nanFNextension + 'RBVA_proteomics_comparison_word_rbatools' + 'mu_0_09'

# fullPlotSavePath = os.path.join(fileLocation, 'figures', 'RBVA_results',
#                                 plotFileName + '.svg')
# fig.figure.savefig(fullPlotSavePath, format='svg')

# fullPlotSavePath = os.path.join(fileLocation, 'figures', 'RBVA_results',
#                                 plotFileName + '.png')
# fig.figure.savefig(fullPlotSavePath, bbox_inches='tight', dpi = 500)