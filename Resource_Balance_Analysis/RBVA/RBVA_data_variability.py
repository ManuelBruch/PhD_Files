# -*- coding: utf-8 -*-
"""
Created on Wed Jan 10 15:24:44 2024

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
setMu = '' # '' or 'mu_0_09_'
setDate = '240109' # '240109' or '240115'


#%% function definitions
def findOccurrences(s, ch):
    return [i for i, letter in enumerate(s) if letter == ch]


#%% load in data
# load in dictionaries with RBVA data
curFolder = os.getcwd()
backSlash = findOccurrences(curFolder, '\\')
RBVAbaseFolder = curFolder[:backSlash[-1]]
fileList = [setDate + '_RBVA_' + setMu + '80mMFA_fullN_ProteinData.csv',
            setDate + '_RBVA_' + setMu + '80mMFA_limN_ProteinData.csv',
            setDate + '_RBVA_' + setMu + '150mMFA_limN_ProteinData.csv']
numerics = ['int16', 'int32', 'int64', 'float16', 'float32', 'float64']

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
    allRBVAdata[dictKey] = pd.read_csv(statFileFullPath, sep='\t',
                                        index_col = 0)
    # log transformation
    for c in [c for c in allRBVAdata[dictKey].columns if allRBVAdata[dictKey][c].dtype in numerics]:
        allRBVAdata[dictKey][c] = np.log10(allRBVAdata[dictKey][c])
    allRBVAdata[dictKey].replace(-np.inf, -50, inplace = True) # replace with large negative number to show area in plot
    allRBVAdata[dictKey]['range'] = allRBVAdata[dictKey]['max'] - allRBVAdata[dictKey]['min']
    allRBVAdata[dictKey]['range'] = allRBVAdata[dictKey]['range'].abs()

# load in data from uniprot search
uniprotIDfiles = ['idmapping_2023_11_01.tsv', 'idmapping_2023_11_01_run02.tsv']
allRBVAidMaps = {}
for f in uniprotIDfiles:
    uniprotIDs = pd.read_csv(f, sep="\t")
    allRBVAidMaps = allRBVAidMaps | uniprotIDs.set_index('From')['Entry'].to_dict()

for key in allRBVAdata:
    allRBVAdata[key]['UniProtIDs'] = allRBVAdata[key].index.map(allRBVAidMaps)
    allRBVAdata[key].dropna(inplace = True)

# load in proteomics data
baseFolder = curFolder[:backSlash[-3]]
fileLocation = os.path.join(baseFolder, 'data', 'proteomics',
                            '20230620', 'data_analysis', 'model_comparisons')
protDataGRRfileName = 'mean_not_BCA_imputed_data_with_GRR_sums.csv'
GRRprotDatafile = os.path.join(fileLocation, protDataGRRfileName)
protData = pd.read_csv(GRRprotDatafile, sep='\t')
for c in [c for c in protData.columns if protData[c].dtype in numerics]:
    protData[c] = np.log10(protData[c])
protData.replace(-np.inf, -50, inplace = True) # replace with large negative number to show area in plot


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
        modelKey = setMu + '80mMFA_fullN'
    else:
        if '150' in cond:
            modelKey = setMu + '150mMFA_limN'
        else:
            modelKey = setMu + '80mMFA_limN'
    modelData = allRBVAdata[modelKey]
    reducedModelData = modelData[modelData['UniProtIDs'].str.contains('nan') == False]
    reducedModelData = reducedModelData[['UniProtIDs', 'min', 'max', 'range']]
    mappingDict = reducedModelData.set_index('UniProtIDs').to_dict()
    proteomicsIDs = list(exprByCond[cond]['UniProtIDs'])
    # create list of nans to write into from the mapping dict
    minResults = np.empty((len(proteomicsIDs)))
    minResults[:] = np.nan
    minResults = list(minResults)
    maxResults = minResults.copy()
    rangeResults = minResults.copy()
    mappingDictOrig = mappingDict['min'].copy()
    for key in mappingDictOrig:
        indices = [i for i, s in enumerate(proteomicsIDs) if key in s]
        for idx in indices:
            # make sure that if multiple IDs are mapped all of the results are captured
            # --> now the mean for that is calculated whilst ignoring nan (np.nanmean())
            if isinstance(minResults[idx], list):
                minResults[idx] = minResults[idx].append(mappingDict['min'][key])
                maxResults[idx] = maxResults[idx].append(mappingDict['max'][key])
                rangeResults[idx] = rangeResults[idx].append(mappingDict['range'][key])
            else:
                minResults[idx] = [mappingDict['min'][key]]
                maxResults[idx] = [mappingDict['max'][key]]
                rangeResults[idx] = [mappingDict['range'][key]]
        if len(indices) != 0:
            mappingDict['min'].pop(key)
            mappingDict['max'].pop(key)
            mappingDict['range'].pop(key)
    for i in range(len(minResults)):
        if isinstance(minResults[i], list):
            minResults[i] = np.nanmean(minResults[i])
            maxResults[i] = np.nanmean(maxResults[i])
            rangeResults[i] = np.nanmean(rangeResults[i])
    exprByCond[cond]['RBVA_min'] = minResults
    exprByCond[cond]['RBVA_max'] = maxResults
    exprByCond[cond]['RBVA_range'] = rangeResults
    remainingMin = pd.DataFrame.from_dict(mappingDict['min'], orient='index',
                                          columns = ['RBVA_min'])
    remainingMin = remainingMin.reset_index()
    remainingMin = remainingMin.rename(columns = {'index':'UniProtIDs'})
    exprByCond[cond] = pd.concat([exprByCond[cond], remainingMin], join = 'outer')
    
    remainingMax = pd.DataFrame.from_dict(mappingDict['max'], orient='index',
                                          columns = ['RBVA_max'])
    remainingMax = remainingMax.reset_index()
    remainingMax = remainingMax.rename(columns = {'index':'UniProtIDs'})
    exprByCond[cond] = pd.concat([exprByCond[cond], remainingMax], join = 'outer')
    
    remainingRange = pd.DataFrame.from_dict(mappingDict['range'], orient='index',
                                          columns = ['RBVA_range'])
    remainingRange = remainingRange.reset_index()
    remainingRange = remainingRange.rename(columns = {'index':'UniProtIDs'})
    exprByCond[cond] = pd.concat([exprByCond[cond], remainingRange], join = 'outer')

# drop nan as only values both in the model and the proteomics data should be
# plotted (for now). Also, multiply values by 1000000 to bring them to nmol/gCDW
redExprByCond = {}
for cond in exprByCond:
    redExprByCond[cond] = exprByCond[cond].dropna()
    # redExprByCond[cond][redExprByCond[cond].select_dtypes(include=['number']).columns] *= 1000000
    # rather take the logarithm --> done above
    redExprByCond[cond].sort_values(by = ['UniProtIDs'], inplace = True, ignore_index = True)



#%% rearrange dictionary for plotting as the same condition across differnet strains are ending up in the same graph
plotDict = {}
dfHeights = []
for cond in exprByCond:
    underScore = findOccurrences(cond, '_')
    newKey = cond[underScore[0] + 1 : ]
    if newKey not in plotDict:
        plotDict[newKey] = {}
    plotDict[newKey][cond[ : underScore[0]]] = redExprByCond[cond]
    if 'fullDF' not in plotDict[newKey]:
        plotDict[newKey]['fullDF'] = redExprByCond[cond]
    else:
        plotDict[newKey]['fullDF'] = pd.merge(plotDict[newKey]['fullDF'],
                                               redExprByCond[cond],
                                               on = ['UniProtIDs',
                                                     'RBVA_min',
                                                     'RBVA_max',
                                                     'RBVA_range'],
                                               how = 'outer')
        plotDict[newKey]['fullDF'].drop_duplicates(inplace = True)
        
    # get height of each dataframe to calculate a vector of x
    dfHeights.append(len(plotDict[newKey]['fullDF']))
        


#%% plot data
fig, axs = plt.subplots(3, 1, figsize = (6, 8))
pltNr = 0

font1 = {'family': 'arial', 'color': 'black', 'size': 10}
font2 = {'family': 'arial', 'color': 'black', 'size': 12, 'weight': 'bold'}
plotOrder = ['80_full', '80_lim', '150_lim']
pltTitles = ['80 C-mM, C/N = 5.3',
             '80 C-mM, C/N = 50',
             '150 C-mM, C/N = 50']

for key in plotOrder:
    curAx = axs[pltNr]
    
    plotDict[key]['fullDF'].sort_values('RBVA_range', inplace = True,
                                        ascending = False)
    plotDict[key]['fullDF'].reset_index(drop = True, inplace = True)
    
    # plot area
    ymin = plotDict[key]['fullDF']['RBVA_min']
    ymax = plotDict[key]['fullDF']['RBVA_max']
    x = np.linspace(0, len(ymin)-1, num = len(ymin))
    curAx.fill_between(x, ymin, ymax, alpha=0.5)
    # possibly try log
    
    # plot lines as borders for the area
    sns.lineplot(data = plotDict[key]['fullDF'],
                 x = plotDict[key]['fullDF'].index, y = 'RBVA_min',
                 ax = curAx, linewidth = 1,
                 color = sns.color_palette()[0])
    rbvaData = sns.lineplot(data = plotDict[key]['fullDF'],
                            x = plotDict[key]['fullDF'].index, y = 'RBVA_max',
                            ax = curAx, linewidth = 1,
                            color = sns.color_palette('colorblind')[0])
    
    # plot scatterplot for proteomics values
    # find column to plot
    for col in list(plotDict[key]['fullDF'].columns):
        if'mean' in col:
            if 'H16' in col:
                H16col = col
            elif 'ALE26' in col:
                ALE26col = col
    
    h16Data = sns.scatterplot(data = plotDict[key]['fullDF'],
                              x = plotDict[key]['fullDF'].index, y = H16col,
                              ax = curAx,
                              color = sns.color_palette('colorblind')[2])
    
    ale26Data = sns.scatterplot(data = plotDict[key]['fullDF'],
                                x = plotDict[key]['fullDF'].index,
                                y = ALE26col,
                                ax = curAx,
                                color = sns.color_palette('colorblind')[3])
    
    curAx.set_title(pltTitles[pltNr], fontdict = font2)
    curAx.set_ylabel('log($c_{protein}$ [mmol/$g_{CDW}$])', fontdict = font1)
    curAx.tick_params(axis = 'x',
                      which = 'both',
                      bottom = True,
                      top = False,
                      labelbottom = False)
    curAx.set_ylim(-12, 0)
    curAx.set_xlim(left = 0)
    
    pltNr += 1
    if pltNr == len(plotOrder):
        curAx.set_xlabel('proteins', fontdict = font1)
        children = curAx.get_children()

fig.legend([children[0], children[5], children[6]],
                 ['RBVA range', 'H16', 'ALE26'], loc = 'right',
                 bbox_to_anchor=(1.15, 0.5))
#%% save the figure
figName = 'RBVA_proteomics_comp_thesis'
fullPlotSavePath = os.path.join(fileLocation, 'figures', figName + '.svg')
fig.figure.savefig(fullPlotSavePath, format='svg')

fullPlotSavePath = os.path.join(fileLocation, 'figures', figName + '.png')
fig.figure.savefig(fullPlotSavePath, bbox_inches='tight', dpi = 1000)