# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 14:55:10 2024

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

# Importing PCA
from sklearn.decomposition import PCA

#%% set parameters
printComparisons = False
plotAll = False
logPlot = False

#%% function definitions
def findOccurrences(s, ch):
    return [i for i, letter in enumerate(s) if letter == ch]


#%% load in data
curFolder = os.getcwd()
backSlash = findOccurrences(curFolder, os.sep)
baseFolder = curFolder[:backSlash[-1]]

fileLocation = os.path.join(baseFolder, 'data', 'proteomics', '20230620',
                            'data_analysis')
figSaveLocation = os.path.join(fileLocation, 'PCA')
figName = 'proteomics_PCA'
fileName = ['export_proteins_exp_in_both_strains_after_imputation.txt',
            'export_H16_only_imputed.txt',
            'export_ALE26_only_imputed.txt']
variant = ['Proteins in Both Strains',
           'Proteins in H16 only',
           'Proteins in ALE26 only']
varIdx = 0

allRawData = []
allProtData = []
allZ = []
allPCAdf = []

for fn in fileName:
    fullFile = os.path.join(fileLocation, fn)
    
    rawProtData = pd.read_csv(fullFile, sep = '\t')
    allRawData.append(rawProtData)
    # remove categorical rows
    protData = rawProtData.drop([0, 1, 2, 3, 4])
    
    # remove unwanted columns
    keepCols = [i for i in list(protData.columns) if any(ele in i for ele in ['H16', 'ALE26', 'Protein IDs'])]
    protData = protData[keepCols]
    
    protData = protData.set_index('Protein IDs')
    # convert all data to float
    for col in protData.columns:
        protData[col] = protData[col].astype(float)
    allProtData.append(protData)
    
    #%% manipulate data for PCA
    transposedData = protData.transpose() # PCA assumes features are in columns and samples in rows
    
    # Data need to be normalised using mean and standard deviation
    transposedData_mean = transposedData.mean()
    transposedData_std = transposedData.std()
    # Standardization
    Z = (transposedData - transposedData_mean) / transposedData_std
    allZ.append(Z)
    
    #%% PCA
    nComponents = 2
    
    pca = PCA(n_components = nComponents)
    pca.fit(Z)
    x_pca = pca.transform(Z)
     
    # Create the dataframe
    df_pca1 = pd.DataFrame(x_pca,
                           columns=['PC{}'.
                           format(i+1)
                            for i in range(nComponents)])
    
    # extend PCA table with some identifiers for the plotting later
    sampleNames = list(Z.index)
    df_pca1['sample'] = sampleNames
    trimmedNames = [i[5:-3] for i in sampleNames]
    df_pca1['trimmedNames'] = trimmedNames
    
    strain = []
    CNratio = []
    carbon = []
    replID = []
    for smpl in sampleNames:
        nameComponents = smpl.split()
        strain.append(nameComponents[1])
        CNratio.append(nameComponents[2])
        carbon.append(nameComponents[3])
        replID.append(nameComponents[4])
        
    df_pca1['strain'] = strain
    df_pca1['CNratio'] = CNratio
    df_pca1['carbon'] = carbon
    df_pca1['replID'] = replID
    
    dataSet = [variant[varIdx]] * len(sampleNames)
    varIdx += 1
    df_pca1['dataSet'] = dataSet
    
    
    allPCAdf.append(df_pca1)

combinedPCAdf = pd.concat(allPCAdf)

#%% simple plot
# giving a larger plot
 
# plt.scatter(x_pca[:, 0], x_pca[:, 1],
#             c = df_pca1['sample'],
#             cmap='plasma')
with sns.plotting_context(rc={"legend.fontsize": 10}):
    g = sns.relplot(data = combinedPCAdf, x = 'PC1', y = 'PC2',
                    hue = 'trimmedNames', col = "dataSet", palette = 'colorblind',
                    height = 2.5, aspect = 0.7)

# modify plot appearance
g._legend.set_title("")

axCount = 0
font1 = {'family': 'arial', 'color': 'black', 'size': 10}
font2 = {'family': 'arial', 'color': 'black', 'weight': 'bold', 'size': 12}
for ax in g.axes.flatten():
    pltTitle = variant[axCount].split()
    pltTitle = pltTitle[-2] + ' ' + pltTitle[-1]
    ax.set_title(pltTitle, fontdict = font2)
    ax.set_xlabel("PC1", fontdict = font1)
    if axCount == 0:
        ax.set_ylabel("PC2", fontdict = font1)
    
    axCount += 1
 
#%% save the figure
figureSavePath = os.path.join(figSaveLocation,
                              figName + '.svg')
g.figure.savefig(figureSavePath, format='svg')

figureSavePath = os.path.join(figSaveLocation,
                              figName + '.png')
g.figure.savefig(figureSavePath, bbox_inches='tight', dpi = 500)