# -*- coding: utf-8 -*-
"""
Created on Mon Jul 31 10:25:49 2023

@author: Manuel Bruch
relative filepaths corrected for all OS
"""

#%% clean up
from IPython import get_ipython
get_ipython().magic('reset -sf') # clears all variables

#%% import statements
import os
import math
import numpy as np
import pandas as pd

#%% function definitions
def findOccurrences(s, ch):
    return [i for i, letter in enumerate(s) if letter == ch]

#%% load in data
curFolder = os.getcwd()
backSlash = findOccurrences(curFolder, os.sep)
baseFolder = curFolder[:backSlash[-1]]

#%% set parameters for multifamily or singlefamily comparisons
desComp = 'singleFam'
dataSet = 'ALE26_only'

if 'both_strains' in dataSet:
    fileNameStart = '2way ANOVA of proteins_in_both_strains_data_after_imputation_'
elif 'H16_only' in dataSet:
    fileNameStart = '2way ANOVA of proteins_in_H16_'
elif 'ALE26_only' in dataSet:
    fileNameStart = '2way ANOVA of proteins_in_ALE26_'

if 'multiFam' in desComp:
    if 'both_strains' in dataSet:
        fileFolder = 'multiFam'
    elif 'H16_only' in dataSet:
        fileFolder = 'H16_only'
    elif 'ALE26_only' in dataSet:
        fileFolder = 'ALE26_only'
    fileLocation = os.path.join(baseFolder, 'data', 'proteomics', 
                                '20230620', 'data_analysis', 'PRISM', fileFolder)
    fileName = fileNameStart + 'multi_family.txt'
    perfTest = "Tukey's multiple comparisons test"
    lowerTableAddStatistic = 'q'
    resultFileNameStart = 'multiFam_'
elif 'singleFam' in desComp:
    if 'both_strains' in dataSet:
        fileFolder = 'singleFam'
    elif 'H16_only' in dataSet:
        fileFolder = 'H16_only'
    elif 'ALE26_only' in dataSet:
        fileFolder = 'ALE26_only'
    fileLocation = os.path.join(baseFolder, 'data', 'proteomics', 
                                '20230620', 'data_analysis', 'PRISM', fileFolder)
    fileName = fileNameStart + 'single_family.txt'
    perfTest = "Šídák's multiple comparisons test"
    lowerTableAddStatistic = 't'
    resultFileNameStart = 'singleFam_'

dataFile = os.path.join(fileLocation, fileName)

rawData = pd.read_csv(dataFile, sep="\t")

#%% formate the data frame
# extract only section of the Tukey results I am interested in
splitRowStart = rawData[rawData['Unnamed: 0']==perfTest].index.tolist()[0]
splitRowEnd = rawData[rawData['Unnamed: 0']=='Test details'].index.tolist()[0]
testResultData = rawData.iloc[splitRowStart + 2 : splitRowEnd]
# remove leading and trailing spaces from strings in cells
testResultData = testResultData.applymap(lambda x: x.strip() if isinstance(x, str) else x)
# make group naming equal for all
testResultData['Unnamed: 0'] = testResultData['Unnamed: 0'].str.replace('H16 C/N = 50 80 mM formic acid','Group B')
# remove columns of all NaNs
testResultData = testResultData.dropna(axis=1,how='all')
# change headers
testResultData = testResultData.set_axis([perfTest,
                                "Predicted (LS) mean diff.",
                                "95.00% CI of diff.", "Below threshold?",
                                "Summary", "Adjusted P Value"], axis=1)
# replace all "greater than" and "less than" values in P-values
mapping = {'>0.9999': 1, '<0.0001': 0.0001}
testResultData = testResultData.replace({'Adjusted P Value': mapping})
# change column datatypes
testResultData["Adjusted P Value"] = pd.to_numeric(testResultData["Adjusted P Value"], errors='coerce')
testResultData["Predicted (LS) mean diff."] = pd.to_numeric(testResultData["Predicted (LS) mean diff."], errors='coerce')

#%% realised I need the lower part of the table as well to calculate the log of the fold change
lowerDF = rawData.iloc[splitRowEnd+2:]
lowerDF = lowerDF.applymap(lambda x: x.strip() if isinstance(x, str) else x)
lowerDF['Unnamed: 0'] = lowerDF['Unnamed: 0'].str.replace('H16 C/N = 50 80 mM formic acid','Group B')
lowerDF = lowerDF.set_axis(["Test details", "Predicted (LS) mean 1", 
                            "Predicted (LS) mean 2",
                            "Predicted (LS) mean diff.", "SE of diff.", "N1",
                            "N2", lowerTableAddStatistic, "DF"], axis=1)
colsToDrop = ["Predicted (LS) mean diff.", "SE of diff.", "N1", "N2",
              lowerTableAddStatistic, "DF"]
lowerDF = lowerDF.drop(columns = colsToDrop)
lowerDF["Predicted (LS) mean 1"] = pd.to_numeric(lowerDF["Predicted (LS) mean 1"], errors='coerce')
lowerDF["Predicted (LS) mean 2"] = pd.to_numeric(lowerDF["Predicted (LS) mean 2"], errors='coerce')

# assign protein IDs as values in rows
lowerDF['ProtID'] = ''
for index, row in lowerDF.iterrows():
    if math.isnan(row["Predicted (LS) mean 1"]):
        curProtID = row["Test details"]
    lowerDF['ProtID'].loc[index] = curProtID

lowerDF = lowerDF.sort_values(by = ["Test details"])
# remove rows that don't contain test results
lowerDF = lowerDF[lowerDF["Test details"].str.contains("Group", na=False)]
# no need to take the log of the values for subtraction as the values already are log transformed
lowerDF["log(fold change)"] = lowerDF["Predicted (LS) mean 1"]-lowerDF["Predicted (LS) mean 2"]

# further trim the dataframe
colsToDrop = ["Predicted (LS) mean 1", "Predicted (LS) mean 2"]
lowerDF = lowerDF.drop(columns = colsToDrop)

#%% format table into multiple dataframes, one per comparison
# assign protein IDs as values in rows
testResultData['ProtID'] = ''
for index, row in testResultData.iterrows():
    if math.isnan(row["Predicted (LS) mean diff."]):
        curProtID = row[perfTest]
    testResultData['ProtID'].loc[index] = curProtID

# sort the dataframe by the comparisons
sortedData = testResultData.sort_values(by = [perfTest])
# only retain rows that belong to a comparison
sortedData = sortedData[sortedData[perfTest].str.contains("Group", na=False)]

# convert p to -log(p)
sortedData["-log(adjusted p)"] = -1*np.log(sortedData["Adjusted P Value"])

# trim the dataframe of unneccessary columns
colsToDrop = ["Predicted (LS) mean diff.", "95.00% CI of diff.",
              "Below threshold?", "Summary", "Adjusted P Value"]
sortedData = sortedData.drop(columns = colsToDrop)
sortedData.rename(columns={perfTest:"Test details"}, inplace=True)

#%% merge the two dataframes
finalDF = pd.merge(sortedData, lowerDF,  how='left', left_on=["Test details", 'ProtID'], right_on = ["Test details", 'ProtID'])

#%% split the dataframe
# get unique names of groups (comparisons)
UniqueNames = finalDF["Test details"].unique()
DataFrameDict = {elem : pd.DataFrame() for elem in UniqueNames}
colsToDrop = ["Test details"]

for key in DataFrameDict.keys():
    DataFrameDict[key] = finalDF[:][finalDF["Test details"] == key]
    DataFrameDict[key] = DataFrameDict[key].drop(columns = colsToDrop)
    DataFrameDict[key] = DataFrameDict[key].sort_values(by = ['ProtID'])
    DataFrameDict[key] = DataFrameDict[key][["log(fold change)", "-log(adjusted p)", "ProtID"]]

#%% save the dataframes
for key in DataFrameDict.keys():
    condName = key
    condName = condName.replace('.','')
    condName = condName.replace(' ','_')
    resultFileName = resultFileNameStart + condName + '.csv'
    DataFrameDict[key].to_csv(os.path.join(fileLocation, resultFileName), sep='\t', index = False)