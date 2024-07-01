# -*- coding: utf-8 -*-
"""
Created on Thu Aug 31 16:31:50 2023

@author: Manuel Bruch
relative filepaths corrected for all OS
"""

#%% clean up
from IPython import get_ipython
get_ipython().magic('reset -sf') # clears all variables

#%% import statements
import os
import numpy as np
import pandas as pd
from functools import reduce

#%% function definitions
def findOccurrences(s, ch):
    return [i for i, letter in enumerate(s) if letter == ch]


#%% calculate factor to convert LFQ intensities into protein masses - roughly
CDWdict = {'ALE26_80_full':   {'CDW_mg': 12.9,
                               'V_CDW_harvested_ml': 50},
           'ALE26_80_lim':    {'CDW_mg': 28.25,
                               'V_CDW_harvested_ml': 100},
           'ALE26_150_lim':   {'CDW_mg': 44.35,
                               'V_CDW_harvested_ml': 80},
           'H16_80_full':     {'CDW_mg': 16.4,
                               'V_CDW_harvested_ml': 80},
           'H16_80_lim':      {'CDW_mg': 18.35,
                               'V_CDW_harvested_ml': 80},
           'H16_150_lim':     {'CDW_mg': 32.65,
                               'V_CDW_harvested_ml': 80}}
for key in CDWdict:
    CDWdict[key]['CDW_g_L'] = CDWdict[key]['CDW_mg'] / CDWdict[key]['V_CDW_harvested_ml']

protDict = {'ALE26_80_full':  {'V_prot_harvested_ml': 50,
                               'V_resuspension_ul': [700, 700, 400, 400],
                               'c_final_ug_ul': [0.3142, 0.3768, 0.3300, 0.3458],
                               'V_final_ul': 15,
                               'V_used_from_BCA_ul': [20.66, 20.66, 16.08, 16.08],
                               'c_BCA_assay_ug_ul': [4.84, 4.84, 6.22, 6.22]},
            
           'ALE26_80_lim':    {'V_prot_harvested_ml': 50,
                               'V_resuspension_ul': [400, 400, 600, 600],
                               'c_final_ug_ul': [0.3300, 0.3458, 0.2206, 0.3300],
                               'V_final_ul': 15,
                               'V_used_from_BCA_ul': [26.50, 26.50, 30.97, 30.97],
                               'c_BCA_assay_ug_ul': [3.77, 3.77, 3.23, 3.23]},
           
           'ALE26_150_lim':   {'V_prot_harvested_ml': 50,
                               'V_resuspension_ul': [400, 400],
                               'c_final_ug_ul': [0.3616, 0.3300],
                               'V_final_ul': 15,
                               'V_used_from_BCA_ul': [26.37, 26.37],
                               'c_BCA_assay_ug_ul': [3.79, 3.79]},
           
           'H16_80_full':     {'V_prot_harvested_ml': 50,
                               'V_resuspension_ul': [400, 400, 400, 400],
                               'c_final_ug_ul': [0.4564, 0.2856, 0.3312, 0.3160],
                               'V_final_ul': 15,
                               'V_used_from_BCA_ul': [27.53, 27.53, 21.32, 21.32],
                               'c_BCA_assay_ug_ul': [3.63, 3.63, 4.69, 4.69]},
           
           'H16_80_lim':      {'V_prot_harvested_ml': 50,
                               'V_resuspension_ul': [400, 400, 400, 400],
                               'c_final_ug_ul': [0.3306, 0.3464, 0.2844, 0.3154],
                               'V_final_ul': 15,
                               'V_used_from_BCA_ul': [28.56, 28.56, 36.95, 36.95],
                               'c_BCA_assay_ug_ul': [3.50, 3.50, 2.71, 2.71]},
           
           'H16_150_lim':     {'V_prot_harvested_ml': 50,
                               'V_resuspension_ul': [400, 400, 700, 700],
                               'c_final_ug_ul': [0.2522, 0.2358, 0.2048, 0.2304],
                               'V_final_ul': 15,
                               'V_used_from_BCA_ul': [12.85, 12.85, 13.37, 13.37],
                               'c_BCA_assay_ug_ul': [7.80, 7.80, 7.48, 7.48]}}

# calculation for ug protein per L in original sample/reactor
# c_final * V_final = m_final [ug]; m_final / V_BCA = c_Trp; c_Trp * V_res = m_resuspension; m_res / V_orig = c_oig
# (((c_final * V_final) / V_BCA) * V_res) / V_orig
# would be in ug/ul which should be the same as g/L??
for key in protDict:
    protDict[key]['c_orig_g_L'] = []
    protDict[key]['c_orig_g_L_BCA'] = []
    for i in range(len(protDict[key]['c_final_ug_ul'])):
        c_final = protDict[key]['c_final_ug_ul'][i]
        V_final = protDict[key]['V_final_ul']
        V_BCA = protDict[key]['V_used_from_BCA_ul'][i]
        V_res = protDict[key]['V_resuspension_ul'][i]
        V_orig = protDict[key]['V_prot_harvested_ml'] * 1000 # conversion to ul
        c_BCA = protDict[key]['c_BCA_assay_ug_ul'][i]
        c_orig = (((c_final * V_final) / V_BCA) * V_res) / V_orig
        c_orig_BCA = c_BCA * V_res / V_orig
        protDict[key]['c_orig_g_L'].append(c_orig)
        protDict[key]['c_orig_g_L_BCA'].append(c_orig_BCA)
    protDict[key]['c_orig_g_L_avg'] = np.average(protDict[key]['c_orig_g_L'])
    protDict[key]['c_orig_g_L_avg_BCA'] = np.average(protDict[key]['c_orig_g_L_BCA'])

convFactors = {}
for key in protDict:
    convFactors[key] = {}
    convFactors[key]['g_prot/g_CDW'] = protDict[key]['c_orig_g_L_avg'] / CDWdict[key]['CDW_g_L']
    convFactors[key]['g_prot_BCA/g_CDW'] = protDict[key]['c_orig_g_L_avg_BCA'] / CDWdict[key]['CDW_g_L']

# the factor calculated with the concentrations estimated from the BCA assay
# performed during the sample prep is closer to Michael's 
# (0.65 g(Protein)/g(CDW)). Likely in the other steps there is too much 
# proteome lost --> use that factor for the conversion

#%% load data
"""
load in data from the imputed tables from Perseus
load in molecular weight data from UniProt search
"""
curFolder = os.getcwd()
backSlash = findOccurrences(curFolder, os.sep)
baseFolder = curFolder[:backSlash[-1]]

allImputedDataLocation = os.path.join(baseFolder, 'data', 'proteomics',
                                '20230620', 'data_analysis', 'allImputedData_dataframe.csv')

if not os.path.exists(allImputedDataLocation):
    excelFile = os.path.join(baseFolder, 'data', 'proteomics', 
                             '20230620', 'data_analysis',
                             '20230725_imputed_proteins_data_prep_for_PRISM.xlsx')
    imputedALEData = pd.read_excel(excelFile, 'export_ALE26_only_imputed')
    imputedWTData = pd.read_excel(excelFile, 'export_H16_only_imputed')
    imputedCommonData = pd.read_excel(excelFile, 'export_proteins_exp_in_both_str')
    
    # prepare master dataframe
    # get groups for later mean calculation
    allConds = list(imputedCommonData.loc[3])
    allConds = set([i for i in allConds if str(i) != 'nan'])
    
    # remove first four rows from dataframes
    imputedALEData.drop([0, 1, 2, 3], inplace = True)
    imputedWTData.drop([0, 1, 2, 3], inplace = True)
    imputedCommonData.drop([0, 1, 2, 3], inplace = True)
    # concatenate dataframes
    allImputedData = pd.concat([imputedCommonData, imputedWTData])
    allImputedData = pd.concat([allImputedData, imputedALEData])
    allImputedData.to_csv(allImputedDataLocation, sep='\t', index = False)
else:
    allImputedData = pd.read_csv(allImputedDataLocation, sep="\t")

# load UniProt. Mass is in Da (g/mol)
uniprotFile = os.path.join(baseFolder, 'data', 'proteomics',
                           'UniProt_idmapping_2023_09_01.tsv')
uniProtData = pd.read_csv(uniprotFile, sep="\t")
massDict = pd.Series(uniProtData.Mass.values, index = uniProtData.From).to_dict()

#%% convert data to g/gCDW
fileName = os.path.join(baseFolder, 'data', 'proteomics',
                        '20230620', 'data_analysis', 'model_comparisons',
                        'allImputedData_mmol_p_g_CDW_BCA.csv')
if not os.path.exists(fileName):
    dataDfColumns = list(allImputedData.columns)
    for i in range(len(dataDfColumns)):
        dataDfColumns[i] = dataDfColumns[i].replace('H16 full 80', 'H16_80_full')
        dataDfColumns[i] = dataDfColumns[i].replace('H16 CN50 80', 'H16_80_lim')
        dataDfColumns[i] = dataDfColumns[i].replace('H16 CN50 150', 'H16_150_lim')
        dataDfColumns[i] = dataDfColumns[i].replace('ALE26 full 80', 'ALE26_80_full')
        dataDfColumns[i] = dataDfColumns[i].replace('ALE26 CN50 80', 'ALE26_80_lim')
        dataDfColumns[i] = dataDfColumns[i].replace('ALE26 CN50 150', 'ALE26_150_lim')
    allImputedData.set_axis(dataDfColumns, axis=1,inplace=True)
    dataDfColumns = [i for i in dataDfColumns if 'Protein' not in i]
    allConditions = set([i[:-3].replace('data ','') for i in dataDfColumns])
    
    # data are currently logarithmic from analysis in Perseus
    # require exponential transformation
    for col in dataDfColumns:
        allImputedData[col] = np.exp(allImputedData[col])
    
    dataColSums = allImputedData.sum()
    
    convertedData = {'g_prot/g_CDW':[],
                     'g_prot_BCA/g_CDW': []}
    for outerKey in convertedData:
        convertedData[outerKey] = allImputedData.copy()
        for col in dataDfColumns:
            # find the correct conversion factor, depending on condition
            for key in convFactors:
                if key in col:
                    convFactor = convFactors[key][outerKey]
                    break
            # divide entire column by its sum (fraction of protein in whole proteome)
            # Question: should I use the non-imputed data here? Otherwise the sum will be too high?
            convertedData[outerKey][col] = convertedData[outerKey][col].div(dataColSums[col])
            
            convertedData[outerKey][col] = convertedData[outerKey][col].multiply(convFactor)
    
    #%% convert to mmol/gCDW row-wise
    for key in convertedData:
        for index, row in convertedData[key].iterrows():
            protID = row['Protein IDs']
            # use average molecular weight in case of multiple IDs
            if ';' in protID:
                protIDs = protID.split(';')
                allMolWeight = []
                for ID in protIDs:
                    try:
                        allMolWeight.append(massDict[ID])
                    except:
                        continue
                molWeight = np.average(allMolWeight)
            else:
                molWeight = massDict[protID]
            for col in dataDfColumns:
                row[col] = (row[col] / molWeight) * 1000 #conversion to mmol
            convertedData[key].loc[index] = row
    
    #%% save the dataframes
    for key in convertedData:
        if 'BCA' in key:
            fileName = os.path.join(baseFolder, 'data', 'proteomics',
                                    '20230620', 'data_analysis', 'model_comparisons',
                                    'allImputedData_mmol_p_g_CDW_BCA.csv')
        else:
            fileName = os.path.join(baseFolder, 'data', 'proteomics',
                                    '20230620', 'data_analysis', 'model_comparisons',
                                    'allImputedData_mmol_p_g_CDW_not_BCA.csv')
        convertedData[key].to_csv(fileName, sep='\t', index = False)
    #%% load the converted dataframes with mmol(protein)/gCDW data
else:
    convertedData = {'g_prot/g_CDW':[],
                     'g_prot_BCA/g_CDW': []}
    for key in convertedData:
        if 'BCA' in key:
            fileName = os.path.join(baseFolder, 'data', 'proteomics',
                                    '20230620', 'data_analysis', 'model_comparisons',
                                    'allImputedData_mmol_p_g_CDW_BCA.csv')
        else:
            fileName = os.path.join(baseFolder, 'data', 'proteomics',
                                    '20230620', 'data_analysis', 'model_comparisons',
                                    'allImputedData_mmol_p_g_CDW_not_BCA.csv')
        convertedData[key] = pd.read_csv(fileName, sep="\t")

#%% calculate averages for each condition
saveFilePath = os.path.join(baseFolder, 'data', 'proteomics',
                            '20230620', 'data_analysis', 'model_comparisons')
meanDFs = {}
for key in convertedData:
    workingDF = convertedData[key]
    dataDfColumns = list(workingDF.columns)
    dataDfColumns = [i for i in dataDfColumns if 'Protein' not in i]
    allConditions = set([i[:-3].replace('data ','') for i in dataDfColumns])
    # split dataframe into grouped dfs
    condDfs = []
    for cond in allConditions:
        usedCols = [i for i in dataDfColumns if cond in i]
        usedCols = ['Protein IDs'] + usedCols
        thisCondDf = workingDF[usedCols]
        meanCondDf = pd.DataFrame(thisCondDf['Protein IDs'])
        meanCondDf[cond + '_mean_mmol_gCDW'] = thisCondDf.mean(axis=1)
        meanCondDf[cond + '_std_mmol_gCDW'] = thisCondDf.std(axis=1)
        condDfs.append(meanCondDf)
    mergedMeanDf = reduce(lambda x,y: pd.merge(x,y, on='Protein IDs', how='outer'), condDfs)
    meanDFs[key] = mergedMeanDf
    
    # save dataframe
    if 'BCA' in key:
        saveFileName = 'mean_allImputedData_mmol_p_g_CDW_BCA.csv'
    else:
        saveFileName = 'mean_allImputedData_mmol_p_g_CDW_not_BCA.csv'
    # meanDFs[key].to_csv(os.path.join(saveFilePath, saveFileName),
    #                     sep='\t', index = False)