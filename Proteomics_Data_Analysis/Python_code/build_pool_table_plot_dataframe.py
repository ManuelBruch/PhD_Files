# -*- coding: utf-8 -*-
"""
Created on Fri Aug 18 14:53:29 2023

@author: Manuel
relative filepaths corrected for all OS
"""
#%% clean up
from IPython import get_ipython
get_ipython().magic('reset -sf') # clears all variables

#%% import statements
import io
import os
import pickle
import re

import pandas as pd
import numpy as np

from Bio.KEGG import REST

#%% function definitions
def findOccurrences(s, ch):
    return [i for i, letter in enumerate(s) if letter == ch]


def string_list_prep(originalList):
    searchList = originalList.copy() # copy to not change the original variable
    allTheLists = []
    multiEntryIndeces = []
    while True: # loop as long as there are semicolons in the list, i.e. there are more multi-entry proteins
        goAgain = False
        backupList = np.empty(len(originalList))
        backupList[:] = np.nan
        backupList = list(backupList)
        curMultiEntryIndeces = []
        for i in range(len(searchList)):
            try:
                if ';' in searchList[i]:
                    semPos = findOccurrences(searchList[i], ';')
                    if len(semPos) > 1:
                        goAgain = True
                    backupList[i] = searchList[i][semPos[0]+1:]
                    searchList[i] = searchList[i][:semPos[0]]
                    curMultiEntryIndeces.append(i)
            except:
                continue
        if i == len(searchList)-1:
            cleanedList = [x for x in searchList if not pd.isnull(x)]
            allTheLists.append(cleanedList)
            multiEntryIndeces.append(curMultiEntryIndeces)
            if not goAgain:
                cleanedList = [x for x in backupList if not pd.isnull(x)]
                allTheLists.append(cleanedList)
                break
            elif goAgain:
                searchList = backupList
    return allTheLists, multiEntryIndeces


def to_df(result):
    return pd.read_table(io.StringIO(result), header=None)


def get_kegg_category(categoryList, categoryIndexList, curProt, keggCategory):
    startIdx = np.nan
    if 'PATHWAY' in keggCategory:
        reStr = 'reh\d*'
    elif 'MODULE' in keggCategory:
        reStr = 'reh_M\d*'
    for j in range(len(categoryList)):
        if keggCategory in categoryList[j]:
            startIdx = categoryIndexList[j]
            endIdx = categoryIndexList[j+1]
            break
    # get the pathways
    if not np.isnan(startIdx):
        kegg_cat = curProt[startIdx:endIdx]
        for k in range(len(kegg_cat)):
            kegg_cat[k] = kegg_cat[k].replace(keggCategory,'')
            kegg_cat[k] = re.sub(reStr, '', kegg_cat[k])
            kegg_cat[k] = kegg_cat[k].strip()
    else:
        kegg_cat = []
    return kegg_cat


def parse_kegg(pathwayDict, subIDList, curProtList):
    for i in range(len(subIDList)):
        curProt = curProtList[i].split('\n')
        categoryList = []
        categoryIndexList = []
        pathwayDict[subIDList[i]] = {}
        for j in range(len(curProt)):
            try:
                if not ' ' in curProt[j][0]:
                    categoryList.append(curProt[j])
                    categoryIndexList.append(j)
            except:
                continue
        pathwayDict[subIDList[i]]['Pathways'] = get_kegg_category(categoryList, categoryIndexList, curProt, 'PATHWAY')
        pathwayDict[subIDList[i]]['Pathway_Module'] = get_kegg_category(categoryList, categoryIndexList, curProt, 'MODULE')
    return pathwayDict


def combine_lists(baseList, baseListIndex, addedItem, listDelimiter, removeDuplicates = False):
    if len(baseList[baseListIndex]) == 0:
        baseList[baseListIndex] = addedItem
    else:
        baseList[baseListIndex] = listDelimiter.join([baseList[baseListIndex], addedItem])
    if removeDuplicates:
        curItem = baseList[baseListIndex].split(listDelimiter)
        uniqueItems = set(curItem)
        uniqueItems = list(uniqueItems)
        baseList[baseListIndex] = listDelimiter.join(uniqueItems)
    return baseList


#%% read data frames with proteomics data
curFolder = os.getcwd()
backSlash = findOccurrences(curFolder, os.sep)
baseFolder = curFolder[:backSlash[-1]]

csvFileLocation = os.path.join(baseFolder, 'data', 'proteomics', 
                            '20230620', 'data_analysis')
dataFile = 'proteomics_condition_count.csv'

dataTable = pd.read_csv(os.path.join(csvFileLocation, dataFile), index_col = 0,
                        sep="\t")
PTdf = dataTable[['proteinIDs',
                  'EntryName',
                  'ProteinName',
                  'geneLocations']].copy()
PTdf.reset_index(inplace = True) # just to make sure that the index is reset
PTdf.rename(columns={'proteinIDs':'UniProtIDs'}, inplace=True)

#%% repackage protein ID list as in differential_expression_analysis.py
allProtIDs = list(PTdf["UniProtIDs"])
allProtSearchLists, allDuplicateEntryPos = string_list_prep(allProtIDs)
# both together
searchParamDict = {}
searchParamDict['searchList'] = allProtSearchLists
searchParamDict['entryPos'] = allDuplicateEntryPos
keggIDs = ["not found"]*len(searchParamDict['searchList'][0])

#%% write ID list into csv file --> required for manual uniprot search
myProtIDsFile = 'all_detected_protein_IDs.csv'
targetFile = os.path.join(baseFolder, 'data', 'proteomics', 
                            myProtIDsFile)
if not os.path.isfile(targetFile):
    for i in range(len(searchParamDict['searchList'])):
        if i == 0:
            with open(targetFile, 'x') as protIDFile:
                for item in searchParamDict['searchList'][0]:
                    protIDFile.write(item+"\n")
        else:
            with open(targetFile, 'a') as protIDFile:
                for item in searchParamDict['searchList'][i]:
                    protIDFile.write(item+"\n")

# read excel table manually downloaded based on submitted IDs
uniProtFile = os.path.join(baseFolder, 'data', 'proteomics',
                           'UniProt_idmapping_2023_08_21.xlsx')
UniProtTable = pd.read_excel(uniProtFile)
UniProtTable["Pathway"] = UniProtTable["Pathway"].str.replace("PATHWAY: ", "")
UniProtTable['KEGG'] = UniProtTable['KEGG'].str.rstrip(';')

# extract KEGG IDs
keggIDs = list(UniProtTable['KEGG'])
keggIDs_split = [item.split(';') for item in keggIDs]
keggIDs_flat = [item for l in keggIDs_split for item in l]

#%% database search for pathways
#%% 1. UniProt
# done manually by downloading an excel table, imported above

#%% 2. KEGG
# create dict with IDs as keys and pathways as values
if os.path.isfile('pathwayDict.pkl'):
    # load the dictionary
    with open('pathwayDict.pkl', 'rb') as f:
        pathwayDict = pickle.load(f)
else:
    pathwayDict = {}
    while True:
        try:
            subIDList = keggIDs_flat[0:10]
            if len(subIDList) == 0:
                break
        except:
            subIDList = keggIDs_flat[0:]
        keggEntries = REST.kegg_get(subIDList).read() # maximum amount of entries that can be retrieved is 10
        # on the laptop I now get this: <urlopen error [SSL: CERTIFICATE_VERIFY_FAILED] certificate verify failed: unable to get local issuer certificate (_ssl.c:1129)>
        curProtList = keggEntries.split('///') # split the list in the individual entries. As the former text ends in the split string, this list will be 11 entries long, the last of which just has a space inside
        pathwayDict = parse_kegg(pathwayDict, subIDList, curProtList)
        try:
            keggIDs_flat = keggIDs_flat[10:]
        except:
            break
    
    # save the dictionary as this ran for about 10 minutes. Also, there are 1761 IDs without associated pathway(s)
    # in the UniProtTable, 2579 entries do not have associated pathways
    
    # Saving the objects:
    with open('pathwayDict.pkl', 'wb') as f:  # Python 3: open(..., 'wb')
        pickle.dump(pathwayDict, f)

#%% time to build the actual table
# table should have the following columns
# UniProtIDs | EntryName | ProteinName | geneLocations | UniProtPathways | keggIDs | keggPathways | keggPathwayModules
# followed by the statistical results, for each comparison one column for p-values and one for log(difference)
# for the pathways: convert lists into simple strings, seperate pathways with ";"
# will need to later worry about sorting through pathways to only have relevant identifiers

# prepare lists to add to the uniprot table from the kegg results
UniProtKeggIdList = UniProtTable['KEGG']
allKeggSearchLists, allDuplicateKeggEntryPos = string_list_prep(UniProtKeggIdList)

# several pathway descriptions should be excluded from the pathway dict from KEGG as they are too general
ignoredPathways = {'Metabolic pathways',
                   'Microbial metabolism in diverse environments',
                   'Biosynthesis of secondary metabolites'}
for key in pathwayDict:
    pathwayList = pathwayDict[key]['Pathways']
    moduleList = pathwayDict[key]['Pathway_Module']
    if len(pathwayList) > 0:
        pathwayList = [i for i in pathwayList if i not in ignoredPathways]
        pathwayDict[key]['Pathways'] = '//'.join(pathwayList)
    else:
        pathwayDict[key]['Pathways'] = ''
    if len(moduleList) > 0:
        pathwayDict[key]['Pathway_Module'] = '//'.join(moduleList)
    else:
        pathwayDict[key]['Pathway_Module'] = ''

# populate a kegg pathway list from the first list in allKeggSearchLists
keggPathways = []
keggPathwayModules = []
for entry in allKeggSearchLists[0]:
    keggPathways.append(pathwayDict[entry]['Pathways'])
    keggPathwayModules.append(pathwayDict[entry]['Pathway_Module'])

# remove the first list from the kegg list of lists of identifiers
allKeggSearchLists = allKeggSearchLists[1:]

# loop through the remaining lists to further populate the keggPathways and keggPathwayModules lists
for i in range(len(allKeggSearchLists)):
    for j in range(len(allKeggSearchLists[i])):
        keggPathways = combine_lists(keggPathways, allDuplicateKeggEntryPos[i][j], pathwayDict[allKeggSearchLists[i][j]]['Pathways'], '//', removeDuplicates = True)
        keggPathwayModules = combine_lists(keggPathwayModules, allDuplicateKeggEntryPos[i][j], pathwayDict[allKeggSearchLists[i][j]]['Pathway_Module'], '//', removeDuplicates = True)

# add these lists as columns to the UniProtTable
UniProtTable['KEGG_Pathways'] = keggPathways
UniProtTable['KEGG_Pathway_Modules'] = keggPathwayModules
# replace all nan values in the uniprot pathways with an empty string akin to the kegg pathways
UniProtTable['Pathway'] = UniProtTable['Pathway'].fillna('')

# create dictionaries for each of the desired dataframe columns
uniProtPathwayDict = pd.Series(UniProtTable['Pathway'].values, index = UniProtTable['From']).to_dict()
keggPathwayDict = pd.Series(UniProtTable['KEGG_Pathways'].values, index = UniProtTable['From']).to_dict()
keggModuleDict = pd.Series(UniProtTable['KEGG_Pathway_Modules'].values, index = UniProtTable['From']).to_dict()
keggIDDict = pd.Series(UniProtTable['KEGG'].values, index = UniProtTable['From']).to_dict()

# then make lists for the PTdf table
uniProtPathways = []
uniprotKeggPathways =[]
uniprotKeggModules = []
uniprotKeggIDs = []

for entry in allProtSearchLists[0]:
    uniProtPathways.append(uniProtPathwayDict[entry])
    uniprotKeggPathways.append(keggPathwayDict[entry])
    uniprotKeggModules.append(keggModuleDict[entry])
    uniprotKeggIDs.append(keggIDDict[entry])

# remove the first list from the kegg list of lists of identifiers
allProtSearchLists = allProtSearchLists[1:]

for i in range(len(allProtSearchLists)):
    for j in range(len(allProtSearchLists[i])):
        if len(uniProtPathwayDict[allProtSearchLists[i][j]]) > 0:
            uniProtPathways = combine_lists(uniProtPathways, allDuplicateEntryPos[i][j], uniProtPathwayDict[allProtSearchLists[i][j]], '///', removeDuplicates = True)
        if len(keggPathwayDict[allProtSearchLists[i][j]]) > 0:
            uniprotKeggPathways = combine_lists(uniprotKeggPathways, allDuplicateEntryPos[i][j], keggPathwayDict[allProtSearchLists[i][j]], '///', removeDuplicates = True)
        if len(keggModuleDict[allProtSearchLists[i][j]]) > 0:
            uniprotKeggModules = combine_lists(uniprotKeggModules, allDuplicateEntryPos[i][j], keggModuleDict[allProtSearchLists[i][j]], '///', removeDuplicates = True)
        if len(keggIDDict[allProtSearchLists[i][j]]) > 0:
            uniprotKeggIDs = combine_lists(uniprotKeggIDs, allDuplicateEntryPos[i][j], keggIDDict[allProtSearchLists[i][j]], '///', removeDuplicates = True)

#%% add pathways to the PTdf table
PTdf['UniProt_Pathways'] = uniProtPathways
PTdf['KEGG_IDs'] = uniprotKeggIDs
PTdf['KEGG_Pathways'] = uniprotKeggPathways
PTdf['KEGG_Pathway_Modules'] = uniprotKeggModules
PTdf.drop(columns = ['index'], inplace = True)

#%% load in the statistical data
generalDataLocation = os.path.join(baseFolder, 'data',
                                   'proteomics', '20230620', 'data_analysis',
                                   'PRISM')
subFolders = ['singleFam', 'ALE26_only', 'H16_only']

# change names of columns for table from original file names
changeDict = {'both':{'Group_A':'H16_80_full',
                      'Group_B':'H16_80_lim',
                      'Group_C':'H16_150_lim',
                      'Group_D':'ALE26_80_full',
                      'Group_E':'ALE26_80_lim',
                      'Group_F':'ALE26_150_lim'},
              'ALE26_only':{'Group_A':'ALE26_80_full',
                            'Group_B':'ALE26_80_lim',
                            'Group_C':'ALE26_150_lim'},
              'H16_only':{'Group_A':'H16_80_full',
                          'Group_B':'H16_80_lim',
                          'Group_C':'H16_150_lim'}}

stats = {}
for subfolder in subFolders:
    if 'Fam' in subfolder:
        dictKey = 'both'
    else:
        dictKey = subfolder
    stats[dictKey] = {}
    fullFolder = os.path.join(generalDataLocation, subfolder)
    for fileName in os.listdir(fullFolder):
        if fileName.endswith('.csv'):
            df = pd.read_csv(os.path.join(fullFolder, fileName), sep="\t")
            subDictKey = fileName.replace(subfolder + '_', '')
            subDictKey = subDictKey.replace('singleFam_', '')
            subDictKey = subDictKey.replace('.csv', '')
            for key in changeDict[dictKey]:
                subDictKey = subDictKey.replace(key, changeDict[dictKey][key])
            subDictKey = subDictKey.replace('_vs_', '-')
            df.rename(columns={'ProtID':'UniProtIDs'}, inplace=True)
            df.rename(columns={'log(fold change)':subDictKey + ' log(fold change)'}, inplace=True)
            df.rename(columns={'-log(adjusted p)':subDictKey + ' -log(adjusted p)'}, inplace=True)
            stats[dictKey][subDictKey] = df

# concatenate dataframes first before merging into PTdf
stats['merged'] = {}
for key in stats['both']:
    if key in stats['H16_only'].keys():
        stats['merged'][key] = pd.concat([stats['both'][key], stats['H16_only'][key]])
    elif key in stats['ALE26_only'].keys():
        stats['merged'][key] = pd.concat([stats['both'][key], stats['ALE26_only'][key]])
    else:
        stats['merged'][key] = stats['both'][key]

# merge with PTdf
PTdfmerged = PTdf.copy()
PTdfboth = PTdf.copy()
PTdfh16_only = PTdf.copy()
PTdfale26_only = PTdf.copy()

for key in stats['merged']:
    PTdfmerged = PTdfmerged.merge(stats['merged'][key], how = 'left', on = 'UniProtIDs')
for key in stats['both']:
    PTdfboth = PTdfboth.merge(stats['both'][key], how = 'left', on = 'UniProtIDs')
for key in stats['ALE26_only']:
    PTdfale26_only = PTdfale26_only.merge(stats['ALE26_only'][key], how = 'left', on = 'UniProtIDs')
for key in stats['H16_only']:
    PTdfh16_only = PTdfh16_only.merge(stats['H16_only'][key], how = 'left', on = 'UniProtIDs')

#%% save the dataframes
saveDataLocation = os.path.join(baseFolder, 'data', 'proteomics',
                                '20230620', 'data_analysis', 'poolTablePlots')

PTdfmerged.to_csv(os.path.join(saveDataLocation, 'poolTable_df_merged.csv'), 
                  sep='\t', index = False)
PTdfboth.to_csv(os.path.join(saveDataLocation, 'poolTable_df_expressed_in_both.csv'), 
                  sep='\t', index = False)
PTdfh16_only.to_csv(os.path.join(saveDataLocation, 'poolTable_df_H16_only.csv'), 
                  sep='\t', index = False)
PTdfale26_only.to_csv(os.path.join(saveDataLocation, 'poolTable_df_ALE26_only.csv'), 
                  sep='\t', index = False)