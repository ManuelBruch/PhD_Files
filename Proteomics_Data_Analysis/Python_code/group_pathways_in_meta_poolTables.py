# -*- coding: utf-8 -*-
"""
Created on Fri Aug 25 11:34:09 2023

@author: Manuel Bruch
relative file paths corrected for all OS
"""

#%% clean up
from IPython import get_ipython
get_ipython().magic('reset -sf') # clears all variables

#%% import statements
import os
import pandas as pd
import json

#%% function definitions
def findOccurrences(s, ch):
    return [i for i, letter in enumerate(s) if letter == ch]


def unpack_brite(brite):
    briteClassMets = {}
    for outerItem in brite['children']:
        briteMetabolism = outerItem['children']
        briteClassMets[outerItem['name'][6:]] = {}
        for item in briteMetabolism:
            if 'regular maps' in item['name']:
                dictKey = 'not assigned'
            else:
                dictKey = item['name'][6:]
            pathwayList = []
            for nestedItem in item['children']:
                pathwayName = nestedItem['name'][6:]
                endPos = findOccurrences(pathwayName, '[')
                try:
                    pathwayList.append(pathwayName[:endPos[0]-1])
                except:
                    pathwayList.append(pathwayName)
            briteClassMets[outerItem['name'][6:]][dictKey] = pathwayList
    return briteClassMets


#%% read data frame with pathway and expression data
curFolder = os.getcwd()
backSlash = findOccurrences(curFolder, os.sep)
baseFolder = curFolder[:backSlash[-1]]

csvFileLocation = os.path.join(baseFolder, 'data', 'proteomics', 
                            '20230620', 'data_analysis', 'poolTablePlots')
# available files:
    # poolTable_df_merged.csv
    # poolTable_df_expressed_in_both.csv
    # poolTable_df_H16_only.csv
    # poolTable_df_ALE26_only.csv
dataFile = 'poolTable_df_merged.csv'

dataTable = pd.read_csv(os.path.join(csvFileLocation, dataFile), sep="\t")

#%% get list of all KEGG pathways and all UniProt pathways
# create sets to get all the unique identifiers and remove nan
keggPws = list(dataTable['KEGG_Pathways'])                                      # 2885 items, 1150 of which are not nan
keggSet = list(set(keggPws))
keggSet = [x for x in keggSet if str(x) != 'nan']                               # 330 items

uniProtPws = list(dataTable['UniProt_Pathways'])                                # 2885 items, 354 of which are not nan
uniProtSet = list(set(uniProtPws))
uniProtSet = [x for x in uniProtSet if str(x) != 'nan']                         # 311 items

# split the strings relating to several pathways in the two sets
flattenedKegg = [x.split('//') for x in keggSet]
flattenedKegg = [item for l in flattenedKegg for item in l]                     # 996 items
flattenedKeggSet = list(set(flattenedKegg))                                     # 122 items

flattenedUniProt = [x.split('//') for x in uniProtSet]
flattenedUniProt = [item for l in flattenedUniProt for item in l]               # 311 items


#%% automatic mapping via json file of KEGG Brite hierarchy to classes
jsonPath = os.path.join(csvFileLocation, 'reh00001.json')
with open(jsonPath) as jsonFile:
    brite = json.load(jsonFile)

# unpack briteMetabolism in usable dictionary {'class1': ['pathway1',,...],...}
briteClassMets = unpack_brite(brite)

# manually add the 7 upcoming pathways into the Brite hierarchy
# do not add in 'Carbon metabolism' or 'Degradation of aromatic compounds' --> should be kept as is in the flattened KeggSet
briteClassMets['Metabolism']['Amino acid metabolism'].append('Biosynthesis of amino acids')
briteClassMets['Metabolism']['Amino acid metabolism'].append('2-Oxocarboxylic acid metabolism')
briteClassMets['Metabolism']['Metabolism of cofactors and vitamins'].append('Biosynthesis of cofactors')
briteClassMets['Metabolism']['Lipid metabolism'].append('Fatty acid metabolism')
briteClassMets['Metabolism']['Nucleotide metabolism'].append('Biosynthesis of nucleotide sugars')

checkPathways = False
if checkPathways:
    # check if all pathways in the flattenedKeggSet are present in the dict
    # flatten dict to list
    realPWlist = list(briteClassMets.values())
    realPWlist = [list(x.values()) for x in realPWlist]
    realPWlist = [x for y in realPWlist for x in y]
    realPWlist = [x for y in realPWlist for x in y]
    
    controlList = []
    for item in flattenedKeggSet:
        if item in realPWlist:
            controlList.append(True)
        else:
            controlList.append(False)
    
    falseIndex = []
    for i in range(len(controlList)):
        if not controlList[i]:
            falseIndex.append(flattenedKeggSet[i])

#%% construct an assignment dictionary from the flattenedKeggSet
keggPWdict = {}
for pw in flattenedKeggSet:
    keggPWdict[pw] = pw
    pwFound = False
    for pwClass in briteClassMets:
        for pwGroup in briteClassMets[pwClass]:
            if pw in briteClassMets[pwClass][pwGroup]:
                keggPWdict[pw] = pwGroup
                pwFound = True
                continue
        if pwFound:
            continue
# this should create 24 groups, including the two non-replaced pathways (see line 88)

#%% make a list for the poolTablePlot table with the respective pathways
keggMetaList = []
for item in keggPws:
    if str(item) == 'nan':
        keggMetaList.append('not annotated')
    else:
        allPws = item.split('//') # will always return a list
        for i in range(len(allPws)):
            allPws[i] = keggPWdict[allPws[i]]
        allPws = set(allPws)
        allPws = '//'.join(allPws)
        keggMetaList.append(allPws)


#%% compare KEGG ant UniProt pathways
diffDict = {}
for i in range(len(uniProtPws)):
    if not isinstance(keggPws[i], str):
        if isinstance(uniProtPws[i], str):
            diffDict[i] = uniProtPws[i]
# there are 26 items in the uniprot pathways that are not annotated in KEGG
# manually replace these
diffDict[13] = 'Carbohydrate metabolism'
diffDict[14] = 'Carbohydrate metabolism'
diffDict[19] = 'Folding, sorting and degradation'
diffDict[248] = 'Glycan biosynthesis and metabolism'
diffDict[283] = 'Biosynthesis of cofactors'
diffDict[389] = 'Glycan biosynthesis and metabolism'
diffDict[723] = 'Fatty acid metabolism'
diffDict[880] = 'Biosynthesis of cofactors'
diffDict[915] = 'Cell growth and death'
diffDict[967] = 'Cell growth and death'
diffDict[1010] = 'Transcription'
diffDict[1011] = 'Transcription'
diffDict[1134] = 'Energy metabolism'
diffDict[1167] = 'Transcription'
diffDict[1335] = 'Fatty acid metabolism'
diffDict[1348] = 'Translation'
diffDict[1738] = 'Amino acid metabolism'
diffDict[2019] = 'Cell growth and death'
diffDict[2163] = 'Cell growth and death'
diffDict[2313] = 'Biosynthesis of cofactors'
diffDict[2336] = 'Transcription'
diffDict[2486] = 'Translation'
diffDict[2497] = 'Translation'
diffDict[2500] = 'Translation'
diffDict[2571] = 'Biosynthesis of cofactors'
diffDict[2862] = 'Cell growth and death'

for key in diffDict:
    keggMetaList[key] = diffDict[key]

dataTable['KEGG_and_UniProt_Meta_Pathways'] = keggMetaList # 1709 unannotated proteins

#%% save the newly created table
saveFile = 'poolTable_df_merged_metaPWs.csv'
dataTable.to_csv(os.path.join(csvFileLocation, saveFile), 
                  sep='\t', index = False)