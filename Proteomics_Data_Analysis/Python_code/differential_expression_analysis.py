# -*- coding: utf-8 -*-
"""
Created on Tue Jul 11 15:01:31 2023

@author: Manuel Bruch
relative filepaths corrected for all OS
"""
# Analysis of proteomics data for strain dependent expression
# importet data tables created from MaxQuant analysis run by Caitriona Scaife
# main purpose: create list of proteins that are only expressed in C. necator 
# H16 orALE26

#%% clean up
from IPython import get_ipython
get_ipython().magic('reset -sf') # clears all variables

#%% import statements
import os
import numpy as np
import pandas as pd
# import requests, sys, json

#%% function definitions
def findOccurrences(s, ch):
    return [i for i, letter in enumerate(s) if letter == ch]

def parse_fasta(fasta):
    proteinDic = {}
    for line in fasta:
        if '>' in line:
            uniqueID, entryName, proteinName, OrgName, OrgID, GeneName, protExist, SeqVer = parse_fasta_header(line)
            proteinDic[uniqueID] = {}
            proteinDic[uniqueID]['EntryName'] = entryName
            proteinDic[uniqueID]['ProteinName'] = proteinName
            proteinDic[uniqueID]['OrganismName'] = OrgName
            proteinDic[uniqueID]['OrganismID'] = OrgID
            proteinDic[uniqueID]['GeneName'] = GeneName
            proteinDic[uniqueID]['ProteinExistence'] = protExist
            proteinDic[uniqueID]['SequenceVersion'] = SeqVer
    return proteinDic


def parse_fasta_header(headerLine):
    # find positions to split the string at
    IDsep = findOccurrences(headerLine, '|')
    spaces = findOccurrences(headerLine, ' ')
    OrgNameStart = str.find(headerLine, 'OS=')
    OrgIDStart = str.find(headerLine, 'OX=')
    GeneNameStart = str.find(headerLine, 'GN=')
    PEStart = str.find(headerLine, 'PE=')
    SeqVerStart = str.find(headerLine, 'SV=')
    # write parts into variables
    uniqueID = headerLine[IDsep[0]+1:IDsep[1]]
    entryName = headerLine[IDsep[1]+1:spaces[0]]
    proteinName = headerLine[spaces[0]+1:OrgNameStart-1]
    OrgName = headerLine[OrgNameStart+3:OrgIDStart-1]
    OrgID = headerLine[OrgIDStart+3:GeneNameStart-1]
    GeneName = headerLine[GeneNameStart+3:PEStart-1]
    protExist = headerLine[PEStart+3:SeqVerStart-1]
    SeqVer = headerLine[SeqVerStart+3:-1]
    return uniqueID, entryName, proteinName, OrgName, OrgID, GeneName, protExist, SeqVer


def tidy_up_Df(df):
    for index, row in df.iterrows():
        firstColEntry = row.iloc[0]
        try:
            if '}' in firstColEntry:
                altVal = firstColEntry[firstColEntry.find('}')+1:]
                df.iat[index,0] = altVal
        except:
            continue
    return df


def extend_df_from_dict(df, myDict, colName, dictValues, omitRows = []):
    searchList = df[colName].tolist()
    if omitRows:
        del searchList[omitRows[0]:omitRows[1]+1]
    allValues = []
    for val in dictValues:
        # curList = [myDict[x] for x in searchList] # error because some have multiple identifiers matched to them
        # as for loops rather than list comprehension because my head hurts
        curList = []
        for item in searchList:
            if ';' in item:
                thisItem = item.split(';')
                thisItem = [x for x in thisItem if ('REV') not in x] # there still might be a reverse sequence discovered --> delete item before list comprehension
                curList2 = [myDict[x][val] for x in thisItem]
                curList.append(curList2)
            else:
                curList.append(myDict[item][val])
        allValues.append(curList)
    # create nan lists to append at the start of the created list based on the number of omitted rows
    nanList = np.empty(omitRows[1]+1)
    nanList[:] = np.nan # I don't think that does anything. np.empty() already creates nan list?
    nanList = list(nanList)
    for i in range(len(allValues)):
        allValues[i] = nanList + allValues[i]
    # add new lists to original dataframe
    for i in range(len(dictValues)):
        df[dictValues[i]] = allValues[i]
    return(df)


# entries in proteomics data might match to multiple uniprot identifiers.
# to search NCBI using Bio.Entrez, a list of identifiers can be provided but I
# am not sure that a nested list would work. Therefore the original list is to 
# be split into individual lists for each "layer" of detected protein
# --> the first list should contain a valid entry for each protein
# --> the second and further lists will contain only entries, if a given protein
#     has been found to have multiple possible DB entries
# the positions of these duplicates in the original list is stored in the
# multiEntryIndeces variable which is a list of lists that has the length
# len(allTheLists) - 1 as in the last loop iteration there are no further ";"
# detected
def ncbi_list_prep(originalList):
    searchList = originalList
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


def get_NCBI_data(protAcc):
    handle = Entrez.efetch(db="protein", id=protAcc, rettype="gb", retmode="text") # xml or text
    plainEntryText = handle.read()
    handle.close()
    return plainEntryText


def map_gene_location(ncbiTextFile):
    individualAccs = ncbiTextFile.split("accession ") # now the first 6 characters should be the accession number
    individualAccs = [x[0:6] for x in individualAccs[1:]]
    xrefNumbers = ncbiTextFile.split("xrefs: ") # first line should be xref numbers
    xrefNumbers = xrefNumbers[1:]
    geneLocations = dict()
    for i in range(len(xrefNumbers)):
        refEnd = str.find(xrefNumbers[i], "xrefs (non-sequence databases)")
        curRef = xrefNumbers[i][:refEnd]
        if "AY305378.1" in curRef:
            geneLocations[individualAccs[i]] = "pHG1"
        elif "AM260479.1" in curRef:
            geneLocations[individualAccs[i]] = "chr1"
        elif "AM260480.1" in curRef:
            geneLocations[individualAccs[i]] = "chr2"
    return geneLocations


#%% parse the fasta file with protein names
# workaround for relative filepath:
curFolder = os.getcwd()
backSlash = findOccurrences(curFolder, os.sep)
baseFolder = curFolder[:backSlash[-1]]
fileLocation = os.path.join(baseFolder, 'data', 'proteomics', 
                            '20230620', 'data_analysis')
fileName = 'cupriavidus_necator_fasta_includeIsoform_true_query__28prot-2023_06_23-12_52_17_72.fasta'
fastaFile = os.path.join(fileLocation, fileName)

with open(fastaFile, 'r') as fasta:
    # proteinData = parse_fasta(fasta)
    proteinDic = parse_fasta(fasta)

#%% read in proteomics table for proteins only expressed in either ALE or H16 strian
# fileName = 'export_proteins_only_exp_in_one_strain.txt'
# dataFile = os.path.join(fileLocation, fileName)
# protData = pd.read_csv(dataFile, sep="\t")

# # after export from Perseus, first few rows have row names that are pushed in to the first cell of the row
# protData = tidy_up_Df(protData)

# # create individual dataframes for H16 and ALE26
# H16Mask = protData.iloc[1].str.contains('H16', na = False)
# H16df = protData[H16Mask.index[H16Mask]]

# ALE26Mask = protData.iloc[1].str.contains('ALE26', na = False)
# ALE26df = protData[ALE26Mask.index[ALE26Mask]]

# # add column with protein identifiers
# protIDs = protData["Protein IDs"]
# H16df = H16df.assign(proteinIDs = protIDs.values)
# ALE26df = ALE26df.assign(proteinIDs = protIDs.values)

# # remove rows where only the protienIDs are not nan
# H16df = H16df.dropna(thresh=2)
# ALE26df = ALE26df.dropna(thresh=2)

# # add columns with protein names and uniprot entry names
# H16df = extend_df_from_dict(H16df, proteinDic, "proteinIDs", ["EntryName", "ProteinName"], omitRows = [0,4])
# ALE26df = extend_df_from_dict(ALE26df, proteinDic, "proteinIDs", ["EntryName", "ProteinName"], omitRows = [0,4])

#%% read in proteomics table for proteins expressed in both strains
fileName = 'export_proteins_exp_in_both_strains.txt'
dataFile = os.path.join(fileLocation, fileName)
protData = pd.read_csv(dataFile, sep="\t")

# after export from Perseus, first few rows have row names that are pushed in to the first cell of the row
protData = tidy_up_Df(protData)
allData = protData.copy()

# there still might be a reverse sequence discovered --> delete item before list comprehension
protIDs = list(allData['Protein IDs'])
for i in range(len(protIDs)):
    try:
        if 'REV' in protIDs[i]:
            protIDs[i] = protIDs[i][:protIDs[i].find('REV')-1] # only leaves IDs before the Rev id
    except:
        continue
allData['Protein IDs'] = protIDs
            

# add column with protein identifiers
allData = extend_df_from_dict(allData, proteinDic, "Protein IDs", ["EntryName", "ProteinName"], omitRows = [0,4])

#%% Basic UniProt search Setup
# Documentation: https://www.uniprot.org/help/api
# WEBSITE_API = "https://rest.uniprot.org/"

# Documentation: https://www.ebi.as.uk/proteins/api/doc/
# PROTEINS_API = "https://www.ebi.ac.uk/proteins/api"

#%% do some uniprot searching
# helper function to download data
# def get_url(url, **kwargs):
#     response = requests.get(url, **kwargs)
    
#     if not response.ok:
#         print(response.text)
#         response.raise_for_status()
#         sys.exit()
    
#     return response
# get the first protein in the H16 dataframe
# protAcc = H16df["EntryName"].iloc[15]
# if type(protAcc) is list:
#     protAccAll = protAcc
#     protAcc = protAcc[0]
# r = get_url(f"{WEBSITE_API}/uniprotkb/{protAcc}?fields=xref_geneid")
# a = r.json()
# crossRefDict = a['uniProtKBCrossReferences']
# a["genes"][0]["geneName"]["value"]
# a["uniProtKBCrossReferences"][0]["geneName"]["value"]
# got to figure this unpacking out
# print(json.dumps(r.json(), indent=2)) # still prints way too much crap

#%% perhaps can search GenBank with the proteinID directly
from Bio import Entrez
Entrez.email = "manuel.bruch@ucdconnect.ie"

"""
curate list of protein IDs from dataframe
only one ID searchable per line directly
"""
# # Cupriavidus necator H16
# allH16Ids = list(H16df["proteinIDs"].iloc[5:])
# H16searchLists, H16duplicateEntryPositions = ncbi_list_prep(allH16Ids)

# # Cupriavidus necator ALE26
# allALE26Ids = list(ALE26df["proteinIDs"].iloc[5:])
# ALE26searchLists, ALE26duplicateEntryPositions = ncbi_list_prep(allALE26Ids)

# all common proteins
allProtIDs = list(allData["Protein IDs"].iloc[5:])
allProtSearchLists, allProtSduplicateEntryPositions = ncbi_list_prep(allProtIDs)

searchParamDict = {}
# # H16 and ALE seperately
# searchParamDict['H16'] = {}
# searchParamDict['H16']['searchList'] = H16searchLists
# searchParamDict['H16']['entryPos'] = H16duplicateEntryPositions
# searchParamDict['ALE26'] = {}
# searchParamDict['ALE26']['searchList'] = ALE26searchLists
# searchParamDict['ALE26']['entryPos'] = ALE26duplicateEntryPositions

# both together
searchParamDict['allProts'] = {}
searchParamDict['allProts']['searchList'] = allProtSearchLists
searchParamDict['allProts']['entryPos'] = allProtSduplicateEntryPositions

# search the protein database
for key in searchParamDict:
    searchList = searchParamDict[key]['searchList']
    entryPos = searchParamDict[key]['entryPos']
    geneLocations = ["not found"]*len(searchList[0])
    for i in range(len(searchList)):
        protAcc = searchList[i]
        plainEntryText = get_NCBI_data(protAcc)
        genomeLocationDict = map_gene_location(plainEntryText)
        if i == 0:
            for j in range(len(protAcc)):
                try:
                    geneLocations[j] = genomeLocationDict[protAcc[j]]
                except:
                    continue
        else:
            for j in range(len(entryPos[i-1])):
                geneLocations[entryPos[i-1][j]] = geneLocations[entryPos[i-1][j]] + '; ' + genomeLocationDict[protAcc[j]]
    # extend the list by nan values for first few rows
    nanList = np.empty(5)
    nanList = list(nanList)
    searchParamDict[key]['geneLocations'] = nanList + geneLocations

#%% alter and save original dataframes

# # H16 and ALE seperately
# H16df['geneLocations'] = searchParamDict['H16']['geneLocations']
# H16FileName = 'H16_only_proteins.csv'
# H16df.to_csv(os.path.join(fileLocation, H16FileName), sep='\t')
# ALE26df['geneLocations'] = searchParamDict['ALE26']['geneLocations']
# ALE26FileName = 'ALE26_only_proteins.csv'
# ALE26df.to_csv(os.path.join(fileLocation, ALE26FileName), sep='\t')

# for common proteins:
allData['geneLocations'] = searchParamDict['allProts']['geneLocations']
commonFileName = 'common_proteins.csv'
allData.to_csv(os.path.join(fileLocation, commonFileName), sep='\t')