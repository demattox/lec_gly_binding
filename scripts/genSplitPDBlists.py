#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug  1 16:23:55 2020

@author: dmattox
"""

import os
import dill

import lec_gly as LecGly
from bSiteResiFeatures import plipFile

os.chdir(LecGly.homeDir)

##########################
outDir = './data/structures/bsites/batchLists/'
if not os.path.exists(outDir):
    os.makedirs(outDir)
##########################

maxBatchLstLength = 50 # generates 28 lists from 1365 PDB IDs, 27 lists of 50 and 1 of 15

##########################

with open(plipFile, "rb") as pickleFH:
    allPLIP = dill.load(pickleFH)

def chunks(lst, n): # Function borrowed from https://stackoverflow.com/questions/312443/how-do-you-split-a-list-into-evenly-sized-chunks
    """Yield successive n-sized chunks from a list"""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]


splitLst = chunks(list(allPLIP.keys()), maxBatchLstLength)
for i,lst in enumerate(splitLst):
    with open(outDir + 'pdbList_' + str(i) + '.txt', 'w') as outFH:
        for pdb in lst:
            outFH.write(pdb + '\n')
    
    
    
    