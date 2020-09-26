#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 13:42:10 2020

@author: dmattox
"""

import os, glob, dill#, time, collections
# import itertools
import numpy as np
# import Bio.PDB
import scipy.spatial

oWD = os.getcwd()
oWD += '/'

import lec_gly as LecGly
from bSiteResiFeatures import plipFile
from pocketAnalysis.generateVoxels import mol2Dir, GRID_POINT_SPACING #, POCKET_THRESHOLDS

os.chdir(LecGly.homeDir)

############################################

BURIED_THRESHOLD = 23

# lineSearchResolution = 0.25 # How far apart to space the points on the line when looking for neighboring pocket points

############################################

pcktSurfDir = './data/structures/bsites/pocketSurfaces/'
if not os.path.exists(pcktSurfDir):
    os.makedirs(pcktSurfDir)

############################################


if __name__ == '__main__':
    BURIED_THRESHOLD += 1 # Add one to account for the point being visited being returned in the query list

    
    with open(plipFile, 'rb') as pickleFH:
        allPLIP = dill.load(pickleFH)

    ############################################
    
    # threshColors = ['magenta', 'red', 'orange', 'yellow']
    # colDict = {}
    # if len(threshColors) == len(POCKET_THRESHOLDS):
    #     colDict = {t : c for t,c in zip(POCKET_THRESHOLDS, threshColors)}
    # else:
    #     print("warning, more pocket thresholds than expected")
    
    runList = [] # List to hold subset of PDB ids to run in batch processing
    files = os.listdir(oWD)
    batchList = [f for f in files if f[:7] == 'pdbList'][0]
    with open(oWD + batchList, 'r') as inFH:
        for line in inFH:
            runList.append(line.strip())
            
    # runLst = list(allPLIP.keys())
    
    for pCnt,pdb in enumerate(runList):
        if pCnt % (round(len(runList)/10, 0)) == 0: print("Processing PBD file number " + str(pCnt) + ' of ' + str(len(runList)))
        
        for bs in allPLIP[pdb].bsites:
            print(pdb, bs)
            
            mol2Files = glob.glob(mol2Dir + pdb + '_' + bs.replace(':','-') + '_*')
            
            if not mol2Files:
                print('WARNING: Skipping',pdb,bs + ', no mol2 files for pocket found')
            else:
                for f in mol2Files:
                    pcktPnts = LecGly.getMol2Pnts(f)
                    binThresh = f.split("_")[-1]
                    binThresh = int(binThresh.split('AngBin')[0])
                    
                    if not pcktPnts:
                        print('WARNING: No volume detected in binding pocket, skipping', pdb, bs, binThresh, 'Ang pocket')
                    else:
                        pcktTree = scipy.spatial.KDTree(pcktPnts) # Implement a kd-tree containing all of the points in the pocket
                        
                        searchRad = np.sqrt(3/4) * GRID_POINT_SPACING * 2 + 0.01
                        surfacePnts = []
                    
                        for i,pnt in enumerate(pcktPnts):
                            nearPnts = pcktTree.query_ball_point(pnt, r =searchRad) # Count how many points are direct neighbors of the point
                            if len(nearPnts) < BURIED_THRESHOLD:
                                surfacePnts.append(pnt)
                    
                        # print(binThresh)
                        print('\t', binThresh, 'Ang thresh:', str(round((len(surfacePnts) * 100)/len(pcktPnts), 2)) + '% surface points')
                        # print(len(buriedPnts))
                        
                        LecGly.coords2Mol2(pcktSurfDir + pdb + '_' + bs.replace(':','-') + '_' + str(binThresh) + '_surfacepnts.mol2', surfacePnts)
                    
                
                
                
                
                
                
                
                