#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 18 11:49:17 2020

@author: dmattox
"""

import os
# import collections
import glob
import time
# import scipy.spatial
import numpy as np

oWD = os.getcwd()
oWD += '/'

import lec_gly as LecGly
from generateVoxels import mol2Dir, GRID_POINT_SPACING

os.chdir(LecGly.homeDir)

#################

RESOLUTION = 64
ORDER = 20

#################

voxDir = './data/structures/bsites/zernVoxels/'
if not os.path.exists(voxDir):
    os.makedirs(voxDir)
zernDir = './data/structures/bsites/zernPolys/'
if not os.path.exists(zernDir):
    os.makedirs(zernDir)

#################

if __name__ == '__main__':
    
    allPcktFiles = glob.glob(mol2Dir+ '*.mol2')
    
    gridPnts = [] # initialize grid pnts 
    for x in range(RESOLUTION):
        for y in range(RESOLUTION):
            for z in range(RESOLUTION):
                gridPnts.append([x,y,z])
    # gridTree = scipy.spatial.KDTree(np.array(gridPnts)) # KDTree for grid
    
    
    runList = [] # List to hold subset of PDB ids to run in batch processing
    files = os.listdir(oWD)
    batchList = [f for f in files if f[:7] == 'pdbList'][0]
    with open(oWD + batchList, 'r') as inFH:
        for line in inFH:
            runList.append(line.strip())
            
    start = time.time()
    for pdb in runList:
        # fullPcktMol2Files = glob.glob(pocketDir + pdb + '_' + bs.replace(':','-') + '_*')
        # f = fullPcktMol2Files[1]
        
        pStart = time.time()
        print('##########\n' + pdb + '\n##########\n')
        
        pcktFiles = [f for f in allPcktFiles if pdb in f]
        
        for f in pcktFiles:
    
            bs = f.split("_")[-2].replace('-',':')
            print(bs)
            
            binThresh = f.split("_")[-1]
            binThresh = int(binThresh.split('AngBin')[0])
        
            pcktPnts = LecGly.getMol2Pnts(f)
            
            if len(pcktPnts) <= 1:
                print('WARNING: No volume detected in binding pocket, skipping ' + pdb + ' ' + bs + ' ' + str( binThresh) + ' Ang pocket')
            else:
                pcktPnts = np.array(pcktPnts)
                # pcktTree = scipy.spatial.KDTree(pcktPnts)

                pcktCent = LecGly.getArrayCentroid(pcktPnts)
                
                gridCenter = np.array([0,0,0])
                transV = gridCenter - pcktCent  # vector to translate points s.t. centroid is at the origin
                pcktPnts += transV # translate points to center at origin
                
                pcktCent = LecGly.getArrayCentroid(pcktPnts)
                # print('pseudoatom newCent, pos=[' + ','.join([str(n) for n in pcktCent]) + ']')
                
                maxDist_from_cent = max([ LecGly.eucDist(pcktCent, p) for p in pcktPnts ])
                scaleFactor = (RESOLUTION*0.6/2.0)/maxDist_from_cent # Scale factor, 1 Ang == scaleFactor voxel units
                pVDW = scaleFactor * GRID_POINT_SPACING/1.5 # pseudoVDW radius for pocket points, scaled
                
                extmult = 0.6/maxDist_from_cent # to scale points between [0,0.6], unit sphere has radius 1 and want to stay within 60% of exterior
                pcktPnts *= extmult
                        
                for i,pnt in enumerate(pcktPnts):
                    pcktPnts[i] = 0.5*(pnt+1)*RESOLUTION # translate points to the center of the grid and scale out to fill 60% of nxnxn grid
                    
                q = [list(map(round, p)) for p in pcktPnts]
                
                vxFile = voxDir + pdb + '_' + bs.replace(':','-') + '_' + str(binThresh) + '.vox'
                with open(vxFile, 'w') as outFH:
                    outFH.write(' '.join([str(RESOLUTION)]*3) + '\n')
                    for pnt in gridPnts:
                        if pnt in q:
                            outFH.write('1\n')
                        else:
                            outFH.write('0\n')
                outFile = zernDir + pdb + '_' + bs.replace(':','-') + '_' + str(binThresh) + '_' + str(ORDER)
                os.system('./scripts/pocketAnalysis/vox ' + ' '.join([vxFile, str(ORDER), outFile, str(0) ]))  # Runs C++ binary to generate invariant Zernike descriptors of the nth (ORDER) order, final option 0 indicates to not reconstruct the original shape based on the Zernike polynomials
                
        pEnd = time.time()
        print(str(round(pEnd - pStart, 2)) + ' seconds for PDB file ' + pdb)
    
    end = time.time()
    print(str(round(end - start, 2)/60.) + ' minutes for ' + str(len(runList)) + ' PDB files')
                
                
                
                