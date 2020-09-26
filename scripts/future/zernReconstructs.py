#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 24 10:08:14 2020

@author: dmattox
"""

import os, collections, glob, time
import numpy as np
import scipy.spatial

import Zernike

np.random.seed(27)

def getMol2Pnts(mol2FH):
    ''' Reads in the mol2 file specified by the argument and returns a list of the cooridantes of all atoms within that mol2 file '''
    out = []
    with open(mol2FH, 'r') as inFH:
        strt = False
        for line in inFH:
            if strt == False: # Stil looking for the start of the points
                if line.strip() == '@<TRIPOS>ATOM':
                    strt = True
                    continue
            else:
                line = line.split('\t')
                out.append([float(c) for c in line[2:5]])
    return out

def getCentroid(atmArr):
    '''Given an array of coordiantes, returns the centroid of their coordinates as an array'''
    out = np.array([0,0,0], dtype = 'float32') # Holds the centroid
    for a in atmArr:
        out += a
    out = out/len(atmArr)
    return out

def eucDist(coord1, coord2):
    ''' Calculate the euclidean distance between a pair of 3D coordinates in separate lists '''
    return np.sqrt((coord1[0]-coord2[0])**2 + (coord1[1] - coord2[1])**2 + (coord1[2] - coord2[2])**2)

#################

res = 64
momNum = 12

clusterRun = False 

if clusterRun:
    path = os.getcwd()
    if path[-1] != '/': path += '/'
    
    pocketDir = '/dartfs-hpc/rc/home/y/f002tsy/cbklab/Mattox/glycans/unilec3d/structures/bSites/bsitePockets/'
    outFile = path + '3DZD_' + str(momNum) + 'ord.csv'
    momentDir = '/dartfs-hpc/rc/home/y/f002tsy/cbklab/Mattox/glycans/unilec3d/voxels/moments'+ str(momNum) +'/'
    structDir = '/dartfs-hpc/rc/home/y/f002tsy/cbklab/Mattox/glycans/unilec3d/structures/'
else:
    pocketDir = '/Users/dmattox/cbk/glycan_binding/data/unilectin/structures/bSites/bsitePockets/'
    outFile = '/Users/dmattox/cbk/glycan_binding/analysis/prelim/prelim3/3DZD_test'+ str(momNum) +'.csv'
    momentDir = '/Users/dmattox/cbk/glycan_binding/analysis/prelim/prelim3/zernikeMoments'+ str(momNum) +'/'
    structDir = '/Users/dmattox/cbk/glycan_binding/data/unilectin/structures/'

if not os.path.exists(momentDir):
    os.makedirs(momentDir)

#################

# pdb = '2CL8'
# bs = 'BGC:A:1247'

allPcktFiles = glob.glob(pocketDir+ '*.mol2')



gridPnts = [] # initialize grid pnts and KDTree for grid
for x in xrange(res):
    for y in xrange(res):
        for z in xrange(res):
            gridPnts.append([x,y,z])
gridTree = scipy.spatial.KDTree(np.array(gridPnts))

invars = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(dict))) # 4 level nested dict to hold invariant moments for each [pdb][bs][thresh][(n,l)]
# allZernObjs = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(dict))) # Holds all the Zernike objects so they don't write over each other

 
runList = ['3LL2','2CL8','4YWA']

runList = ['3LL2']

start = time.time()
for pdb in runList:
    # fullPcktMol2Files = glob.glob(pocketDir + pdb + '_' + bs.replace(':','-') + '_*')
    # f = fullPcktMol2Files[1]
    
    pStart = time.time()
    print '##########\n' + pdb + '\n##########\n'
    
    pcktFiles = [f for f in allPcktFiles if pdb in f]
    
    for f in pcktFiles:

        bs = f.split("_")[-2].replace('-',':')
        
        binThresh = f.split("_")[-1]
        binThresh = int(binThresh.split('AngBin')[0])
        
        if binThresh == 6:
            pcktPnts = getMol2Pnts(f)
            
            if len(pcktPnts) <= 1:
                print 'WARNING: No volume detected in binding pocket, skipping ' + pdb + ' ' + bs + ' ' + str( binThresh) + ' Ang pocket'
            else:
                pcktPnts = np.array(pcktPnts)
                pcktTree = scipy.spatial.KDTree(pcktPnts)
                
                randPnts = [pcktPnts[i] for i in np.random.choice(range(len(pcktPnts)), size = 5)] # Pick 5 random points out of the grid use in finding the grid spacing used to construct this pocket
                gridSpacing = pcktTree.query(x = randPnts, k =3) # Find the 3 nearest neighbors for each of these points (will include the original point)
                gridSpacing = np.min(gridSpacing[0][gridSpacing[0] != 0]) # Find the closes distance between any neighbors for any of the 5 random points. All of these points should have direct neighbors on the grid and the smallest non-zero value should be the gridspacing
                if round(gridSpacing, 4) != round((0.5)**(1./3), 4):
                    print("Warning: inconsistent grid spacing detected in PDB " + pdb)
                
                pcktCent = getCentroid(pcktPnts)
                
                gridCenter = np.array([0,0,0])
                transV = gridCenter - pcktCent  # translate points s.t. centroid is at the origin
                pcktPnts += transV
                
                pcktCent = getCentroid(pcktPnts)
                # print('pseudoatom newCent, pos=[' + ','.join([str(n) for n in pcktCent]) + ']')
                
                maxDist_from_cent = max([ eucDist(pcktCent, p) for p in pcktPnts ])
                scaleFactor = (res*0.6/2.0)/maxDist_from_cent # Scale factor, 1 Ang == scaleFactor voxel units
                pVDW = scaleFactor * gridSpacing/1.5 # pseudoVDW radius for pocket points, scaled
                
                extmult = 0.6/maxDist_from_cent # to scale points between [0,0.6], unit sphere has radius 1 and want to stay within 60% of exterior
                pcktPnts *= extmult
                        
                for i,pnt in enumerate(pcktPnts):
                    pcktPnts[i] = 0.5*(pnt+1)*res # translate points to the center of the grid
                pcktCent = getCentroid(pcktPnts)
                # print('pseudoatom newCent, pos=[' + ','.join([str(n) for n in pcktCent]) + ']')
                
                pcktTree = scipy.spatial.KDTree(pcktPnts)
                #################
                
                reload(Zernike) # Need to reload iteratively because Zernike module does not clear existing variables on initializing (deepcopy doesn't work, neither does importing repeatedly)
                
                vox = Zernike.Voxels()
                vox.SetResolution(res)
                
                q = pcktTree.query_ball_tree(gridTree, r = pVDW)
                
                hits = []
                for l in q:
                    for i in l:
                        if i not in hits:
                            hits.append(i)
                            x,y,z = gridPnts[i]
                            vox.SetVoxel(x,y,z,1.0)
                
                vox.Grid2DX(structDir + 'vox_' + pdb + '_' + bs.replace(':','-') + '_ord' + str(momNum) + '.dx')
                        
                zern = Zernike.Zernike()
                
                moments = zern.CalculateMoments(vox,momNum)
                
                moments.SaveMoments(momentDir + pdb + '_' + bs.replace(':','-') + '.mom')
                
                moments.CalcInvariants()
                
                out = [m for m in moments.invariants if m[2] != 0]
                
                for n,l,v in out:
                    invars[pdb][bs][binThresh][(n,l)] = v
                
                # zern.InitialiseReconstruction(moments, momNum, 32)
                # test = zern.ReconstructAll()
                
                # test.Grid2DX(structDir + 'recon_' + pdb  + '_' + bs.replace(':','-') + '_ord' + str(momNum) + '.dx')
        
    pEnd = time.time()
    print str(round(pEnd - pStart, 2)) + ' seconds for PDB file ' + pdb
    
end = time.time()
print str(round(end - start, 2)/60.) + ' minutes for ' + str(len(runList)) + ' PDB files'
    

############
reload(Zernike) 
new_mom = Zernike.Moments()
new_mom.LoadMoments(momentDir + pdb + '_' + bs.replace(':','-') + '.mom')

new_zern = Zernike.Zernike()
new_zern.InitialiseReconstruction(new_mom, momNum, 32)
new_zern.ReconstructAll().Grid2DX(structDir + 'recon_newTest.dx')
############





threshLst = [4,6,8,10]

allTermKeys = set()
for pdb in invars.keys():
    for bs in invars[pdb].keys():
        for thresh in invars[pdb][bs].keys():
            for p in invars[pdb][bs][thresh].keys():
                allTermKeys.add(p)
allTermKeys = sorted(list(allTermKeys))
# print(len(allTermKeys))


with open(outFile,'w') as outFH:
    outFH.write('bsite')
    for k in threshLst:
        outFH.write(',' +','.join(['F'+str(n)+ '.' +str(l)+'_'+str(k)+'Ang' for n,l in allTermKeys]))
    outFH.write('\n')
    for pdb in invars.keys():
        for bs in invars[pdb].keys():
            outFH.write(pdb+ '_' + bs.replace('-',':'))
            for thresh in threshLst:
                if thresh in invars[pdb][bs].keys():
                    outFH.write(',' + ','.join([ str(invars[pdb][bs][thresh][p]) for p in allTermKeys ]) )
                else:
                    outFH.write(',' + ','.join([ 'NA' for p in allTermKeys ]) )
            outFH.write('\n')
    




        
    
    
    











