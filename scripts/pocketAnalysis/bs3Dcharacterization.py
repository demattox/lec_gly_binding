#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 17 10:44:57 2020

@author: dmattox
"""

import os, glob, dill, time
# import itertools, collections
import numpy as np
import Bio.PDB
# import pliptool.plip.modules.plipxml as plipxml
from scipy.spatial import ConvexHull
from sklearn.cluster import DBSCAN

# from checkPLIPligands import res2plipID, eucDist, getLigAtoms

import sys, math # for running on the cluster
sys.path.insert(0,'/dartfs-hpc/rc/lab/C/CBKlab/Mattox/glycans/unilec3d/voxels')


############################################

def getCentroid(atmLst):
    '''Given a list of BioPDB atoms, returns the centroid of their coordinates as an array'''
    out = np.array([0,0,0], dtype = 'float32') # Holds the centroid
    for a in atmLst:
        out += a.coord
    out = out/len(atmLst)
    return out

def getPntsFromWRL(wrlFH):
    '''Given a list of BioPDB atoms, returns the centroid of their coordinates as an array'''
    out = []
    with open(wrlFH, 'r') as inFH:
        start = False
        for line in inFH:
            line = line.strip().split(' ')
            if start:
                if len(line) != 3:
                    return out # once the formatting changes, done collecting points
                else:
                    line = [float(l.replace(',','')) for l in line]
                    out.append(line)
            elif line[0] == 'point':
                start = True

def coords2Mol2(mol2FH, coordArray):
    '''
    Writes an array of coordinates to a TRIPOS mol2 file with each cooridnate point assigned to an unbonded oxygen atom.
    Saved to file handle indicated in mol2FH argument
    '''
    with open(mol2FH, 'w') as outFH:
        outFH.write('@<TRIPOS>MOLECULE\n')
        outFH.write('grid\n')
        outFH.write(str(len(coordArray)) + ' 0 0 0\t0\n')
        outFH.write('SMALL\n')
        outFH.write('NO_CHARGES\n')
        outFH.write('\n\n@<TRIPOS>ATOM\n')
        for i,pnt in enumerate(coordArray):
            i += 1
            outFH.write(str(i) + '\tO' + str(i) + '\t' )
            outFH.write('\t'.join([str(p) for p in pnt]) + '\t')
            outFH.write('O.2\t1\tGRID\t0.000\n')

# Manually importing fuctions from checkPLIPligands to avoid installing additional dependencies on cluster
def eucDist(coord1, coord2):
    ''' Calculate the euclidean distance between a pair of 3D coordinates in separate lists '''
    return math.sqrt((coord1[0]-coord2[0])**2 + (coord1[1] - coord2[1])**2 + (coord1[2] - coord2[2])**2)

def res2plipID(bioPDBresi):
    ''' Generate a unique id from Bio.PDB residue following the PLIP format, automated residue vs hetatm detection '''
    if bioPDBresi.full_id[3][0][:2].strip() == 'H_':  # If hetatm
        out = bioPDBresi.full_id[3][0][2:].strip() + ':' + bioPDBresi.full_id[2] + ':' + str(bioPDBresi.full_id[3][1])
    else:
        out = bioPDBresi.get_resname() + ':' + bioPDBresi.full_id[2] + ':' + str(bioPDBresi.full_id[3][1])
    return out

def getLigAtoms(bioPDBModel, plipXML_bsite):
    ''' For a given PLIP binding site (eg allPLIP["5D9Z"].bsites["BMA:A:201"]), returns all heavy ligand atoms in the ligand '''
    out = []
    for res in bioPDBModel.get_residues():
        resID = res2plipID(res)
        if resID in plipXML_bsite.members:
            for a in res.get_atoms():
                if a.element != 'H':
                    out.append(a)
    return out

##########


if __name__ == '__main__':
    
    # os.chdir('/Users/dmattox/cbk/glycan_binding/')
    
    ############################################
    
    # pdbDir = './data/unilectin/structures/holo/rawPDBs/' # Downloaded straight from RCSB with wget, original numbering including icodes
    pdbDir = '/dartfs-hpc/rc/home/y/f002tsy/cbklab/Mattox/glycans/unilec3d/structures/holo/' # for running on cluster
    
    # wrlDir = './data/unilectin/structures/bSites/bsiteSurface/' # wrl representation of surface of binding site residues
    wrlDir = '/dartfs-hpc/rc/home/y/f002tsy/cbklab/Mattox/glycans/unilec3d/structures/bSites/bsiteSurface/wrlFiles/' # for running on cluster
    
    # mol2Dir = './data/unilectin/structures/bSites/bsitePockets/'
    mol2Dir = '/dartfs-hpc/rc/home/y/f002tsy/cbklab/Mattox/glycans/unilec3d/structures/bSites/bsitePockets/' # for running on cluster
    
    # plipDat = './data/unilectin/cleanPLIPrecords.p'
    plipDat = '/dartfs-hpc/rc/home/y/f002tsy/cbklab/Mattox/glycans/unilec3d/voxels/cleanPLIPrecords.p'  # for running on cluster
    
    # logFile = './data/unilectin/structures/logFile_3Dbsites.txt'
    logFile = 'logFile_3Dbsites.txt' # for running on cluster
    
    ############################################
    
    minCubeEdge = 20 # define the minimum length of the cube edge
    gridPntSpacing = (0.5)**(1/3) # how far apart to place the grid points
    proAtmThresh = 2.5 # Threshold distance to a heavy protein atom to consider a point within the protein surface
    maxDistFromLig = 10 # Threshold for maximum allowed distance for any grid point from a heavy ligand atom
    binThresholds = [4,6,8] # Distance values for binning points prior to clustering
    binThresholds.append(maxDistFromLig)
    
    
    DBSCANeps = np.sqrt(3/4) * gridPntSpacing * 2 + 0.01 # The maximum distance between two samples for one to be considered as in the neighborhood of the other. Setting to the max length between directly neighboring points, allowing for 26 possible nieghbors 
    DBSCANmin_samps = 13 # The number of samples (or total weight) in a neighborhood for a point to be considered as a core point.
    clusterSizeThresh = 15 # Prune clusters of pocket points with less than 15 Ang**3 of volume
    clusterSizeThresh = clusterSizeThresh/(gridPntSpacing**3) # convert from Ang**2 to number of points based on spacing
    pntsNearLig = 0.1 # Prune clusters with fewer than 10% of its points within 2 Ang of any heavy ligand atom
    
    
    
    ############################################
    
    pdbFiles = glob.glob(pdbDir + '*.pdb')
    
    with open(plipDat, 'rb') as pickleFH:
        allPLIP = dill.load(pickleFH)
    
    # with open("./data/unilectin/icodeResi2PLIPmapping.p", "rb")as pickleFH:
    #     zeroResis = dill.load(pickleFH)
    
    ############################################
    
    globalStart = time.time()
        
    plipCnt = 0
    bsCnt = 0
    
    runList = [] # List to hold subset of PDB ids to run in batch processing
    files = os.listdir()
    batchList = [f for f in files if f[:7] == 'pdbList'][0]
    with open(batchList, 'r') as inFH:
        for line in inFH:
            runList.append(line.strip())
    
    # runList = list(allPLIP.keys())

    with open(logFile, 'w') as logFH:
        for k,pdb in enumerate(runList):
            
            plipCnt +=1
            if k % (round(len(runList)/10, 0)) == 0: print("Processing PBD file number " + str(k) + ' of ' + str(len(runList)))
            logFH.write('\n########\n'+pdb+'\n########\n')
            
            start = time.time()
            
            pdbIN = [f for f in pdbFiles if f[-8:-4] == pdb][0]
            parser = Bio.PDB.PDBParser(QUIET=True)
            model = parser.get_structure(pdb, pdbIN)[0] # grab the first model in PDB file
            
            searcher = Bio.PDB.NeighborSearch([a for a in model.get_atoms() if (a.get_parent().full_id[3][0][:2].strip() == '' and a.element != 'H')]) # Prepare bioPDB searcher to eliminate resiudes within the protein surface (<2.5 Ang from a heavy protein atom)
            
            for bs in allPLIP[pdb].bsites:
                
                bsCnt += 1
                logFH.write('\n' + bs + '\n')
                            
                # Identify ligand and binding site atoms
                ligAtms = getLigAtoms(model, allPLIP[pdb].bsites[bs])
                
                # bsResIDs = [res['aa'] + ':' + res['reschain'] + ':' + str(res['resnr']) for res in allPLIP[pdb].bsites[bs].bs_res]
                # if pdb in zeroResis.keys(): # Check for zero numbered resiudes and correct the res num to match the original number including the i_code
                #     if bs in zeroResis[pdb]:
                #         for i,resID in enumerate(bsResIDs):
                #             resID = resID.split(':')
                #             rID = resID[2]+resID[1] # resnum from PLIP with chain ID
                #             if rID in zeroResis[pdb][bs].keys():
                #                 resID[2] = zeroResis[pdb][bs][rID] # Replace PLIP resnum with true resnum with icode
                #                 bsResIDs[i] = ':'.join(resID)
                # proAtms = [] # Holds all heavy atoms in residues within the binding site
                # for res in model.get_residues():
                #     resID = res2plipID(res)
                #     if res.full_id[3][2] != ' ': # If the residue has an i_code, append it to the resID to match the modified PLIP resIDs
                #         resID = resID + res.full_id[3][2]
                #     if resID in bsResIDs:
                #         for a in res.get_atoms():
                #             if a.element != 'H':
                #                 proAtms.append(a)
                
                CL = getCentroid(ligAtms) # Ligand centroid
                # CP = getCentroid(proAtms) # Binding sites residues centroid
                # cDist = eucDist(CL,CP) # Distance between ligand centroid and binding site residues centroid
                
                ligAtmDists = [eucDist(CL, a.coord) for a in ligAtms] # Distance between ligand centroid and all atoms in the ligand
                cubeEdge = minCubeEdge if max(ligAtmDists) <= (minCubeEdge/2.5) else round(2.5*max(ligAtmDists)) # Set the size of the cube to be at least 20 Ang along each edge, or larger for larger ligands. Round to keep grid points even
                # apprxCL = [round(p) for p in CL] # rounded coordinates of ligand centroid
                
                cubeLimits = [[d + cubeEdge/2, d - cubeEdge/2] for d in CL]
                # cubeVerts = list(itertools.product(*cubeLimits))
                # coords2Mol2('data/unilectin/structures/testcubeverts.mol2', cubeVerts)
                
                gridPnts = [] # hold the grid points within the cube defined above centered on the ligand centroid
                for x in np.arange(min(cubeLimits[0]),max(cubeLimits[0])+gridPntSpacing, step = gridPntSpacing):
                    for y in np.arange(min(cubeLimits[1]),max(cubeLimits[1])+gridPntSpacing, step = gridPntSpacing):
                        for z in np.arange(min(cubeLimits[2]),max(cubeLimits[2])+gridPntSpacing, step = gridPntSpacing):
                            gridPnts.append([x,y,z])
                
                # coords2Mol2('data/unilectin/structures/testgridInit.mol2', gridPnts)                
                
                wrlFile = wrlDir + pdb + '_' + bs.replace(':','-') + '.wrl'
                proCoords = getPntsFromWRL(wrlFile) # Points representing the van der Waals surface of the binding site residues
                
                coords2Mol2('data/unilectin/structures/' + bs+'surface.mol2', proCoords)
                
                ########################################################
                # Begin excluding points          
                      
                # Exclude points within protein surface
                nonBuriedPnts = []
                for p in gridPnts:
                    closeProAtms = searcher.search(p,radius = proAtmThresh, level = 'A')
                    if len(closeProAtms) == 0:
                        nonBuriedPnts.append(p)
                # coords2Mol2('data/unilectin/structures/testgridCloseInhullNotburied.mol2', nonBuriedPnts)
                logFH.write('From ' + str(len(gridPnts)) + ' points, ' + str(len(nonBuriedPnts)) + ' points found outside of the protein surface\n')
                # print('From ' + str(len(gridPnts)) + ' points, ' + str(len(nonBuriedPnts)) + ' points found outside of the protein surface\n')
                

                # Exclude points far from the ligand
                binPnts = [[] for b in binThresholds] # Holds points in each bin for each binning distance threshold lower than the max threshold
                for p in nonBuriedPnts:
                    dists = [eucDist(p,a.coord) for a in ligAtms]
                    if min(dists) <= maxDistFromLig:
                        d = min(dists)
                        for i,b in enumerate(binThresholds):
                            if d <= b:
                                binPnts[i].append(p)
                            
                            
                # coords2Mol2('data/unilectin/structures/testgridClose.mol2', binPnts[-1])
                logFH.write('From ' + str(len(nonBuriedPnts)) + ' points, ' + str(len(binPnts[-1])) + ' points found within the max threshold distance of a ligand atom\n')
                # print('From ' + str(len(nonBuriedPnts)) + ' points, ' + str(len(binPnts[-1])) + ' points found within the threshold distance of a ligand atom\n')

        
                # Exclude points outside of convex hull
                initHull = ConvexHull(proCoords)
                hullVerts = initHull.points[initHull.vertices]
                hull = ConvexHull(points = hullVerts, incremental=True)
                inHullPnts = []
                for p in binPnts[-1]:
                    newVertInd = len(hull.points) # If a new vert shows up with this index, it is outside of the convex hull
                    hull.add_points(np.array([p]))
                    if newVertInd in hull.vertices: # if point is outside convex hull
                        hull = ConvexHull(points = hullVerts, incremental=True) # reset convex hull
                    else:
                        inHullPnts.append(p)
                hull.close()
                # coords2Mol2('./data/unilectin/structures/bSites/pntsb4DBSCAN/' + pdb + '_' + bs.replace(':','-') + '.mol2', inHullPnts)
                logFH.write('From ' + str(len(binPnts[-1])) + ' points, ' + str(len(inHullPnts)) + ' points found within the convex hull\n\n')
                # print('From ' + str(len(binPnts[-1])) + ' points, ' + str(len(inHullPnts)) + ' points found within the convex hull\n')
                
                # Drop points outside of the convex hull from the binned points as well
                for i in range(len(binPnts)):
                    binPnts[i] = [p for p in binPnts[i] if p in inHullPnts]
                
                # for b in binPnts: print(len(b))
                
                # Cluster points together to identify individual pockets and eliminate pockets far from the ligand
                
                prevHit = False # Store whether or not a set of points with a lower distance threshold had a cluster pass the minimum percentage threshold for points within 2 Ang of a heavy ligand atom
                for i,pntSet in enumerate(binPnts):
                    pntArr = np.array(pntSet)
                    clusters = DBSCAN(eps = DBSCANeps, min_samples= DBSCANmin_samps).fit(pntArr) #DBSCAN to cluster non-excluded points together into pockets
                    # print(i, set(clusters.labels_))
                    
                    # for ind in set(clusters.labels_):
                    #     print('Cluster ind',str(ind))
                    #     print('\t',str(sum(clusters.labels_ == ind)), 'points')
                    #     coords2Mol2('data/unilectin/structures/testClust'+ '_' + str(i) + '_' +str(ind) +'.mol2', pntArr[clusters.labels_ == ind])
                
                    goodClusts = []
                    failed_bestInd = ''
                    failed_bestPercentage = 0 # If a pocket has clusters that pass the threshold initially but later fails to hit the 10% (pntsNearLig) threshold as the pocket expands, for each threshold that fails keep the cluster with the most points near the ligand atoms
                    for ind in set(clusters.labels_):
                        if ind != -1: # Don't consider the noise cluster
                            clusPnts =  pntArr[clusters.labels_ == ind]
                            if len(clusPnts) > clusterSizeThresh: # Filter out clusters with too few points
                                closePntCnt = 0
                                for p in clusPnts:
                                    dists = [eucDist(p, a.coord) for a in ligAtms]
                                    if min(dists) <= 2:
                                        closePntCnt += 1
                                if closePntCnt >= len(clusPnts)*pntsNearLig:
                                    goodClusts.append(clusPnts)
                                    prevHit = True
                                elif prevHit == True and len(goodClusts) == 0: # If one of the clusters from the tighter bin passed the percentage threshold, take the cluster with the most points within 2 Ang of a heavy ligand atom
                                    if (closePntCnt/len(clusPnts)) > failed_bestPercentage:
                                        failed_bestInd = ind
                                        failed_bestPercentage = (closePntCnt/len(clusPnts))
                    if failed_bestInd != '' and len(goodClusts) == 0: # No clusters pass the threshold but clusters from smaller bins did pass, & any clusters have any points within 2 Ang of a heavy ligand atom
                        goodClusts.append(pntArr[clusters.labels_ == failed_bestInd])
                                    
                    
                    allPocketPnts = [p for c in goodClusts for p in c] # Flatten separate pockets
                                        
                    logFH.write(str(binThresholds[i]) + ' Ang binning threshold:\n')
                    logFH.write('\t'+str(round(len(allPocketPnts) * gridPntSpacing**3,2)) + ' Ang^3 of binding pocket volume from ')
                    if len(goodClusts) == 1: logFH.write('1 binding pocket\n')
                    else: logFH.write(str(len(goodClusts)) + ' binding pockets\n')
                    
                    # Write bin points to a mol2 file
                    coords2Mol2(mol2Dir + pdb + '_' + bs.replace(':','-') + '_' + str(binThresholds[i]) + 'AngBin' + '.mol2', allPocketPnts) # Save out all remaining binding pockets points to a mol2 file, named by the PDB ID, the parent ligand PLIP ID with ":" swapped for "-", and the binning threshold  used

            
            
            end = time.time()
            # print(end-start)
            logFH.write('\n\n' + str(round(end-start, 2)) + ' seconds to process ' + str(len(allPLIP[pdb].bsites)) + ' binding site(s) from pdb ' + pdb + '\n\n\n')
            logFH.flush()
    
        
    globalEnd = time.time()
    print(round((globalEnd-globalStart)/60, 2), 'minutes to process',bsCnt, 'binding sites from',plipCnt,'PDB files')
        
        
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
