#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 17 10:44:57 2020

@author: dmattox
"""

import os
import glob
import dill
import time
# import itertools, collections
import numpy as np
import Bio.PDB
# import pliptool.plip.modules.plipxml as plipxml
from scipy.spatial import ConvexHull
from sklearn.cluster import DBSCAN

oWD = os.getcwd()
oWD += '/'

import lec_gly as LecGly
import checkPLIPligands as LigCheck
from bSiteResiFeatures import plipFile
from pocketAnalysis.pdb2wrl_surface import wrlDir

os.chdir(LecGly.homeDir)

############################################

mol2Dir = './data/structures/bsites/bsiteVoxels/'
if not os.path.exists(mol2Dir):
    os.makedirs(mol2Dir)

logFile = oWD + 'generateVoxels_logfile.txt' # Generates a log file in the directory from which the script is called (thoght to be the subdirectory containing a batch processing list of PDBs to run)

############################################

MIN_CUBE_EDGE = 20 # define the minimum length of the cube edge
GRID_POINT_SPACING = (0.5)**(1/3) # how far apart to place the grid points
PRO_ATM_THRESH = 2.5 # Threshold distance to a heavy protein atom to consider a point within the protein surface
MAX_DIST_FROM_LIG = 10 # Threshold for maximum allowed distance for any grid point from a heavy ligand atom
POCKET_THRESHOLDS = [4,6,8] # Distance values for binning points prior to clustering
POCKET_THRESHOLDS.append(MAX_DIST_FROM_LIG)


DBSCAN_EPS = np.sqrt(3/4) * GRID_POINT_SPACING * 2 + 0.01 # The maximum distance between two samples for one to be considered as in the neighborhood of the other. Setting to the max length between directly neighboring points, allowing for 26 possible nieghbors 
DBSCAN_MIN_SAMPS = 13 # The number of samples (or total weight) in a neighborhood for a point to be considered as a core point.
CLUSTER_SIZE_THRESH = 15 # Prune clusters of pocket points with less than 15 Ang**3 of volume
CLUSTER_SIZE_THRESH = CLUSTER_SIZE_THRESH/(GRID_POINT_SPACING**3) # convert from Ang**2 to number of points based on spacing
PNTS_NEAR_LIG = 0.1 # Prune clusters with fewer than 10% of its points within 2 Ang of any heavy ligand atom

############################################

if __name__ == '__main__':
    
    pdbFiles = glob.glob(LigCheck.originalPDBsDir + '*.pdb')
    
    with open(plipFile, 'rb') as pickleFH:
        allPLIP = dill.load(pickleFH)
    
    ############################################
    
    globalStart = time.time()
        
    plipCnt = 0
    bsCnt = 0
    
    runList = [] # List to hold subset of PDB ids to run in batch processing
    files = os.listdir(oWD)
    batchList = [f for f in files if f[:7] == 'pdbList'][0]
    with open(oWD + batchList, 'r') as inFH:
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
                ligAtms = LecGly.getLigAtoms(model, allPLIP[pdb].bsites[bs])
                
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
                
                CL = LecGly.getCentroid(ligAtms) # Ligand centroid
                # CP = LecGly.getCentroid(proAtms) # Binding sites residues centroid
                # cDist = LecGly.eucDist(CL,CP) # Distance between ligand centroid and binding site residues centroid
                
                ligAtmDists = [LecGly.eucDist(CL, a.coord) for a in ligAtms] # Distance between ligand centroid and all atoms in the ligand
                cubeEdge = MIN_CUBE_EDGE if max(ligAtmDists) <= (MIN_CUBE_EDGE/2.5) else round(2.5*max(ligAtmDists)) # Set the size of the cube to be at least 20 Ang along each edge, or larger for larger ligands. Round to keep grid points even
                # apprxCL = [round(p) for p in CL] # rounded coordinates of ligand centroid
                
                cubeLimits = [[d + cubeEdge/2, d - cubeEdge/2] for d in CL]
                # cubeVerts = list(itertools.product(*cubeLimits))
                # coords2Mol2('data/unilectin/structures/testcubeverts.mol2', cubeVerts)
                
                gridPnts = [] # hold the grid points within the cube defined above centered on the ligand centroid
                for x in np.arange(min(cubeLimits[0]),max(cubeLimits[0])+GRID_POINT_SPACING, step = GRID_POINT_SPACING):
                    for y in np.arange(min(cubeLimits[1]),max(cubeLimits[1])+GRID_POINT_SPACING, step = GRID_POINT_SPACING):
                        for z in np.arange(min(cubeLimits[2]),max(cubeLimits[2])+GRID_POINT_SPACING, step = GRID_POINT_SPACING):
                            gridPnts.append([x,y,z])
                
                # coords2Mol2('data/unilectin/structures/testgridInit.mol2', gridPnts)                
                
                wrlFile = wrlDir + pdb + '_' + bs.replace(':','-') + '.wrl'
                proCoords = LecGly.getPntsFromWRL(wrlFile) # Points representing the van der Waals surface of the binding site residues
                
                # LecGly.coords2Mol2('data/unilectin/structures/' + bs+'surface.mol2', proCoords)
                
                ########################################################
                # Begin excluding points          
                      
                # Exclude points within protein surface
                nonBuriedPnts = []
                for p in gridPnts:
                    closeProAtms = searcher.search(p,radius = PRO_ATM_THRESH, level = 'A')
                    if len(closeProAtms) == 0:
                        nonBuriedPnts.append(p)
                # coords2Mol2('data/unilectin/structures/testgridCloseInhullNotburied.mol2', nonBuriedPnts)
                logFH.write('From ' + str(len(gridPnts)) + ' points, ' + str(len(nonBuriedPnts)) + ' points found outside of the protein surface\n')
                # print('From ' + str(len(gridPnts)) + ' points, ' + str(len(nonBuriedPnts)) + ' points found outside of the protein surface\n')
                

                # Exclude points far from the ligand
                binPnts = [[] for b in POCKET_THRESHOLDS] # Holds points in each bin for each binning distance threshold lower than the max threshold
                for p in nonBuriedPnts:
                    dists = [LecGly.eucDist(p,a.coord) for a in ligAtms]
                    if min(dists) <= MAX_DIST_FROM_LIG:
                        d = min(dists)
                        for i,b in enumerate(POCKET_THRESHOLDS):
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
                for p in binPnts[-1]: # Exclude points from largest cpocket threshold that are exterior to convex hull first, then check for those points prescence in the smaller pocket threshold point sets
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
                
                # prevHit = False
                prevNum_Vol = [0,0] # If a pocket from a lower threshold had 1+ clusters pass the size and ligand proximity threshold, store the number of clusters prevNum_Vol[0] and the total volume of the pocket prevNum_Vol[1]
                for i,pntSet in enumerate(binPnts):
                    pntArr = np.array(pntSet)
                    clusters = DBSCAN(eps = DBSCAN_EPS, min_samples= DBSCAN_MIN_SAMPS).fit(pntArr) #DBSCAN to cluster non-excluded points together into pockets
                    # print(i, set(clusters.labels_))
                    
                    # for ind in set(clusters.labels_):
                    #     print('Cluster ind',str(ind))
                    #     print('\t',str(sum(clusters.labels_ == ind)), 'points')
                    #     coords2Mol2('data/unilectin/structures/testClust'+ '_' + str(i) + '_' +str(ind) +'.mol2', pntArr[clusters.labels_ == ind])
                    
                    ### First pass method, take acceptable clusters as you go and if you get 0 clusters, take the next best one if you got any clusters at a smaller pocket threshold
                    # goodClusts = []
                    # failed_bestInd = ''
                    # failed_bestPercentage = 0 # If a pocket has clusters that pass the threshold initially but later fails to hit the 10% (PNTS_NEAR_LIG) threshold as the pocket expands, for each threshold that fails keep the cluster with the most points near the ligand atoms
                    # for ind in set(clusters.labels_):
                    #     if ind != -1: # Don't consider the noise cluster
                    #         clusPnts =  pntArr[clusters.labels_ == ind]
                    #         if len(clusPnts) > CLUSTER_SIZE_THRESH: # Filter out clusters with too few points
                    #             closePntCnt = 0
                    #             for p in clusPnts:
                    #                 dists = [LecGly.eucDist(p, a.coord) for a in ligAtms]
                    #                 if min(dists) <= 2:
                    #                     closePntCnt += 1
                    #             if closePntCnt >= len(clusPnts)*PNTS_NEAR_LIG:
                    #                 goodClusts.append(clusPnts)
                    #                 prevHit = True
                    #             elif prevHit == True and len(goodClusts) == 0: # If one of the clusters from the tighter bin passed the percentage threshold, take the cluster with the most points within 2 Ang of a heavy ligand atom
                    #                 if (closePntCnt/len(clusPnts)) > failed_bestPercentage:
                    #                     failed_bestInd = ind
                    #                     failed_bestPercentage = (closePntCnt/len(clusPnts))
                    # if failed_bestInd != '' and len(goodClusts) == 0: # No clusters pass the threshold but clusters from smaller bins did pass, & any clusters have any points within 2 Ang of a heavy ligand atom
                    #     goodClusts.append(pntArr[clusters.labels_ == failed_bestInd])
                    
                    largeClusInds = []
                    largeClusPercents = []
                    goodClusts = []
                    for ind in set(clusters.labels_):
                        if ind != -1: # Don't consider the noise cluster
                            clusPnts =  pntArr[clusters.labels_ == ind]
                            if len(clusPnts) > CLUSTER_SIZE_THRESH: # Filter out clusters with too few points
                                largeClusInds.append(ind)
                                closePntCnt = 0
                                for p in clusPnts:
                                    dists = [LecGly.eucDist(p, a.coord) for a in ligAtms]
                                    if min(dists) <= 2:
                                        closePntCnt += 1
                                largeClusPercents.append(closePntCnt / len(clusPnts))
                                
                    # Take the top clusters of acceptable size that pass the 10% PNTS_NEAR_LIG cutoff. Recovers clusters that drop below the 10% cutoff as the pocket threshold grows
                    largeClusts = sorted(zip(largeClusInds,largeClusPercents), reverse = True, key = lambda pair: pair[1])
                    if prevNum_Vol[0] == 0:
                        goodClusts = [pntArr[clusters.labels_ == ind] for ind,percent in largeClusts if percent >= PNTS_NEAR_LIG]
                        prevNum_Vol = [len(goodClusts), sum(map(len, goodClusts))] # count the number of clusters and the number of points in all the clusters combined
                    else:
                        goodClusts = [pntArr[clusters.labels_ == ind] for ind,percent in largeClusts if percent >= PNTS_NEAR_LIG]
                        newNum_Vol = [len(goodClusts), sum(map(len, goodClusts))]
                        if newNum_Vol[0] < prevNum_Vol[0] and newNum_Vol[1] < prevNum_Vol[1]: # Catch cases where a cluster from a smaller threshold drops below the 10% cutoff while growing, take the next closest cluster in this case
                            missedClustCnt = prevNum_Vol[0] - newNum_Vol[0]
                            for j in range(newNum_Vol[0], newNum_Vol[0] + missedClustCnt):
                                goodClusts.append(pntArr[clusters.labels_ == largeClusts[j][0]])
                            prevNum_Vol = [len(goodClusts), sum(map(len, goodClusts))]
                        else:
                            prevNum_Vol = newNum_Vol[:]
                    allPocketPnts = [p for c in goodClusts for p in c] # Flatten separate pockets
                                        
                    logFH.write(str(POCKET_THRESHOLDS[i]) + ' Ang binning threshold:\n')
                    logFH.write('\t'+str(round(len(allPocketPnts) * GRID_POINT_SPACING**3,2)) + ' Ang^3 of binding pocket volume from ')
                    if len(goodClusts) == 1: logFH.write('1 binding pocket\n')
                    else: logFH.write(str(len(goodClusts)) + ' binding pockets\n')
                    
                    # Write thresholded points to a mol2 file
                    LecGly.coords2Mol2(mol2Dir + pdb + '_' + bs.replace(':','-') + '_' + str(POCKET_THRESHOLDS[i]) + 'AngBin' + '.mol2', allPocketPnts) # Save out all remaining binding pockets points to a mol2 file, named by the PDB ID, the parent ligand PLIP ID with ":" swapped for "-", and the binning threshold  used

            
            
            end = time.time()
            # print(end-start)
            logFH.write('\n\n' + str(round(end-start, 2)) + ' seconds to process ' + str(len(allPLIP[pdb].bsites)) + ' binding site(s) from pdb ' + pdb + '\n\n\n')
            logFH.flush()
    
        
    globalEnd = time.time()
    print(round((globalEnd-globalStart)/60, 2), 'minutes to process',bsCnt, 'binding sites from',plipCnt,'PDB files')
        
        
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
