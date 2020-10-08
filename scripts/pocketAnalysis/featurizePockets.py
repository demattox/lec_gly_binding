#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  5 14:29:43 2020

@author: dmattox
"""

import os
import glob
import dill
import time
import collections
# import itertools
import numpy as np
# import scipy.spatial
# import scipy.stats

import matplotlib.pyplot as plt
import seaborn as sns

oWD = os.getcwd()
oWD += '/'

import lec_gly as LecGly
from bSiteResiFeatures import plipFile
from pocketAnalysis.generateVoxels import mol2Dir, GRID_POINT_SPACING, POCKET_THRESHOLDS
from pocketAnalysis.getPocketSurfacePoints import pcktSurfDir

os.chdir(LecGly.homeDir)


############################################

# np.random.seed(27)

# lineSearchResolution = 0.25 # How far apart to space the points on the line when looking for neighboring pocket points

# def orthProj(norm, planePnt, mvingPnt):
#     ''' Orthogonally projects mvingPnt onto the plane containing the planePnt with the normal vector norm '''
#     t = (np.dot(norm, planePnt) - np.dot(norm, mvingPnt)) / np.linalg.norm(norm)        
#     return mvingPnt + t*norm

NUM_SCALED_BINS = 40

############################################

if __name__ == '__main__':
    
    ############################################
      
    distPlotsDir = './analysis/bSitePlots/allPairwiseDists/'
    if not os.path.exists(distPlotsDir):
        os.makedirs(distPlotsDir)
    
    distributionFeatsOut = oWD + 'd2DistributionFeatures.csv'
    allDistsOut = oWD + 'allD2Distances.csv'
    scaledDistsOut = oWD + 'scaledBins_D2Distances.csv'
        
    ############################################
        
    with open(plipFile, 'rb') as pickleFH:
        allPLIP = dill.load(pickleFH)
    
    # with open("./data/unilectin/icodeResi2PLIPmapping.p", "rb")as pickleFH:
    #     zeroResis = dill.load(pickleFH)
   
    
    ############################################
    
    sns.set(rc={'axes.facecolor':'lightgrey', 'figure.facecolor':'lightgrey'})
    
    threshColors = ['magenta', 'red', 'orange', 'yellow']
    colDict = {}
    if len(threshColors) == len(POCKET_THRESHOLDS):
        colDict = {t : c for t,c in zip(POCKET_THRESHOLDS, threshColors)}
    else:
        print("warning, more pocket thresholds than expected")
    
    allDistsDict = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(list))) # three level dict to store lists of all pairwise distances [PDB][BS][Threshold] : [distances]
    binnedAllDists = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(lambda: collections.Counter()))) # Store binned D2 measurments
    scaled_binnedAllDists = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(dict))) # Scale bins for D2 measurments s.t. there are always 40 bins evenly divided b/w the min & max observed distances
    
    pocketMeasures = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(dict)))
    
    # clDistsDict = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(list)))
    # binnedCLdists = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(lambda: collections.Counter())))
    
    # runList = list(allPLIP.keys())
    
    runList = [] # List to hold subset of PDB ids to run in batch processing
    files = os.listdir(oWD)
    batchList = [f for f in files if f[:7] == 'pdbList'][0]
    with open(oWD + batchList, 'r') as inFH:
        for line in inFH:
            runList.append(line.strip())

    start = time.time()
    
    header = ['vol','pcntSurf','var', 'med', 'q1', 'q3', 'skew', 'majorMaxima', 'minorMaxima']
    with open(distributionFeatsOut, 'w') as outFH:
        outFH.write('bsite')
        for k in colDict.keys():
            outFH.write(',' + ','.join([t + '_'+str(k)+'Ang' for t in header]))
        outFH.write('\n')
         
        for pCnt,pdb in enumerate(runList):
            # if len(runList) >=10:
            #     if pCnt % (round(len(runList)/10, 0)) == 0: print("Processing PBD file number " + str(pCnt) + ' of ' + str(len(runList)), flush = True)
            
            print(pdb)
            
            # pdb = '2CL8'
            # bs = 'BGC:A:1247'
            # pdb = '3LL2'
            # bs='MAN:B:1'
            # pdb = '1LES'
            # bs='GLC:A:205'
        
               
            for bs in allPLIP[pdb].bsites:
                 
                # pdbIN = [f for f in pdbFiles if f[-8:-4] == pdb][0]
                # parser = Bio.PDB.PDBParser(QUIET=True)
                # model = parser.get_structure(pdb, pdbIN)[0] # grab the first model in PDB file
        
                # ligAtms = getLigAtoms(model, allPLIP[pdb].bsites[bs])
                # CL = getCentroid(ligAtms) # Ligand centroid
        
                
                fullPcktMol2Files = glob.glob(mol2Dir + pdb + '_' + bs.replace(':','-') + '_*')
                surfFiles = glob.glob(pcktSurfDir + pdb + '_' + bs.replace(':','-') + '_*')
                
                if not fullPcktMol2Files:
                    print('WARNING: Skipping',pdb,bs + ', no mol2 files for pocket found')
                    outFH.write(pdb + '_' + bs)
                    for b in colDict.keys():
                        outFH.write(',' + ','.join([ str(0) for t in header ]) )
                    outFH.write('\n')

                else:
                    for f in fullPcktMol2Files:
                        pcktPnts = LecGly.getMol2Pnts(f)
                        binThresh = f.split("_")[-1]
                        binThresh = int(binThresh.split('AngBin')[0])
                        
                        if not pcktPnts:
                            print('WARNING: No volume detected in binding pocket, skipping', pdb, bs, binThresh, 'Ang pocket')
                            for meas in header:
                                pocketMeasures[pdb][bs][binThresh][meas] = 0
                        else:
                            
                            # surfF = surfFiles[0]
                            surfF = [fh for fh in surfFiles if pdb + '_' + bs.replace(':','-') + '_' + str(binThresh) in fh][0]
                            surfacePnts = LecGly.getMol2Pnts(surfF)
                            
                            # pcktTree = scipy.spatial.KDTree(pcktPnts) # Implement a kd-tree containing all of the points in the pocket
                            
                            pocketMeasures[pdb][bs][binThresh]['vol'] = len(pcktPnts) * round(GRID_POINT_SPACING**3,5) # Volume of the pocket in Ang^3
                            pocketMeasures[pdb][bs][binThresh]['pcntSurf'] = round(len(surfacePnts)/len(pcktPnts), 5)
                                                    
                            # for i,p1 in enumerate(surfacePnts):
                            #     d = eucDist(p1,CL)
                            #     clDistsDict[pdb][bs][binThresh].append(d)
                            #     binD = float(d // 1)
                            #     if d%1 >= 0.5:
                            #         binD += 0.5
                            #     binnedCLdists[pdb][bs][binThresh][binD]+= 1
                                                        
                            for i,p1 in enumerate(surfacePnts):
                                for j in range(i+1, len(surfacePnts)):
                                    d = LecGly.eucDist(p1,surfacePnts[j])
                                    allDistsDict[pdb][bs][binThresh].append(d)
                                    binD = float(d // 1)
                                    if d%1 >= 0.5:
                                        binD += 0.5
                                    binnedAllDists[pdb][bs][binThresh][binD]+= 1
                            
                            binVals = np.linspace(min(allDistsDict[pdb][bs][binThresh]),  max(allDistsDict[pdb][bs][binThresh])+0.01, NUM_SCALED_BINS + 1)  # Create scaled bin values to have the same number of bins for each shape
                            allDists = np.array(allDistsDict[pdb][bs][binThresh])  # Create a new array object form the list of all D2 measurements for faster binning
                            for i in range(len(binVals)-1):
                                scaled_binnedAllDists[pdb][bs][binThresh][i] = sum((binVals[i] <= allDists) & (allDists < binVals[i+1]))
                            # print(sum(scaled_binnedAllDists[pdb][bs][binThresh].values()))
                            
                            pocketMeasures[pdb][bs][binThresh]['var'] = np.var(allDistsDict[pdb][bs][binThresh])
                            pocketMeasures[pdb][bs][binThresh]['med'] = np.median(allDistsDict[pdb][bs][binThresh])
                            pocketMeasures[pdb][bs][binThresh]['q1'] = np.percentile(allDistsDict[pdb][bs][binThresh], 25, interpolation='midpoint')
                            pocketMeasures[pdb][bs][binThresh]['q3'] = np.percentile(allDistsDict[pdb][bs][binThresh], 75, interpolation='midpoint')
                            pocketMeasures[pdb][bs][binThresh]['skew'] = np.mean(allDistsDict[pdb][bs][binThresh]) - np.median(allDistsDict[pdb][bs][binThresh]) 
                    
                    revLst = sorted(list(pocketMeasures[pdb][bs].keys()))[::-1]
                    for binThresh in revLst:
                        x_d = list(binnedAllDists[pdb][bs][binThresh].keys())
                        if x_d:
                            x_d.sort()
                            
                            y_d = [binnedAllDists[pdb][bs][binThresh][d] for d in x_d]
                            if len(y_d) > 9:
                                y_d_s9 = LecGly.smooth(np.array(y_d), window_len = 9)
                                y_d_s5 = LecGly.smooth(np.array(y_d), window_len = 5)
                                pocketMeasures[pdb][bs][binThresh]['majorMaxima'] = LecGly.cntLocalMaxima(y_d_s9)
                                pocketMeasures[pdb][bs][binThresh]['minorMaxima'] = LecGly.cntLocalMaxima(y_d_s5)
                            elif len(y_d) > 3:
                                print(pdb,bs,binThresh, 'Ang bin only has', len(y_d), 'distance bins')
                                pocketMeasures[pdb][bs][binThresh]['majorMaxima'] = LecGly.cntLocalMaxima(y_d)
                                pocketMeasures[pdb][bs][binThresh]['minorMaxima'] = LecGly.cntLocalMaxima(y_d)
                            else:
                                print(pdb,bs,binThresh, 'Ang bin only has', len(y_d), 'distance bins')
                                pocketMeasures[pdb][bs][binThresh]['majorMaxima'] = 1
                                pocketMeasures[pdb][bs][binThresh]['minorMaxima'] = 1
                            
                            #########
                            # Plot y_d as a histogram, then s9 and s5 over the bars
                            plt.bar(x_d, y_d, alpha=0.3, color = colDict[binThresh], width = 0.5 )
                            
                            plt.plot(x_d, y_d_s9, alpha=0.8, color = colDict[binThresh])
                            plt.plot(x_d, y_d_s5, alpha=0.8, color = colDict[binThresh])
                            # plt.scatter(x_d, y_d, alpha=0.8, color = colDict[binThresh])
                        else:
                            pocketMeasures[pdb][bs][binThresh]['majorMaxima'] = 0
                            pocketMeasures[pdb][bs][binThresh]['minorMaxima'] = 0
                        
                    plt.title('All pairwise distances b/w surface points for ' + pdb + ' ' + bs)
                    plt.xlabel('Distance (Angstroms)')
                    plt.ylabel('Count')
                    
                    plt.savefig(distPlotsDir + pdb + '_' + bs.replace(':','-') +'.jpg', dpi = 240)
                    plt.close()
                    
                    outFH.write(pdb + '_' + bs)
                    for binThresh in colDict.keys():
                        outFH.write(',' + ','.join([ str(pocketMeasures[pdb][bs][binThresh][t]) for t in header ]) )
                    outFH.write('\n')
                    
            
    end = time.time()
    print(round(end-start,2),'seconds for', len(runList), 'PDB files')
    
    with open(scaledDistsOut, 'w') as outFH:
        outFH.write('bsite')
        for k in colDict.keys():
            outFH.write(',' +','.join([str(d)+'_'+str(k)+'Ang' for d in range(NUM_SCALED_BINS)]))
        outFH.write('\n')
        for pdb in scaled_binnedAllDists.keys():
            for bs in scaled_binnedAllDists[pdb].keys():
                outFH.write(pdb + '_' + bs)
                for binThresh in colDict.keys():
                    if binThresh in scaled_binnedAllDists[pdb][bs].keys():
                        outFH.write(',' + ','.join([ str(scaled_binnedAllDists[pdb][bs][binThresh][d]) for d in range(NUM_SCALED_BINS) ]) )
                    else:
                        outFH.write(',' + ','.join([ str(0) for d in range(NUM_SCALED_BINS) ]) )
                outFH.write('\n')
    
    # with open("binnedAllDists.p", "wb")as pickleFH:
    #     dill.dump(binnedAllDists, pickleFH)
    
    topAll = 0
    for pdb in binnedAllDists.keys():
        for bs in binnedAllDists[pdb].keys():
            for binThresh in binnedAllDists[pdb][bs].keys():
                if binnedAllDists[pdb][bs][binThresh]:
                    if topAll < max(binnedAllDists[pdb][bs][binThresh].keys()):
                        topAll = max(binnedAllDists[pdb][bs][binThresh].keys())

    maxDist = topAll
    binDists = np.arange(0.5,maxDist+0.5, 0.5)
    with open(allDistsOut, 'w') as outFH:
        outFH.write('bsite')
        for k in colDict.keys():
            outFH.write(',' +','.join([str(d)+'_'+str(k)+'Ang' for d in binDists]))
        outFH.write('\n')
        
        for pdb in binnedAllDists.keys():
            for bs in binnedAllDists[pdb].keys():
                outFH.write(pdb + '_' + bs)
                for binThresh in colDict.keys():
                    if binThresh in binnedAllDists[pdb][bs].keys():
                        outFH.write(',' + ','.join([ str(binnedAllDists[pdb][bs][binThresh][d]) for d in binDists ]) )
                    else:
                        outFH.write(',' + ','.join([ str(0) for d in binDists ]) )
                outFH.write('\n')
     

        
        
        
        
        
