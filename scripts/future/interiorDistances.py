#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 21 16:42:36 2020

@author: dmattox
"""

import os, glob, time, collections
import scipy.spatial
import numpy as np

import seaborn as sns
import matplotlib.pyplot as plt

from checkPLIPligands import eucDist
from getPocketSurfacePoints import getMol2Pnts

pocketDir = '/Users/dmattox/cbk/glycan_binding/data/unilectin/structures/bSites/bsitePockets/'
pocketSurfaceDir = '/Users/dmattox/cbk/glycan_binding/data/unilectin/structures/bSites/pocketSurfacePoints/'
# plipDat = '/Users/dmattox/cbk/glycan_binding/data/unilectin/cleanPLIPrecords.p'

pdb = '3LL2'
bs = 'MAN:A:122'

fullPcktMol2Files = glob.glob(pocketDir + pdb + '_' + bs.replace(':','-') + '_*')
surfFiles = glob.glob(pocketSurfaceDir + pdb + '_' + bs.replace(':','-') + '_*')

f = '/Users/dmattox/cbk/glycan_binding/data/unilectin/structures/bSites/bsitePockets/3LL2_MAN-A-122_4AngBin.mol2'

pcktPnts = getMol2Pnts(f)

pcktTree = scipy.spatial.KDTree(pcktPnts) # Implement a kd-tree containing all of the points in the pocket

f = '/Users/dmattox/cbk/glycan_binding/data/unilectin/structures/bSites/pocketSurfacePoints/3LL2_MAN-A-122_4_surfacepnts.mol2'

surfPnts = getMol2Pnts(f)


gridSpace = (0.5)**(1/3)
step = gridSpace
search = gridSpace

start = time.time()

dists = []
for i,p1 in enumerate(surfPnts):
    for j in range(i+1, len(surfPnts)):
        
        p2 = surfPnts[j]
        v = np.array(p2)-np.array(p1)
        v = v/np.linalg.norm(v) # unit vector between pnts
        
        interior = True
        p = p1 + step*v
        while eucDist(p,p2) > step and interior:
            # print('pseudoatom p, pos=[' + ','.join([str(n) for n in p]) + ']')
            if len(pcktTree.query_ball_point(x = p, r = search)) == 0:
                interior = False
            p += step*v
        
        if interior == True:
            dists.append(eucDist(p1,p2))
        
end = time.time()
print(round(end-start, 2), 'seconds')


binned = collections.Counter()
for d in dists:
    binD = float(d // 1)
    if d%1 >= 0.5:
        binD += 0.5
    binned[binD]+= 1

x_d = list(binned.keys())
if x_d:
    x_d.sort()
    y_d = [binned[d] for d in x_d]
    
    if len(y_d) > 9:
        y_d_s9 = smooth(np.array(y_d), window_len = 9)
        y_d_s5 = smooth(np.array(y_d), window_len = 5)

plt.bar(x_d, y_d, alpha=0.3, color = 'magenta', width = 0.5 )
plt.plot(x_d, y_d_s9, alpha=0.8, color = 'magenta')
plt.plot(x_d, y_d_s5, alpha=0.8, color = 'magenta')
plt.title('INTERIOR pairwise distances b/w surface points for 3LL2')
plt.xlabel('Distance (Angstroms)')
plt.ylabel('Count')

oldDists = collections.Counter()
for i,p1 in enumerate(surfPnts):
    for j in range(i+1, len(surfPnts)):
        d = eucDist(p1,surfPnts[j])
        binD = float(d // 1)
        if d%1 >= 0.5:
            binD += 0.5
        oldDists[binD]+= 1

x_d = list(oldDists.keys())
if x_d:
    x_d.sort()
    y_d = [oldDists[d] for d in x_d]
    
    if len(y_d) > 9:
        y_d_s9 = smooth(np.array(y_d), window_len = 9)
        y_d_s5 = smooth(np.array(y_d), window_len = 5)

plt.bar(x_d, y_d, alpha=0.3, color = 'magenta', width = 0.5 )
plt.plot(x_d, y_d_s9, alpha=0.8, color = 'magenta')
plt.plot(x_d, y_d_s5, alpha=0.8, color = 'magenta')
plt.title('ALL pairwise distances b/w surface points for 3LL2')
plt.xlabel('Distance (Angstroms)')
plt.ylabel('Count')



print('pseudoatom p1, pos=[' + ','.join([str(n) for n in p1]) + ']')
print('pseudoatom p2, pos=[' + ','.join([str(n) for n in p2]) + ']')

print('pseudoatom p, pos=[' + ','.join([str(n) for n in p]) + ']')
