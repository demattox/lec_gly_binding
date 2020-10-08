#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 28 07:16:01 2020

@author: dmattox
"""

import os
import glob
import collections
from string import ascii_lowercase, ascii_uppercase

####################


####
os.chdir('/Users/dmattox/cbk/lec_gly_binding')
####


cdhitChainDir = './data/structures/holo/seqs/allSeqsCDHIT90/'

outFile = './data/structures/holo/seqs/nrSeqs/allSeqs.fst'

#####################

clust_lsts = sorted(glob.glob(cdhitChainDir + '*_nr.fst.clstr'))
clust_seqs = sorted(glob.glob(cdhitChainDir + '*_nr.fst'))

[c[:-6] for c in clust_lsts] == clust_seqs  # File lists match

chainCompare = {letter: i for i, letter in enumerate(ascii_uppercase+ascii_lowercase)}  # Assign numerical values to the first 52 chain indices to compare order

#####################
with open(outFile, 'w') as outFH:
    for clustF, seqF in zip(clust_lsts, clust_seqs):
        out = ''  # hold the final sequence to write to the file
        
        seqs = {}  # stores the representative sequence for each cluster, indexed by chain ID
        with open(seqF, 'r') as seqIn:
            readSeq = False
            for line in seqIn:
                if not readSeq and line[0] == '>':
                    chainID = line[1]
                    readSeq = True
                else:
                    s = line.strip()
                    seqs[chainID] = s
                    readSeq = False
        
        if len(seqs) == 1:  # if only one cluster, save it out
            out = [v for v in seqs.values()][0]
        else:
            clustMembers = collections.defaultdict(list)  # Store the chainIDs for all the chains in a cluster
            clustReps = {}  # Store the chainIDs of the representative sequence for each cluster
            with open(clustF, 'r') as clustIn:
                for line in clustIn:
                    if line[0] == '>':
                        clustID = line[1:].strip().split(' ')[-1]  # Get the cluster index at the end of a string like: '>Cluster 3'
                    else:
                        member = line.strip().split()
                        mID = member[2]  # Get the name of the cluster member
                        clustMembers[clustID].append(mID[1])  # Get the chainID from the name, formatted as ">A..."
                        if member[3] == '*':  # Store the cluster representative
                            clustReps[clustID] = mID[1]
            for clustID in clustMembers.keys():  # Order the clusters by the chain IDs closest to A
                if len(clustMembers[clustID]) == 1:
                    clustMembers[clustID] = clustMembers[clustID][0]
                else:
                    memberLst = [chainCompare[m] if m in chainCompare.keys() else m for m in clustMembers[clustID] if m in chainCompare.keys()]
                    clustMembers[clustID] = clustMembers[clustID][min(range(len(memberLst)), key=memberLst.__getitem__)]  # get the chainID closest to 'A' for each cluster
            clustOrder = [itm[0] for itm in sorted(clustMembers.items(), key=lambda x: x[1])]  # order the dictionary by the values and extract the keys from the list of tuples
            
            for clustID in clustOrder:
                out += seqs[clustReps[clustID]]
        
        pdbID = seqF.split('/')[-1][:4]
        outFH.write('>'+pdbID+'\n')
        outFH.write(out+'\n')
        
        
        
        
        
        
        
        
        
        
        
        