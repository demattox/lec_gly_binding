#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  7 13:32:13 2020

@author: dmattox

"""


import os
import glob
# import math
import itertools
import collections
import dill
import copy
import numpy as np
import matplotlib.pyplot as plt
#import sys
from rdkit import Chem
from rdkit.Chem import Crippen
from rdkit import RDLogger
import pliptool.plip.modules.plipxml as plipxml
import Bio.PDB
from Bio.PDB.DSSP import make_dssp_dict

import lec_gly as LecGly

os.chdir(LecGly.homeDir)

########################################################################
# Paths and file handles
# manualFilter  = False 


plipDir = './data/plip/PLIPreports/'
PLIPpdbDir = './data/structures/holo/plipFixed/'
originalPDBsDir = './data/structures/holo/origPDB/'
unilecDataFH = './data/unilectin/boundUnilectinData.tsv'
dsspDir = './data/dssp/dsspOut/'
qcDir = './analysis/QC/'
# outTSV = './analysis/bsiteResiFeatures.tsv'

# if manualFilter:
#     manualDropList = './data/unilectin/manuallyFilteredBSites.csv'

########################################################################
# Thresholds and defaults
COV_THRESH = 1.7 # 1.7 angstrom covalent bond threshold, drop to ~1.43 for C-O bond length?
CENTROID_THRESH = 1.7 # If the centroids of the residues potentially missing a covalent link are less than 1.7 Ang apart, don't inlcude
CA_THRESH = 3.0 # "Max. distance between metal ion and interacting atom (Harding, 2001)" from plip/plip/modules/config.py
CHAIN_LEN_THRESH =  15 # If a glycan is covalently linked to a protein chain less than this threshold, glycan treated as ligand, otherwise treated as glycosylation and removed from analysis
potentialGlycosylationResis = ['SER', 'THR', 'ASN', 'TRP', 'CYS'] # List of residues that glycans might be covalently attached to in a glycosylation event

cifLst = ['1BOS'] # List of PDB IDs that were originally downloaded in mmCIF format and convert PDB format does not load into biopython properly



########################################################################
# Load data

uniprotMap = collections.defaultdict(list) # Dictionary to store PDB IDs associated with each UniProt ID
unilecData = {} # dictionary for importing metadata for each lectin provided by UniLectin
unilecCols = [] # for storing column headers of interest from

with open(unilecDataFH, 'r') as dataIn:
    for line in dataIn:
        line = line.split('\t')
        if (line[0] == 'pdb'):
                unilecCols = [line[1]] + line[2:5] + line[7:12] # Get column names of interest
        else:
            unilecData[line[0]] = [line[1]] + line[2:5] + line[7:12] # Get data from columns of interest
            uniprotMap[line[1]].append(line[0]) # Get uniprot-PDB associations

plipFiles = glob.glob(plipDir + '*.xml')
pdbFiles = glob.glob(PLIPpdbDir + '*.pdb')
cifFiles = glob.glob(originalPDBsDir + '*.cif')
dsspFiles = glob.glob(dsspDir + '*.dssp') # 

plipFiles = [f for f in plipFiles if f[-8:-4] in unilecData.keys()] # Remove XML reports not in set of Unilectin data for analysis
pdbFiles = [f for f in pdbFiles if f[-8:-4] in unilecData.keys()] # Same but with PDB files

if (len(plipFiles) != len(pdbFiles)):
    print("Warning, different numbers of PDB files and PLIP reports for consideration")
elif (len(plipFiles) != len(unilecData)):
    print("Warning, numbers of PLIP reports and PDB files do not match number of structures in dataframe")

pdbIDs = list(unilecData.keys())

smilesDict = {} # Holds mappings from 3 letter PDB heteroatom IDs to SMILES stings
with open('./data/supp/Components-smiles-stereo-oe.tsv', 'r') as compIN:
    for line in compIN:
        line = line.split('\t')
        if len(line) == 3:
            smilesDict[line[1]] = line[0]
                
########################################################################

if __name__ == '__main__':
    RDLogger.DisableLog('rdApp.*') # Suppress rdkit valencey warnings (ligand structures not always clean/correct)

    
    # if manualFilter:
    #     dropList = collections.defaultdict(lambda: collections.defaultdict(list))
    #     with open(manualDropList, 'r') as inFH:
    #         for line in inFH:
    #             line = line.split(',')
    #             line = [l.strip() for l in line]
    #             if len(line) == 1:
    #                 reason = line[0]
    #             else:
    #                 pdb = line.pop(0)
    #                 dropList[reason][pdb].extend(line)
    
    ##############################
    # Verify covalent linkages reported by PLIP
    
    
    
    #pdb = '1HJV'
    
    # Counts of different ligands at each stage of merging and filtration
    allBS = collections.Counter()
    postMergeCnt = collections.Counter()
    postGlyFilt = collections.Counter()
    postProLinkFilt = collections.Counter()
    postCentFilt = collections.Counter()
    # postContactFilt  = collections.Counter()
    
    
    # Dictionaries to store information about cleaning events for each PDB
    covMissing = collections.defaultdict(list) # dictionary listing PDB files with a missing covalent bond between glycan ligands
    caMissing = collections.defaultdict(list) # dictionary listing PDB files with a coordinating Ca 2+ ion in a separate binding site, with the values as the HETATM residues within CA_THRESH distance
    proLinkMissing = collections.defaultdict(list) # dictionary listing PDB files with a missing covalent links between a glycan and the protein, values as a list of tuples of the glycan and the connected residue
    residueOverlap = collections.defaultdict(list) # dictionary listing PDB files with multiple ligands that have overlapping centroids [alpha & beta anomeric isomers listed as separate ligands]
    glycoPepChains = collections.defaultdict(list) # Dictionary to indicate PDB files in which a glycan ligand is covalently attached to a peptide chain listed as a standard protein chain
    
    ligCount = collections.defaultdict(list) # Dictionary with pdb ids as keys, and values as a list of the number of ligands (distinct binding sites) at each step [first count, post merge, glyInMol filter, N/O-link filter, overlap filter]
    
    allPLIP = {} # Dictionary to hold all revised and altered PLIP records
    
    
    # manual overrides for special cases
    
    # 1LGB glycan presented on separate protein fragment longer than 15 aa
    glycoPepChains['1LGB'].append('C')
    
    
    
    for i,pdb in enumerate(pdbIDs):
    
        if i % 100 == 0:
            print("Processing PBD file number " + str(i) + ' of ' + str(len(pdbIDs)))
        
        
        plipIN = [f for f in plipFiles if f[-8:-4] == pdb]
        pdbIN = [f for f in pdbFiles if f[-8:-4] == pdb]
    
        if (len(plipIN) != 1) or (len(pdbIN) != 1):
            print("ERROR: ambigious or missing match to provided files for PDB ID in UniLectin data: " + pdb)
        else:
            plipIN = plipIN[0]
            pdbIN = pdbIN[0]
    
        # Read in PLIP data with plipxml parser
        plipDat = plipxml.PLIPXML(plipIN)
    
        # Read in PDB file
        if pdb in cifLst:
            parser = Bio.PDB.MMCIFParser(QUIET=True)
            pdbIN = [f for f in cifFiles if f[-8:-4] == pdb][0]
            model = parser.get_structure(pdb, pdbIN)[0] # grab the first model in PDB file
        else:
            parser = Bio.PDB.PDBParser(QUIET=True)
            model = parser.get_structure(pdb, pdbIN)[0] # grab the first model in PDB file
        
        #Manual fix of PLIP records in one(few) select edge cases
        if pdb == '4X0C': # PLIP PDB-SMILES atom mapping off by one
            for lig in ['FUC:C:2', 'NDG:C:1']:
                for k,v in plipDat.bsites[lig].mappings['smiles_to_pdb'].items():
                    plipDat.bsites[lig].mappings['smiles_to_pdb'][k] = v+1
                plipDat.bsites[lig].mappings['pdb_to_smiles'] = {v:k for k,v in plipDat.bsites[lig].mappings['smiles_to_pdb'].items()}
        
        if pdb == '1QNU': # Chain D ligand included in chain B binding site & chain E ligand included with chain A binding site
            for bigLig in ['EMB:A:193', 'EMB:B:293']:
                if bigLig == 'EMB:A:193':
                    bChain = ['A','F']
                    newLig = 'EMB:E:593'
                    nChain = ['E','J']
                else:
                    bChain = ['B','G']
                    newLig = 'EMB:D:493'
                    nChain = ['D','I']
                
                plipDat.bsites[newLig] = copy.deepcopy(plipDat.bsites[bigLig]) # copy double binding site to new instance and remove data from the donor binding site
                
                bigBS = plipDat.bsites[bigLig]
                newBS = plipDat.bsites[newLig]
                
                bigBS.members = [ m for m in bigBS.members if m.split(':')[1] in bChain]
                newBS.members = [ m for m in newBS.members if m.split(':')[1] in nChain]
                newBS.longname = '-'.join([ m.split(':')[0] for m in newBS.members])
                bigBS.longname = '-'.join([ m.split(':')[0] for m in bigBS.members])
                
                newBS.bsid = newLig
                newBS.uniqueid = ':'.join([newBS.pdbid, newBS.bsid])
                newBS.position = int(newLig.split(':')[2])
                newBS.chain = newLig.split(':')[1]
                
                # separete binding site residue lists
                bsResLst = [r['aa']+':'+r['reschain']+':'+str(r['resnr']) for r in bigBS.bs_res]
                newBS.bs_res = []
                newLigResis = [ res for res in model.get_residues() if LecGly.res2plipID(res) in newBS.members ]
                for res in model.get_residues():
                    resID = LecGly.res2plipID(res)
                    if resID in bsResLst:
                        r_ind = bsResLst.index(resID)
                        d = LecGly.newMinDist(res, newLigResis)
                        if (bigBS.bs_res[r_ind]['min_dist'] - 0.1) <= d <= (bigBS.bs_res[r_ind]['min_dist'] + 0.1): # If the min dist for a bs_res matches (accounting for rounding errors) the min dist to the correct residues in the new ligand, it belongs to that pocket
                            newBS.bs_res.append(bigBS.bs_res[r_ind])
                to_del = []
                for i,r in enumerate(bigBS.bs_res):
                    if r in newBS.bs_res:
                        to_del.append(i)
                to_del.sort(reverse = True) # Reverse order to preseve indexing while deleting
                for i in to_del:
                    del bigBS.bs_res[i]
                
                bigBS.interacting_chains = list(set([ r['reschain'] for r in bigBS.bs_res if r['contact'] == True ]))
                newBS.interacting_chains = list(set([ r['reschain'] for r in newBS.bs_res if r['contact'] == True ]))
                                
                # Drop interaction records from each site
                bigResIDs = [r['aa']+':'+r['reschain']+':'+str(r['resnr']) for r in bigBS.bs_res]
                newResIDs = [r['aa']+':'+r['reschain']+':'+str(r['resnr']) for r in newBS.bs_res]
                intTypes = ['hydrophobics', 'hbonds', 'wbridges', 'sbridges', 'pi_stacks', 'pi_cations', 'halogens', 'metal_complexes']
                for intType in intTypes:
                    for bs,resLst in [(bigBS,bigResIDs), (newBS, newResIDs)]:
                        offSet = 0
                        to_del = []
                        for i,interaction in enumerate(getattr(bs,intType)):
                            resID = getattr(interaction, 'restype') + ':' + getattr(interaction, 'reschain') + ':' + str(getattr(interaction, 'resnr') )
                            if resID not in resLst:
                                to_del.append(i)
                                offSet -= 1
                            else:
                                setattr(interaction, 'id', str(offSet + int(getattr(interaction, 'id')))) # Update interaction id number based on the number of interactions already present
                        to_del.sort(reverse = True)
                        for i in to_del:
                            del getattr(bs,intType)[i]            
                bigBS.counts = bigBS.get_counts()
                bigBS.num_contacts = bigBS.counts['total']
                bigBS.has_interactions = bigBS.num_contacts > 0
                newBS.counts = newBS.get_counts()
                newBS.num_contacts = newBS.counts['total']
                newBS.has_interactions = newBS.num_contacts > 0
                
            plipDat.num_bsites = len(plipDat.bsites)
            
        ligs = [] # List of lists of members of each PLIP binding site
        for site in plipDat.bsites:
            ligs.append(plipDat.bsites[site].members)
            allBS[plipDat.bsites[site].longname] += 1 # Count the number of occurences of each hetatm
        allLigs = [lig for lst in ligs for lig in lst]
        ligCount[pdb].append(len(ligs)) # 1st ele
        
        chainLens = {} # Number of amino acids per chain
        for chain in model:
            chainLens[chain.id] = len([res for res in chain if (res.full_id[3][0][:2] != 'H_' and res.id[0] != 'W') ]) # Number of non-water & non-HETATM residues per chain
        
        # Check for "false links" in PLIP covlink records and remove them (ligands recorded as covalently bonded but stored as separate ligands)
        ## Otherwise, the merge filter step will assume they are compostie ligands
        linksToDrop = set()
        for l1,l2 in plipDat.covlinks:
            for lig in plipDat.bsites:
                siteMembers = plipDat.bsites[lig].members
                if ((l1 in siteMembers) != (l2 in siteMembers)) and l1 in allLigs and l2 in allLigs: # If both linked resis are ligands, and one of the linked resis is in a binding site but the other isn't, remove the link
                    linksToDrop.add((l1,l2))
        if linksToDrop:
            for li in linksToDrop:
                plipDat.covlinks.remove(li)
        
    
        #### Flag potential issues in PLIP interaction profile
    
        searcher = Bio.PDB.NeighborSearch([a for a in model.get_atoms() if a.get_parent().full_id[3][0][:2] == 'H_']) # Searcher with all hetatms (no waters)
        proSearcher = Bio.PDB.NeighborSearch([a for a in model.get_atoms() if (a.get_parent().id[0] != 'W' and a.get_parent().full_id[3][0][:2] != 'H_' and a._assign_atom_mass() > 2.5)]) # Searcher with all heavy protein atoms (no waters, hydrogens/deuts [assume MM < 2.5], or HETATMs)
    
        linkPairs = [] # Hold pairs of hetatoms from PDB within at the threshold distance
        linkCentDists = {} # Store the distances between centroids for pairs of hetatms within the threshold distance (looking for overlapping anomers)
        proteinLinks = [] # Stores tuples of potential N/O linked glycans stored as ligands, and a list of bioPDB residues they are within the threshold distance of
        resiCentroids = [] # Stores tuples of potential centroid overlaps among HETATMs
    
        for res in model.get_residues():
            if res.full_id[3][0][:2] == 'H_': # For all HETATM residues in structure
                pID1 = LecGly.res2plipID(res) # unique residue ID following PLIP format (hetid:chain:resnum)
                bioPDBresi = res
                if pID1.split(':')[0] == 'CA': # Check for coordinating calciums
                    nbrs = searcher.search(res['CA'].coord, CA_THRESH, level = 'R')
                    if nbrs:
                        nbrIDs = [LecGly.res2plipID(nbr) for nbr in nbrs if (LecGly.res2plipID(nbr) in allLigs)] # Keep all neighboring residues within the threshold distance if they are implicated in a PLIP binding site
                        for bsMembers in ligs:
                            if pID1 in bsMembers and not all(n in bsMembers for n in nbrIDs): # If the Ca2+ ion is in a PLIP binding site but the neighboring residues aren't --> MISSING
                                caMissing[pdb].append(nbrIDs)
    
                else: # Check for residue overlap & missing covalent links
                    pID1cent = LecGly.resCentroid(res)
                    nbrs = [searcher.search(atm.coord, COV_THRESH, level = 'R') for atm in res.get_atoms()]
                    if nbrs:
                        nbrs = [nbr for lst in nbrs for nbr in lst if (nbr != res)] # Collapse nested list into single list
                        nbrIDs = [LecGly.res2plipID(nbr) for nbr in nbrs] # Generate PLIP IDs
                        for n,nID in zip(nbrs,nbrIDs):
                            linkCentDists[(pID1,nID)] = linkCentDists[(nID,pID1)] = LecGly.eucDist(pID1cent, LecGly.resCentroid(n))
                            if linkCentDists[(pID1,nID)] > CENTROID_THRESH: # Exclude potential missing bonds if the centroids of the residues are very close
                                if (pID1,nID) not in linkPairs:
                                    linkPairs.append((pID1,nID)) # Add new covalent links to linkPairs, mirrored incase flipped in PLIP output
                                    linkPairs.append((nID,pID1))
    
                            elif (pID1 in allLigs) and (nID in allLigs): # Potential anomeric isomer overlap between ligands
                                lig1 = [bsLig for bsLig in plipDat.bsites if pID1 in plipDat.bsites[bsLig].members][0] # Get the parent ligand ID of each
                                lig2 = [bsLig for bsLig in plipDat.bsites if nID in plipDat.bsites[bsLig].members][0]
                                try:
                                    if LecGly.glyInMolID(plipDat.bsites[lig1].longname, smilesDict) and LecGly.glyInMolID(plipDat.bsites[lig2].longname, smilesDict): # Ligands both meet minimal criteria for glycan
                                        if ((pID1, nID) not in resiCentroids) and ((nID, pID1) not in resiCentroids): # Not stored yet
                                            if lig1 != lig2:
                                                resiCentroids.append((pID1, nID))
                                except:
                                    print('***Warning: glyInMolID function error for ' + lig1 + ' & ' + lig2 + ' in pdb ' + pdb + '***\n\tDefaulting to using PLIP-provided SMILES string for glycan content check.')
                                    if LecGly.glyInMol(plipDat.bsites[lig1].smiles) and LecGly.glyInMol(plipDat.bsites[lig2].smiles): # Ligands both meet minimal criteria for glycan
                                        if ((pID1, nID) not in resiCentroids) and ((nID, pID1) not in resiCentroids): # Not stored yet
                                            if lig1 != lig2:
                                                resiCentroids.append((pID1, nID))
    
                    # Check for potential protein links
                    if pID1 in allLigs:
                        bsLig = [bsLig for bsLig in plipDat.bsites if pID1 in plipDat.bsites[bsLig].members][0]
                        try:
                            if LecGly.glyInMolID(plipDat.bsites[bsLig].longname, smilesDict):
                                proNbrs = [proSearcher.search(atm.coord, COV_THRESH, level = 'R') for atm in res.get_atoms() if atm._assign_atom_mass() > 2.5]
                                proNbrs = [nbr for lst in proNbrs for nbr in lst]
                                if proNbrs:
                                    if len(proNbrs) == 1:
                                        proteinLinks.append((pID1,proNbrs[0]))
                                    elif len(proNbrs) < 1:
                                        closestNbr,dist = LecGly.findClosestResPair(res, proNbrs)
                                        proteinLinks.append((pID1,closestNbr))
                        except:
                            print('***Warning: glyInMolID function error for ' + bsLig + ' in pdb ' + pdb + '***\n\tDefaulting to using PLIP-provided SMILES string for glycan content check.')
                            if LecGly.glyInMol(plipDat.bsites[bsLig].smiles):
                                proNbrs = [proSearcher.search(atm.coord, COV_THRESH, level = 'R') for atm in res.get_atoms() if atm._assign_atom_mass() > 2.5]
                                proNbrs = [nbr for lst in proNbrs for nbr in lst]
                                if proNbrs:
                                    if len(proNbrs) == 1:
                                        proteinLinks.append((pID1,proNbrs[0]))
                                    elif len(proNbrs) < 1:
                                        closestNbr,dist = LecGly.findClosestResPair(res, proNbrs)
                                        proteinLinks.append((pID1,closestNbr))
    
        # Filter down potential missing covalent links between ligands to confirmed missing links between ligands with glycans
        silent = [] # Store recorded missing linkages in both orders to not repeat
        missing = {} # Store a single instance of each missing linkage, as well as the indices of the atoms to join
        for link in linkPairs:
            if (link not in plipDat.covlinks) and (link[::-1] not in plipDat.covlinks): # If the pair (in either order) is noted as bonded in PLIP, don't add
                if link not in silent: # Check if already recorded
                    silent.append(link); silent.append(link[::-1])
                    if link[0] in allLigs and link[1] in allLigs: # If PLIP doesn't recognize the ligand or it's not bound, there's no binding site info to merge, skip it
                        addFlag = True # Flag variable, stores boolean of if link should be added to missing based on glycan status of each residue
                        for l in link:
                            linkRes = [bsLig for bsLig in plipDat.bsites if l in plipDat.bsites[bsLig].members][0] # Find the binding site identifying residue for a link residue if it is a member of a composite ligand (returns orignal ID if not composite)
                            if linkRes in plipDat.bsites:
                                try:
                                    if not LecGly.glyInMolID(plipDat.bsites[linkRes].longname, smilesDict):
                                        addFlag = False
                                except:
                                    print('***Warning: glyInMolID function error for ' + linkRes + ' in pdb ' + pdb + '***\n\tDefaulting to using PLIP-provided SMILES string for glycan content check.')
                                    if not LecGly.glyInMol(plipDat.bsites[linkRes].smiles):
                                        addFlag = False
                        if addFlag == True:
                            missing[link] = LecGly.findClosestAtmPair(link, biopdbModel = model, thresh = COV_THRESH)
                            covMissing[pdb].append(link)
    
    
        #### Merge information for composite ligands split across binding sites
        if missing or (pdb in caMissing.keys()):
            # Find overlap between missing covalent bpnds between mutliple ligands
            toMerge = [] # List of lists of PLIP binding sites to merge
            bsIDs2links = {} # Connect binding sites to merge with the actual residues missing covalent bonds, assumes one missing bond per two binding sites
            for link in missing.keys():
                l1resID = [bsLig for bsLig in plipDat.bsites if link[0] in plipDat.bsites[bsLig].members][0] # Get the parent ligand ID
                l2resID = [bsLig for bsLig in plipDat.bsites if link[1] in plipDat.bsites[bsLig].members][0]
                linkParents = (l1resID, l2resID) #  the residues used by PLIP to identify the binding sites
                bsIDs2links[linkParents] = link
                for mergeLst in toMerge:
                    if any(l in pair for pair in mergeLst for l in linkParents): # If either binding site is already indicated to be merged, group this
                        mergeLst.append(linkParents)
                if not any(linkParents in lst for lst in toMerge): # If the either binding site is not involved in another other mergers yet
                    toMerge.append([linkParents])
        
            # Add missing calcium ions involved in coordinating a ligand
            if pdb in caMissing.keys():
                for link in caMissing[pdb]:
                    l1resID = [bsLig for bsLig in plipDat.bsites if link[0] in plipDat.bsites[bsLig].members][0] # Get the parent ligand ID
                    l2resID = [bsLig for bsLig in plipDat.bsites if link[1] in plipDat.bsites[bsLig].members][0]
                    linkParents = (l1resID, l2resID) #  the residues used by PLIP to identify the binding sites
                    for mergeLst in toMerge:
                        if any(l in pair for pair in mergeLst for l in linkParents): # If either binding site is already indicated to be merged, group this
                            mergeLst.append(linkParents)
                    if not any(linkParents in lst for lst in toMerge): # If the either binding site is not involved in another other mergers yet
                        toMerge.append([linkParents])
        
            # Go through and integrate info from the merged binding sites
            for mergeLst in toMerge: 
                # Join non-calcium ligands together in RDkit to generate new smiles, inchikey, and logp
                newSmiles = ''
                newInchikey = ''
                molsToLink = list(set(bs for l in mergeLst for bs in l if bs.split(':')[0] != 'CA'))
                if len(molsToLink) > 1: # Would be only one if only merger invloves a CA
                    bsToLink = [plipDat.bsites[bs] for bs in molsToLink]
                    mols = []
                    for bs in bsToLink: # Create rdkit molecules, sanitized if possible
                        mols.append(Chem.MolFromSmiles(bs.smiles, sanitize = False))
                    atmIdxAdj = [0] # Number to use in adjusting the smiles index for the combined molecule (after combining, the smiles indices pick up where previous index count stops)
                    atmIdxAdj.extend([a for a in itertools.accumulate([m.GetNumAtoms() for m in mols[:-1]])]) # Get the running sum of atoms that came before for each molecule (first molecule keeps original indices)
                    atmIdxAdj = [a-1 for a in atmIdxAdj] # SMILES indices not 0-indexed in map from PLIP (1-indexed)
        
                    mergeMol = mols[0] # Recursively build combined molecule (function only accepts 2 arguements)
                    for m in mols[1:]:
                        mergeMol = Chem.CombineMols(mergeMol, m)
        
                    linksToBuild = []
                    for pair in mergeLst:
                        pdbIdxs = missing[bsIDs2links[pair]] # Get the PDB atom numbers for the atoms to link
                        linksToBuild.append(tuple(plipDat.bsites[pair[i]].mappings['pdb_to_smiles'][pdbIdxs[i]] + atmIdxAdj[molsToLink.index(pair[i])] for i in range(2))) # Get the smiles index from the PLIP mapping and adjust by corresponding number from atmIdxAdj
        
                    edMerge = Chem.EditableMol(mergeMol)
                    for l in linksToBuild:
                        edMerge.AddBond(l[0], l[1], order = Chem.rdchem.BondType.SINGLE) # Add single covalent bonds between each pair of atoms
                    mergeMol = edMerge.GetMol()
                    # Chem.FastFindRings(mergeMol) # Clean up rings to save image
                    # mergeMol
        
                    newSmiles = Chem.MolToSmiles(mergeMol)
                    newInchikey = Chem.inchi.MolToInchiKey(mergeMol)
                    try:
                        newLogP = Crippen.MolLogP(mergeMol)
                    except:
                        newLogP = ''
                        print("Warning: logP calculation error for PDB ID " + pdb)
                        
                
                for link in mergeLst:
                    plipDat.covlinks.append(link) # Add missing covalent bonds to list of linkages
        
                bsIDs = list(set([bs for tup in mergeLst for bs in tup])) # Collapse list to tuples to a list of unique res
                bsIDs = [bsid for _,bsid in sorted(zip([bsid.split(':')[1:] for bsid in bsIDs], bsIDs))] # Sorted PLIP binding site IDs by chain membership & residue number
                accID = bsIDs.pop(0) # Pop off the binding site to serve as the acceptor for the information from the others
        
                accBS = plipDat.bsites[accID] # PLIP binding site to hold all merged information
        
                accBS.composite = True
        
                for bs in bsIDs: # merge all ligand info
                    accBS.longname += '-' + plipDat.bsites[bs].longname
                    if newSmiles: # If molecules were merged in RDkit (ligands to merge other than Ca2+)
                        accBS.smiles = newSmiles
                        accBS.inchikey = newInchikey
                        accBS.logp = newLogP # Note: not quite the same as logp values generated in PLIP report, MolLogP function acting buggy,
                    accBS.members.extend(plipDat.bsites[bs].members)
                    
                    # Sum ligand properties
                    ligProps = ['heavy_atoms', 'hbd', 'unpaired_hbd', 'hba', 'unpaired_hba', 'hal', 'unpaired_hal', 'molweight', 'rotatable_bonds'] # note: rotatable bond count ignores newly introduced bonds
                    for prop in ligProps:
                        setattr(accBS, prop, getattr(plipDat.bsites[bs], prop) + getattr(accBS, prop)) # accBS.prop += plipDat.bsites[bs].prop
                    
                    # Update interactions
                    accBS.interacting_chains.extend([chain for chain in plipDat.bsites[bs].interacting_chains if chain not in accBS.interacting_chains]) # Add new interacting chains
                    intTypes = ['hydrophobics', 'hbonds', 'wbridges', 'sbridges', 'pi_stacks', 'pi_cations', 'halogens', 'metal_complexes']
                    for intType in intTypes:
                        maxAccId = len(getattr(accBS, intType)) # Find the last previous interaction id number
                        for interaction in getattr(plipDat.bsites[bs],intType):
                            setattr(interaction, 'id', str(maxAccId + int(getattr(interaction, 'id')))) # Update interaction id number based on the number of interactions already present
                            getattr(accBS, intType).append(interaction)
                    accBS.num_contacts += plipDat.bsites[bs].num_contacts
                    accBS.has_interactions = accBS.num_contacts > 0
                    accBS.counts = accBS.get_counts()
                    
                    # Update binding site residues
                    accResis = [resi['aa'] + ':' + resi['reschain'] + ':' + str(resi['resnr']) for resi in accBS.bs_res]
                    for resi in plipDat.bsites[bs].bs_res:
                        resiID = resi['aa'] + ':' + resi['reschain'] + ':' + str(resi['resnr'])
                        if resiID not in accResis: # New residue, append information
                            accBS.bs_res.append(resi)
                        else: # Update resi info, contact = True overwrites false, lowest min_dist overwrites previous
                            if resi['contact']:
                                accBS.bs_res[accResis.index(resiID)]['contact'] = True
                            if resi['min_dist'] < accBS.bs_res[accResis.index(resiID)]['min_dist']:
                                accBS.bs_res[accResis.index(resiID)]['min_dist'] = resi['min_dist']
                    
                    del plipDat.bsites[bs] # Remove redundant donor binding site
                
                plipDat.num_bsites = len(plipDat.bsites)
        
        # Update lists of ligands after merging and count occurances
        ligs = [] # List of lists of members of each PLIP binding site
        for site in plipDat.bsites:
            ligs.append(plipDat.bsites[site].members)
            postMergeCnt[plipDat.bsites[site].longname] += 1 # Count the number of occurences of each long form ligand name
        allLigs = [lig for lst in ligs for lig in lst]
        ligCount[pdb].append(len(ligs)) # 2nd ele
    
    
        #### Filter extraneous ligands
        # Check for minimal glycan content
        delLst = [] # Store binding sites to delete after the loop
        for site in plipDat.bsites:
            try:
                if not LecGly.glyInMolID(plipDat.bsites[site].longname, smilesDict):
                    delLst.append(site)
            except:
                print('***Warning: glyInMolID function error for ' + site + ' in pdb ' + pdb + '***\n\tDefaulting to using PLIP-provided SMILES string for glycan content check.')
                if not LecGly.glyInMol(plipDat.bsites[site].smiles):
                    delLst.append(site)
        for site in delLst:
            del plipDat.bsites[site]
        
        # Update lists of ligands and counters
        plipDat.num_bsites = len(plipDat.bsites) # update binding site counter
        ligs = [] # List of lists of members of each PLIP binding site
        for site in plipDat.bsites:
            ligs.append(plipDat.bsites[site].members)
            postGlyFilt[plipDat.bsites[site].longname] += 1 # Count the number of occurences of each hetatm
        allLigs = [lig for lst in ligs for lig in lst]
        ligCount[pdb].append(len(ligs)) # 3rd ele
        
        #### Filter for ligands covalently attached to the protein
        if proteinLinks:
            for proLink in proteinLinks:
                glyParent = [bsLig for bsLig in plipDat.bsites if proLink[0] in plipDat.bsites[bsLig].members][0] # get the PLIP binding site identifying residue
                res = proLink[1] # extract biopython residue to appropriately named variable
                resChain = res.get_parent().id
                if pdb in glycoPepChains:
                    if resChain in glycoPepChains[pdb]:
                        continue # Check if resChain is manually set as glycopeptide (ex: pdb 1LGB), if so, don't remove but move on to filtering interactions from glycopeptide chain
                if chainLens[resChain] > CHAIN_LEN_THRESH: # If the residue's chain is long, treat as a potential glycosylation event to remove from analysis
                    if res.get_resname() in potentialGlycosylationResis: # If the attached residue is a Ser, Thr, Asn, Trp, Cys
                        proLinkMissing[pdb].append((proLink[0], LecGly.res2plipID(res)))
                        del plipDat.bsites[glyParent]
                else: # If the chain length is below the threshold, treat glycan as a ligand attached to a peptide fragment, keep ligand but remove interactions with short chain from all binding sites in PLIP record
                    if resChain not in glycoPepChains[pdb]:
                        glycoPepChains[pdb].append(resChain)
                        
            if pdb in glycoPepChains.keys(): # If there are glycopeptides present, remove all interactions with them from the PLIP binding sites
                for bs in plipDat.bsites: # for all ligands
                    
                    glyBS = plipDat.bsites[bs]
                    
                    for resChain in glycoPepChains[pdb]: # remove all contacts and interactions with any glycopeptide chains
                        
                        #if resChain in glyBS.interacting_chains:
    
                        # Drop bs residues on resChain
                        toDrop = [i for i,bsr in enumerate(glyBS.bs_res) if bsr['reschain'] == resChain]
                        for i in sorted(toDrop, reverse=True): # Go in reverse order to not change the index of later elements
                            del glyBS.bs_res[i]
                            
                        # Drop interactions with residues on resChain
                        if resChain in glyBS.interacting_chains:
                            glyBS.interacting_chains.remove(resChain)
                            
                            intTypes = ['hydrophobics', 'hbonds', 'wbridges', 'sbridges', 'pi_stacks', 'pi_cations', 'halogens', 'metal_complexes']
                            for intType in intTypes:
                                toDrop = [i for i,interaction in enumerate(getattr(glyBS, intType)) if interaction.reschain == resChain] # Find interactions with residues from the resChain
                                for i in sorted(toDrop, reverse=True):
                                    del getattr(glyBS, intType)[i] # Delete identified residues is reverse order to preserve indices
                                for i,interaction in enumerate(getattr(glyBS, intType)): # renumber existing interactions in increasing intergers from 1
                                    interaction.id = i + 1
                                    
                            glyBS.has_interactions = glyBS.num_contacts > 0 # Set to false if all interactions removed
                            glyBS.counts = glyBS.get_counts() # Update counts table
                            glyBS.num_contacts = glyBS.counts['total'] # update number of interactions
            
        # Update lists of ligands and counters
        plipDat.num_bsites = len(plipDat.bsites) # update binding site counter
        ligs = [] # List of lists of members of each PLIP binding site
        for site in plipDat.bsites:
            ligs.append(plipDat.bsites[site].members)
            postProLinkFilt[plipDat.bsites[site].longname] += 1 # Count the number of occurences of each hetatm
        allLigs = [lig for lst in ligs for lig in lst]
        ligCount[pdb].append(len(ligs)) # 4th ele
        
        #### Filter for centroid overlap
        if resiCentroids:
            residueOverlap[pdb].extend(resiCentroids)
            
            # pick one ligand to keep in each case
            for (r1,r2) in resiCentroids:
                lig1 = plipDat.bsites[r1]; lig2 = plipDat.bsites[r2]
                if lig1.composite and not lig2.composite: # Keep composite ligand over non-composite
                    drop = r2
                elif not lig1.composite and  lig2.composite:
                    drop = r1
                else:
                    if lig1.counts['total'] > lig2.counts['total']: # drop ligand with most interactions
                        drop = r2
                    elif lig1.counts['total'] < lig2.counts['total']:
                        drop = r1
                    else:
                        if len(lig1.bs_res) > len(lig2.bs_res): # drop ligand with most interactions
                            drop = r2
                        elif len(lig1.bs_res) < len(lig2.bs_res):
                            drop = r1
                        else: # otherwise, pick one randomly (seems to be the case for PLIP usually)
                            drop = r2
                            
                del plipDat.bsites[drop]
                
            plipDat.num_bsites = len(plipDat.bsites)
            
        # Update lists of ligands and counters
        plipDat.num_bsites = len(plipDat.bsites) # update binding site counter
        ligs = [] # List of lists of members of each PLIP binding site
        for site in plipDat.bsites:
            ligs.append(plipDat.bsites[site].members)
            postCentFilt[plipDat.bsites[site].longname] += 1 # Count the number of occurences of each hetatm
        allLigs = [lig for lst in ligs for lig in lst]
        ligCount[pdb].append(len(ligs)) # 5th ele
        
        
        # #### Filter out ligands from manual check of ligands with minimal interactions
        # if manualFilter:
        #     if pdb in dropList['Missed_glycosylation'].keys():
        #         for bs in dropList['Missed_glycosylation'][pdb]:
        #             del plipDat.bsites[bs]
            
        #     if pdb in dropList['Missed_covalent_link'].keys():
        #         for bs in dropList['Missed_covalent_link'][pdb]:
        #             del plipDat.bsites[bs]
        #         plipDat.num_bsites = len(plipDat.bsites)        
            
        #     if pdb in dropList['Other'].keys():
        #         for bs in dropList['Other'][pdb]:
        #             del plipDat.bsites[bs]
                    
        #     plipDat.num_bsites = len(plipDat.bsites) # update binding site counter
        
        
        
        allPLIP[pdb] = plipDat
    
    ########## Generate QC report in ./analysis/QC
    
     # Highlight cases where all ligands are removed
    flagCases = []
    
    initCnts = []
    mergeCnts = []
    glyCnts = []
    prolinkCnts = []
    centroidCnts = []
    for pdb in ligCount:
        initCnts.append(ligCount[pdb][0])
        mergeCnts.append(ligCount[pdb][1])
        glyCnts.append(ligCount[pdb][2])
        prolinkCnts.append(ligCount[pdb][3])
        centroidCnts.append(ligCount[pdb][4])
        if ligCount[pdb][4] == 0:
            print('WARNING: after post-processing, PDB ID', pdb, 'ended up with 0 PLIP ligands. The initial ligand count was', ligCount[pdb][0])
            flagCases.append(pdb)
    
    recordDelLst = []
    for pdb in allPLIP.keys():
        if len(allPLIP[pdb].bsites) == 0:
            recordDelLst.append(pdb)
    for pdb in recordDelLst:
        del allPLIP[pdb]
    
    if not os.path.exists(qcDir):
        os.makedirs(qcDir)
    if not os.path.exists(qcDir+'plots/'):
        os.makedirs(qcDir+'plots/')
    
    with open(file = qcDir+'ligandCheckReport.txt', mode = 'w') as qcFH:
        qcFH.write('\n#################################\nLigand checker script quality control report, generated from checkPLIPligands.py\n#################################\n\n')
        if flagCases:
            zeroByGly = []
            zeroByPro = []
            qcFH.write(str(len(flagCases)) + ' PDB files with no ligands left over after all filter steps\n\n')
            qcFH.write('Ligand counts at each step of cleaning and filtering for ligandless PDBs:\n')
            qcFH.write('pdb \t       init \t     merge \t    glycan \t    protein \t     cent\n')
            for pdb in flagCases:
                qcFH.write(pdb+' :\t\t ' + '\t\t'.join([str(cnt) for cnt in ligCount[pdb]]) + '\n')
                if ligCount[pdb][2] == 0:
                    zeroByGly.append(pdb)
                else:
                    zeroByPro.append(pdb)
                    
            qcFH.write('Zeroed out by glycan content filter:\n\t' + str(zeroByGly) + '\n')
            qcFH.write('Zeroed out by protein-linkage filter:\n\t'+ str(zeroByPro) + '\n')
            qcFH.write('\nCheck post-processing reports for each pdb to verify each cases, add to override exceptions if appropriate\n')
            qcFH.write('In PDB 4P14:\n\tLigand is adenine w/o ribose, no glycan content.\n')
            qcFH.write('In PDBs 4B15, 5Y96, 4Z8S, 4ZA3, & 4EB2:\n\tOnly glycans present are glycosylation events.\n\n\n')
    
    
    
        # Overview of ligand counts at each step of ligand filtration
        cnts = [initCnts, mergeCnts, glyCnts, prolinkCnts, centroidCnts]
        plt.boxplot(cnts,  labels = ['Initial', 'Lig. Merge', 'Gly. Check', 'Pro. Linked', 'Cent. Overlap'])
        plt.title('Distribution of ligand counts per PDB file following each filter')
        plt.savefig(qcDir+'plots/'+'ligCntBoxplots.jpg', dpi = 240)
        plt.close()
        
        qcFH.write('Total ligand counts at each step of ligand filtering, accompanies figure in ./analysis/QC/plots/ligCntBoxplots.jpg:\n')
    
        qcFH.write('Initial number of binding sites across all PDB files:\n')
        qcFH.write('\t'+str(sum(initCnts)) + ', ' + str(len(allBS)) + ' unique ligands\n')
        qcFH.write('Number of binding sites across all PDB files after merging glycan ligands with missing links and Ca2+ coord ions:\n\t'+str(sum(mergeCnts)) + ', ' + str(len(postMergeCnt)) + ' unique ligands\n')
        qcFH.write('Number of binding sites across all PDB files after filtering for minimal glycan content:\n\t'+str(sum(glyCnts)) + ', ' + str(len(postGlyFilt)) + ' unique ligands\n')
        qcFH.write('Number of binding sites across all PDB files after filtering for protein glycosylation:\n\t'+str(sum(prolinkCnts)) + ', ' + str(len(postProLinkFilt)) + ' unique ligands\n')
        qcFH.write('Number of binding sites across all PDB files after filtering glycan ligands with overlapping centroids:\n\t'+str(sum(centroidCnts)) + ', ' + str(len(postCentFilt)) + ' unique ligands\n')
    
    
        # Check glycosylation occurences
        proLinkUnis = set()
        glycosylationCnts = []
        proLinkResTypes = collections.Counter()
        
        for pdb in proLinkMissing:
            glycosylationCnts.append(len(proLinkMissing[pdb]))
            for uni in uniprotMap:
                if pdb in uniprotMap[uni]:
                    proLinkUnis.add(uni)
            for link in proLinkMissing[pdb]:
                proLinkResTypes[link[1][:3]] += 1
        
        
        qcFH.write('\n\n\n{} PLIP ligands detected within covalent bond distance of a protein residue in {} PDB files from {} unique UniProt IDs.\n'.format(sum(glycosylationCnts),len(proLinkMissing), len(proLinkUnis)))
        if glycosylationCnts:
            qcFH.write('\tMedian {0} {1}\n\tMean {0} {2}'.format('glycosylation occurrences when glycosylation observed:',np.median(glycosylationCnts), round(np.mean(glycosylationCnts), 2)))
        
            plt.hist(glycosylationCnts, bins = 20, range = [0,20], edgecolor = 'black', color = 'blue') 
            plt.title('Hist of number of glycosylation occurrences per PDB file')
            plt.xlabel('Number of ligands within colvalent bond distances of protein')
            plt.savefig(qcDir+'plots/'+'glycosylationBarplot.jpg', dpi = 240)
            plt.close()
            
            qcFH.write('\n\nThe closest types of protein residue when a PLIP ligand glycan is within {} Angstroms, along with the number of occurrences\n\t'.format(COV_THRESH))
            qcFH.write(str(proLinkResTypes.most_common())+'\n')
    
        # Brief stats on other ligand filter steps
        qcFH.write('\n\n\n'+str(len(covMissing)) + ' pdb files with missing covalent bonds in ligands detected & fixed:\n')
        for k in covMissing:
            qcFH.write(k + '\n')
            for l in covMissing[k]:
                qcFH.write('\t'+str(l) + '\n')
        
        qcFH.write('\n\n\n'+str(len(caMissing)) + ' pdb files with coordination Ca2+ ions missing from the compound ligand:\n')
        for k in caMissing:
            qcFH.write(k + '\n')
            for l in caMissing[k]:
                qcFH.write('\t'+str(l) + '\n')
    
        qcFH.write('\n\n\n'+str(len(residueOverlap)) + ' pdb files with residues from separate ligands that have overlapping centroids:\n')
        for k in residueOverlap:
            qcFH.write(k + '\n')
            for l in residueOverlap[k]:
                qcFH.write('\t'+str(l) + '\n')
    
        qcFH.write('\n\n\n'+str(len(glycoPepChains)) + ' pdb files with glycan ligands detected on the listed glycopeptide chains.\nInteractions wih the listed chains excluded from the analysis but ligands are kept.\n\n')
        for k in glycoPepChains:
            qcFH.write(k + '\n')
            qcFH.write('\t'+str(glycoPepChains[k]) + '\n')
    
        # Bar plot of ligand frequencies at each step of filtering for all unique ligands within the top 10 most common from any step
        topLigs = set()
        for k,v in allBS.most_common(10): topLigs.add(k)
        for k,v in postMergeCnt.most_common(10): topLigs.add(k)
        for k,v in postGlyFilt.most_common(10): topLigs.add(k)
        for k,v in postProLinkFilt.most_common(10): topLigs.add(k)
        for k,v in postCentFilt.most_common(10): topLigs.add(k)
        
        labels = topLigs
        allbs = [allBS[l] if l in allBS.keys() else 0 for l in topLigs]
        postmerge = [postMergeCnt[l] if l in postMergeCnt.keys() else 0 for l in topLigs]
        postgly = [postGlyFilt[l] if l in postGlyFilt.keys() else 0 for l in topLigs]
        postpro = [postProLinkFilt[l] if l in postProLinkFilt.keys() else 0 for l in topLigs]
        postcent = [postCentFilt[l] if l in postCentFilt.keys() else 0 for l in topLigs]
        
        x = np.arange(len(labels))  # the label locations
        width = 0.15  # the width of the bars
        fig, ax = plt.subplots()
        
        rects1 = ax.bar(x,
                        allbs,
                        width,
                        label='allBS')
        rects2 = ax.bar(x + width,
                        postmerge,
                        width,
                        label='mergeLigs')
        rects3 = ax.bar(x + 2*width,
                        postgly, 
                        width, 
                        label='glyFilter')
        rects4 = ax.bar(x + 3*width, 
                        postpro, 
                        width, 
                        label='proteinLinked')
        rects5 = ax.bar(x + 4*width, 
                        postcent, 
                        width, 
                        label='centOverlap')
        # Add some text for labels, title and custom x-axis tick labels, etc.
        ax.set_ylabel('Ligand count')
        ax.set_title('Ligand counts across all PDBs by filter step')
        ax.set_xticks([p + 1.5*width for p in x])
        ax.set_xticklabels(labels, rotation='vertical')
        ax.legend()
        fig.tight_layout()
        plt.savefig(qcDir+'plots/'+'commonLigandCounts.jpg', dpi = 240)
        plt.close()
        
        pdbCnt = len(allPLIP)
        bsCnt = 0
        for pdb in allPLIP.keys():
            for bs in allPLIP[pdb].bsites:
                bsCnt += 1
        qcFH.write('\n#####\nThere are ' + str(bsCnt) + ' binding sites from ' + str(pdbCnt) + ' structures\n#####\n')
        
    
    
    ##################  Expand binding sites, add dssp info, look for hetatms near ligand
    
    #pdb = '3LL2'
    
    
    dsspErrs = collections.defaultdict(list)
    bsresCnt = 0
    
    negIndDict = {4294967295: -1, # PLIP can't handle negative residue numbers, this maps the values back to the original index value
                  4294967294: -2} 
    
    renumPDBLst = ['1RV0', '1RVT', '3HTP', '3HTT', '3M5H', '3M5I', '3M5J', '3UBE', '3UBJ', '3UBN', '3UBQ'] # PDB IDs that give a key error from the dssp data for a zero numbered residue in the binding site
    zeroResis = collections.defaultdict(lambda: collections.defaultdict(dict)) # nested default dicts to store 0-numbered resi occurences by PDB, then by binding site
    for pdb in renumPDBLst:
        for bs in allPLIP[pdb].bsites:
                for res in allPLIP[pdb].bsites[bs].bs_res:
                    if res['resnr'] == 0:
                        zeroResis[pdb][bs][str(res['resnr']) + res['reschain']] = ''
                        
    
    zeroResis['1RV0']['DAN:L:701']['0L'] = '133A' # Manually map identities of 0-numbered resis to original numbering (Assumes same chain)
    
    zeroResis['1RVT']['GAL:A:1']['0H'] = '133A'
    zeroResis['1RVT']['BGC:B:1']['0J'] = '133A'
    zeroResis['1RVT']['GAL:C:1']['0L'] = '133A'
    
    zeroResis['3HTP']['NAG:D:1']['0A'] = '132A'
    
    zeroResis['3HTT']['GLC:D:1']['0A'] = '132A'
    
    zeroResis['3M5H']['NAG:G:1']['0A'] = '158A'
    zeroResis['3M5H']['GAL:I:1']['0C'] = '158A'
    zeroResis['3M5H']['NAG:K:1']['0E'] = '158A'
    
    zeroResis['3M5I']['GAL:G:1']['0A'] = '158A'
    zeroResis['3M5I']['GAL:I:1']['0C'] = '158A'
    zeroResis['3M5I']['GAL:L:1']['0E'] = '158A'
    
    zeroResis['3M5J']['GAL:G:1']['0A'] = '158A'
    zeroResis['3M5J']['NAG:H:1']['0C'] = '158A'
    zeroResis['3M5J']['NAG:I:1']['0E'] = '158A'
    
    zeroResis['3UBE']['GAL:M:1']['0A'] = '133A'
    zeroResis['3UBE']['GAL:N:1']['0C'] = '133A'
    zeroResis['3UBE']['GAL:N:1']['0K'] = '116C'
    zeroResis['3UBE']['NAG:O:1']['0E'] = '133A'
    zeroResis['3UBE']['SIA:G:401']['0G'] = '133A'
    zeroResis['3UBE']['NAG:Q:1']['0I'] = '133A'
    zeroResis['3UBE']['GAL:R:1']['0K'] = '133A'
    
    zeroResis['3UBJ']['GAL:N:1']['0C'] = '133A'
    zeroResis['3UBJ']['GAL:P:1']['0G'] = '133A'
    zeroResis['3UBJ']['GAL:Q:1']['0I'] = '133A'
    zeroResis['3UBJ']['GAL:Q:1']['0E'] = '116C'
    zeroResis['3UBJ']['GAL:R:1']['0K'] = '133A'
    
    zeroResis['3UBN']['NAG:M:1']['0A'] = '133A'
    zeroResis['3UBN']['NAG:N:1']['0C'] = '133A'
    zeroResis['3UBN']['NAG:N:1']['0K'] = '116C'
    zeroResis['3UBN']['NAG:O:1']['0E'] = '133A'
    zeroResis['3UBN']['GAL:P:1']['0G'] = '133A'
    zeroResis['3UBN']['NAG:R:1']['0I'] = '133A'
    zeroResis['3UBN']['NAG:S:1']['0K'] = '133A'
    
    zeroResis['3UBQ']['NAG:N:1']['0C'] = '133A'
    zeroResis['3UBQ']['GAL:P:1']['0I'] = '133A'
    zeroResis['3UBQ']['GAL:P:1']['0E'] = '116C'
    zeroResis['3UBQ']['GAL:Q:1']['0K'] = '133A'
    
        
    for pdb in allPLIP.keys():
        
        pdbIN = [f for f in pdbFiles if f[-8:-4] == pdb][0]
        
        
        if pdb in cifLst:
            parser = Bio.PDB.MMCIFParser(QUIET=True)
            pdbIN = [f for f in cifFiles if f[-8:-4] == pdb][0]
            model = parser.get_structure(pdb, pdbIN)[0] # grab the first model in PDB file
        else:
            parser = Bio.PDB.PDBParser(QUIET=True)
            model = parser.get_structure(pdb, pdbIN)[0] # grab the first model in PDB file
        
        searcher = Bio.PDB.NeighborSearch([a for a in model.get_atoms() if a.get_parent().full_id[3][0][:2] == 'H_']) # Searcher with all hetatms (no waters)
        
        # chainLens = {} # Number of amino acids per chain
        # for chain in model:
        #     chainLens[chain.id] = len([res for res in chain if (res.full_id[3][0][:2] != 'H_' and res.id[0] != 'W') ]) # Number of non-water & non-HETATM residues per chain
        
        # Expand PLIP binding site residues
        
        for bs in allPLIP[pdb].bsites: # expand binding sites to include residues +/- 2 residues from interaction resdiues (ie xxXXIXXxx, where X = BS residues, I = interacting residue, and x = non-BS residues)
               
            intRes = [] # Store all residues with interactions with ligand
            bsResis = [] # Store all residues in binding site regardless of interactions
            ligResis = [res for res in model.get_residues() if LecGly.res2plipID(res) in allPLIP[pdb].bsites[bs].members]
            for res in allPLIP[pdb].bsites[bs].bs_res:
                resID = str(res['resnr']) + res['reschain']
                bsResis.append(resID)
                if res['contact']: intRes.append(resID)
            for resID in intRes: # expand binding sites to include residues +/- 2 residues from interacting resdiues (ie xxXXIXXxx, where X = BS residues, I = interacting residue, and x = non-BS residues)
                resn = int(''.join(list(filter(str.isdigit, resID))))
                chain = resID[len(str(resn)):]
                toAdd = [resn-2, resn-1, resn+1, resn+2]
                toAdd = [str(rn)+chain for rn in toAdd if ((str(rn)+chain not in bsResis))]
                for addResID in toAdd:
                    resn = int(''.join(list(filter(str.isdigit, addResID))))
                    chain = addResID[len(str(resn)):]
                    try:
                        res = model[chain][resn] # Get bioPDD residue if it exists, otherwise skip it and move on to next one
                        aa = res.resname 
                        minDist = LecGly.newMinDist(res, ligResis)
                        minDist = round(minDist, 1)
                        allPLIP[pdb].bsites[bs].bs_res.append({'resnr': resn, 'reschain': chain, 'aa': aa, 'contact': False, 'min_dist': minDist})
                        bsResis.append(addResID)
                    except:
                        continue
            
            allBSsplit = [(int(''.join(list(filter(str.isdigit, resID)))), resID[len(list(filter(str.isdigit, resID))):]) for resID in bsResis] # create list of tuples of binding site residues with (resnr, chain)
            toAdd = []
            for resID in allBSsplit: # expand bs resis to include residues that are between two bs resis in sequence space
                if (resID[0]+2, resID[1]) in allBSsplit and (resID[0]+1, resID[1]) not in allBSsplit:
                   toAdd.append(str(resID[0]+1)+resID[1])
            for resID in toAdd:
                resn = int(''.join(list(filter(str.isdigit, resID))))
                chain = resID[len(str(resn)):]
                try:
                    res = model[chain][resn] # Get bioPDD residue if it exists, otherwise skip it and move on to next one
                    aa = res.resname 
                    minDist = LecGly.newMinDist(res, ligResis)
                    minDist = round(minDist, 1)
                    allPLIP[pdb].bsites[bs].bs_res.append({'resnr': resn, 'reschain': chain, 'aa': aa, 'contact': False, 'min_dist': minDist})
                    bsResis.append(resID)
                except:
                    continue
                
    
        # Add DSSP info
        dsspIN = [f for f in dsspFiles if f[-9:-5] == pdb]
    
        if (len(dsspIN) != 1):
            print("ERROR: ambigious or missing match to dssp file for PDB ID in UniLectin data: " + pdb)
            for bs in allPLIP[pdb].bsites:
                for res in allPLIP[pdb].bsites[bs].bs_res: # For all residues in the PLIP binding site, add missing DSSP info placeholders
                    res['struct'] = ''; res['acc'] = ''; res['phi'] = ''; res['psi'] = '.'
        else:
            dsspIN = dsspIN[0]
            dsspDat = make_dssp_dict(dsspIN) # Parse DSSP file with Bio.PDB DSSP parser
            dsspDat = LecGly.reformatDSSP(dsspDat[0]) # extract dictionary from tuple and reformat key
            toFix = set() # List of residues that are overwritten in the dsspDat dict by the dssp parser
            
            for bs in allPLIP[pdb].bsites:
                for res in allPLIP[pdb].bsites[bs].bs_res: # For all residues in the PLIP binding site, update with struct, acc, phi, and psi from DSSP
                    bsresCnt += 1
                    if res['resnr'] in negIndDict.keys(): # If the residue number is the 4 billions, it's probably the neg index bug. This fixes the PLIP record
                        res['resnr'] = negIndDict[res['resnr']]
                    resID = str(res['resnr']) + res['reschain']
                    if pdb in zeroResis.keys():
                        if resID in zeroResis[pdb][bs]:
                            dsspInfo = LecGly.manualDSSPpull(dsspIN, residueID = zeroResis[pdb][bs][resID][:-1] + res['reschain'], insertionCode = zeroResis[pdb][bs][resID][-1]) # Read in single line from DSSP file matching inserted residue
                            res['struct'] = dsspInfo[0] # one letter DSSP secondary structure code
                            res['acc'] = dsspInfo[1] # solvent accesibility
                            res['phi'] = dsspInfo[2] # phi angle
                            res['psi'] = dsspInfo[3] # psi angle
                            
                            resID = zeroResis[pdb][bs][resID][:-1] + res['reschain'] # DSSP info for inserted, 0-numbered resis parsed into entry for previous resis (ie 133AL (PLIP:0L), DSSP info stored under 133L)
                            toFix.add(resID) # Mark the overwritten resi (ie 133L) to extract info directly from the DSSP file if in the binding site
                            continue
                    try:
                        res['struct'] = dsspDat[resID][1] # one letter DSSP secondary structure code
                        res['acc'] = dsspDat[resID][2] # solvent accesibility
                        res['phi'] = dsspDat[resID][3] # phi angle
                        res['psi'] = dsspDat[resID][4] # psi angle
                    except KeyError: # If the residue isn't in the DSSP record
                        dsspErrs[pdb].append(resID)
                        res['struct'] = ''; res['acc'] = ''; res['phi'] = ''; res['psi'] = ''
            for missingResi in toFix:
                for bs in allPLIP[pdb].bsites:
                    for res in allPLIP[pdb].bsites[bs].bs_res:
                        if res['resnr'] == int(missingResi[:-1]) and res['reschain'] == missingResi[-1]:
                            dsspInfo = LecGly.manualDSSPpull(dsspIN, missingResi)
                            print('Fixing',pdb,res)
                            res['struct'] = dsspInfo[0] # one letter DSSP secondary structure code
                            res['acc'] = dsspInfo[1] # solvent accesibility
                            res['phi'] = dsspInfo[2] # phi angle
                            res['psi'] = dsspInfo[3] # psi angle
                    
                    
    
    if dsspErrs:
        n=0
        for k,v in dsspErrs.items():
            for r in v:
                n += 1
        print('Warning: DSSP records missing at least one residue from ' + str(len(dsspErrs)) + ' different PDB ids')
        print('\tMissing records for ' + str(n) + ' total binding site residues')
        print('DSSP info for absent residues recorded as empty string ("")')
    
  
    
  
    ###################
    # Serialize cleaned records to file
    
    for pdb in allPLIP.keys(): # Remove original XML parsed through xml.etree. Data already extracted and can't pickle xml.etree._Element objects
        if hasattr(allPLIP[pdb], 'doc'):
            del allPLIP[pdb].doc
        for lig in allPLIP[pdb].bsites:
            if hasattr(allPLIP[pdb].bsites[lig], 'bindingsite'):
                del allPLIP[pdb].bsites[lig].bindingsite
    
    with open("./data/unilectin/PLIPrecords.p", "wb") as pickleFH:
        dill.dump(allPLIP, pickleFH)
    
    with open("./data/unilectin/icodeResi2PLIPmapping.p", "wb")as pickleFH:
        dill.dump(zeroResis, pickleFH)
    
    # To load in again:
    # allPLIP = pickle.load(open("./data/unilectin/PLIPrecords.p", "rb"))
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
