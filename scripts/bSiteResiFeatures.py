#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 11 16:41:37 2020

@author: dmattox

Extract features for each lectin binding site from the PLIP output and PDB file and save to file in a table format
"""

import os
import glob
import collections
import dill
import numpy as np
import Bio.PDB
# from Bio.PDB.DSSP import make_dssp_dict

import lec_gly as LecGly
import checkPLIPligands as LigCheck

os.chdir(LecGly.homeDir)


########################################################################    

plipFile = './data/unilectin/cleanPLIPrecords.p'
zeroResMapping = './data/unilectin/icodeResi2PLIPmapping.p'

# pdbDir = './data/structures/holo/origPDB/'
# cifDir = './data/structures/holo/origPDB/'
unilecDataFH = './data/unilectin/boundUnilectinData.tsv'
# dsspDir = './data/dssp/dsspOut/'
outTSV = './analysis/bsiteResiFeatures.tsv'

########################################################################

pdbFiles = glob.glob(LigCheck.originalPDBsDir + '*.pdb')


if __name__ == '__main__':
    
    # Read in serialized PLIP recordspost cleaning and filtering, along with mapping from PLIP 0-numbered residues to original numbering
    with open(plipFile, "rb") as pickleFH:
        allPLIP = dill.load(pickleFH)
    
    with open(zeroResMapping, "rb")as pickleFH:
        zeroResis = dill.load(pickleFH)
    
    ########################################################################
    # Check for other ions or small molecules that might be present in the binding site, including other ligands or N/O-linked 
    
    hetAtmTracker = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.Counter()))  # [pdb][bs]: {glycan_lig: cnt, glycosylation: cnt, CA: cnt, other: cnt}  Potentially add other ions later? MG: cnt, NA: cnt}  
    
    for i,pdb in enumerate(allPLIP.keys()):       
        
        pdbIN = [f for f in pdbFiles if f[-8:-4] == pdb][0]

        if pdb in LigCheck.cifLst:
            parser = Bio.PDB.MMCIFParser(QUIET=True)
            pdbIN = [f for f in LigCheck.cifFiles if f[-8:-4] == pdb][0]
        else:
            parser = Bio.PDB.PDBParser(QUIET=True)
            
        model = parser.get_structure(pdb, pdbIN)[0] # grab the first model in PDB file
        
        searcher = Bio.PDB.NeighborSearch([a for a in model.get_atoms() if a.get_parent().full_id[3][0][:2] == 'H_']) # Searcher with all hetatms (no waters
                
        ligs = [] # List of lists of members of each PLIP binding site
        for site in allPLIP[pdb].bsites:
            ligs.append(allPLIP[pdb].bsites[site].members)
        allLigs = [lig for lst in ligs for lig in lst]
        
        for bs in allPLIP[pdb].bsites:
            # Look for neighobring hetAtms
            
            ligAtms = LecGly.getLigAtoms(model, allPLIP[pdb].bsites[bs])
            ligAtms = [la for la in ligAtms if la.element not in ['CA','MN','MG','NA']]
            neighbors = [ n for lst in [searcher.search(atm.coord, LigCheck.CA_THRESH, level = 'R') for atm in ligAtms] for n in lst ] # return a flattened list of neighboring, non-water het atm residues
            neighbors = [LecGly.res2plipID(nbr) for nbr in  list(set(neighbors))]
            
            hetAtmTracker[pdb][bs]['CA'] += len([nbrID for nbrID in neighbors if nbrID.split(':')[0] == 'CA']) # Count ions before dropping self hetatms
            hetAtmTracker[pdb][bs]['MN'] += len([nbrID for nbrID in neighbors if nbrID.split(':')[0] == 'MN']) # Count ions before dropping self hetatms
            hetAtmTracker[pdb][bs]['NA'] += len([nbrID for nbrID in neighbors if nbrID.split(':')[0] == 'NA']) # Count ions before dropping self hetatms
            hetAtmTracker[pdb][bs]['MG'] += len([nbrID for nbrID in neighbors if nbrID.split(':')[0] == 'MG']) # Count ions before dropping self hetatms
            
            neighbors = [n for n in neighbors if n not in allPLIP[pdb].bsites[bs].members]
            
            if neighbors:
                for nbr in neighbors:
                    nbr_het_id = nbr.split(':')[0]
                    if nbr in allLigs: # if the neighboring het atom is a member of another glycan binding site
                        hetAtmTracker[pdb][bs]['glycan_lig'] = 1
                    elif not LecGly.glyInMolID(nbr_het_id, LigCheck.smilesDict): # exclude anomeric isomers (alpha & beta glycans on top of each other)
                        hetAtmTracker[pdb][bs]['other'] += 1
    
 
    
    #####################
    # Save out features of bindingsites a dataframe (tsv)
    
    dsspCodes = ['H', 'B', 'E', 'G', 'I', 'T', 'S', '-']
    aaGroups = ['nonpolar','polar','posCharge','negCharge','aromatic']
    res2group = {
        'GLY': 'nonpolar',
        'ALA': 'nonpolar',
        'VAL': 'nonpolar',
        'LEU': 'nonpolar',
        'MET': 'nonpolar',
        'MSE': 'nonpolar',
        'ILE': 'nonpolar',
        'SER': 'polar',
        'THR': 'polar',
        'CYS': 'polar',
        'PRO': 'polar',
        'ASN': 'polar',
        'GLN': 'polar',
        'LYS': 'posCharge',
        'ARG': 'posCharge',
        'HIS': 'posCharge',
        'ASP': 'negCharge',
        'GLU': 'negCharge',
        'PHE': 'aromatic',
        'TYR': 'aromatic',
        'TRP': 'aromatic'}
    resis = list(res2group.keys()) + ['ALT'] # Count non-standard amino acids as "ALT" for alternate
    intNames = ['hydrophobics', 'hbonds', 'wbridges', 'sbridges', 'pistacks', 'pications', 'halogens', 'metal', 'hbond_back', 'hbond_nonback', 'total']
    hetAtmResiCnts = ['glycan_lig', 'CA', 'MN', 'MG', 'NA', 'other']
    
    # Binning distances
    bin1 = 3.5 # (0,3.5] 
    bin2 = 4.5 # (3.5,4.5]
    bin3 = 6.5 # (4.5,6.5]
    bin4 = 8 # (6.5,8]
    
    binColnames = dsspCodes + aaGroups + resis
    binColnames = [ele+'_bin'+str(i) for i in range(1,5) for ele in binColnames]
    
    with open(outTSV, 'w') as outFH:
        outFH.write('bSite\tpdb\t' + '\t'.join(LigCheck.unilecCols) + '\tnumBsites\tlongname\tbsiteSequence\tnumBSresis_bin1\tnumBSresis_bin2\tnumBSresis_bin3\tnumBSresis_bin4\t' + '\t'.join(intNames) + '\t' + '\t'.join(binColnames) + '\t' + '\t'.join(hetAtmResiCnts) + '\n') # write column names
        for pdb in allPLIP.keys():
            plipDat = allPLIP[pdb]
            numBS = len(plipDat.bsites)
            for lig in plipDat.bsites:
                bs = plipDat.bsites[lig]
                
                # Get seqeunces of continuous pieces of protein involved in binding
                bsResis = [] # list to hold binding site residues formatted as a list of lists [chain, index, i_code, amino acid]]
                for r in bs.bs_res:
                    resID = [r['reschain'], r['resnr'], '', r['aa']]
                    if pdb in zeroResis.keys(): # Check for zero numbered resiudes and correct the res num to match the original number including the i_code
                        if lig in zeroResis[pdb]:
                            rID = str(resID[1])+resID[0] # resnum from PLIP with chain ID
                            if rID in zeroResis[pdb][lig].keys():
                                resID[1] = int(zeroResis[pdb][lig][rID][:-1]) # Replace PLIP resnum with true resnum
                                resID[2] = zeroResis[pdb][lig][rID][-1] # Add i_code
                    bsResis.append(resID)
                bsResis.sort()
                bsSeqs = LecGly.three2one[bsResis[0][3]] if bsResis[0][3] in LecGly.three2one.keys() else 'X'# stores the sequences in order, adding a split "|" when a piece of continuous seqeunce ends
                for i in range(1,len(bsResis)):
                    prevRes = bsResis[i-1]
                    if bsResis[i][3] in LecGly.three2one.keys():
                        aaLetter = LecGly.three2one[bsResis[i][3]]
                    else:
                        aaLetter = 'X' # Record non-standard amino acids as X
                    if ((bsResis[i][0] == prevRes[0]) and ((bsResis[i][1] == prevRes[1]) or (bsResis[i][1] == prevRes[1]+ 1))):
                        bsSeqs += aaLetter
                    else:
                        bsSeqs += '|' + aaLetter
                
                # Store bin-level information about id, type, & secondary struct. of each residue in the binding site
                bsResiCnts = [0]*4 # Hold the number of residues in each bin, 1st ele corresponds to 1st bin, 2nd -> 2nd bin, ...
                bin1stats = collections.Counter() # Hold the number of each type of ss or res type occurences in each bin, as well as aa identity
                bin2stats = collections.Counter()
                bin3stats = collections.Counter()
                bin4stats = collections.Counter()
                
                bin1Resis = collections.defaultdict(list) # Save chain and res number of binned residues
                bin2Resis = collections.defaultdict(list) # Save chain and res number of binned residues
                bin3Resis = collections.defaultdict(list) # Save chain and res number of binned residues
                bin4Resis = collections.defaultdict(list) # Save chain and res number of binned residues
                
                for res in bs.bs_res:
                    minD = res['min_dist']
                    if minD <= bin1:
                        bsResiCnts[0] += 1
                        if res['aa'] in res2group.keys():
                            bin1stats[res2group[res['aa']]] += 1 # record res type if one of 20 standard amino acids
                            bin1stats[res['aa']] += 1 # Record the amino acid identity
                        else:
                            bin1stats['ALT'] += 1 # Count non-standard amino acids as "ALT" for alternate
                            
                        bin1stats[res['struct']] += 1 # Record secondary structure
                        bin1Resis[res['reschain']].append(res['resnr'])
                        
                    elif minD <= bin2:
                        bsResiCnts[1] += 1
                        if res['aa'] in res2group.keys():
                            bin2stats[res2group[res['aa']]] += 1 # record res type if one of 20 standard amino acids
                            bin2stats[res['aa']] += 1 # Record the amino acid identity
                        else:
                            bin2stats['ALT'] += 1 # Count non-standard amino acids as "ALT" for alternate
                            
                        bin2stats[res['struct']] += 1 # Record secondary structure
                        bin2Resis[res['reschain']].append(res['resnr'])
                        
                    elif minD <= bin3:
                        bsResiCnts[2] += 1
                        if res['aa'] in res2group.keys():
                            bin3stats[res2group[res['aa']]] += 1 # record res type if one of 20 standard amino acids
                            bin3stats[res['aa']] += 1 # Record the amino acid identity
                        else:
                            bin3stats['ALT'] += 1 # Count non-standard amino acids as "ALT" for alternate
                            
                        bin3stats[res['struct']] += 1 # Record secondary structure
                        bin3Resis[res['reschain']].append(res['resnr'])
                        
                    elif minD <= bin4: # bin 4
                        bsResiCnts[3] += 1
                        if res['aa'] in res2group.keys():
                            bin4stats[res2group[res['aa']]] += 1 # record res type if one of 20 standard amino acids
                            bin4stats[res['aa']] += 1 # Record the amino acid identity
                        else:
                            bin4stats['ALT'] += 1 # Count non-standard amino acids as "ALT" for alternate
                            
                        bin4stats[res['struct']] += 1 # Record secondary structure
                        bin4Resis[res['reschain']].append(res['resnr'])
        
                # # Print bsite resiudes in proper format for pymol
                # # bin 1
                # print(lig)
                # out = 'sele '
                # for chainKey in bin1Resis.keys():
                #     if out == 'sele ':
                #         out += '(resi '
                #     else:
                #         out += ' or (resi '
                #     for r in bin1Resis[chainKey]:
                #         if out[-1] != ' ':
                #             out += '+'
                #         out += str(r)
                #     out += ' and chain ' + chainKey + ')'
                # print(out)
                # print('set_name sele, bin1')
                
                # # bin 2
                # out = 'sele '
                # for chainKey in bin2Resis.keys():
                #     if out == 'sele ':
                #         out += '(resi '
                #     else:
                #         out += ' or (resi '
                #     for r in bin2Resis[chainKey]:
                #         if out[-1] != ' ':
                #             out += '+'
                #         out += str(r)
                #     out += ' and chain ' + chainKey + ')'
                # print(out)
                # print('set_name sele, bin2')
                
                # # bin 3
                # out = 'sele '
                # for chainKey in bin3Resis.keys():
                #     if out == 'sele ':
                #         out += '(resi '
                #     else:
                #         out += ' or (resi '
                #     for r in bin3Resis[chainKey]:
                #         if out[-1] != ' ':
                #             out += '+'
                #         out += str(r)
                #     out += ' and chain ' + chainKey + ')'
                # print(out)
                # print('set_name sele, bin3')
                
                # # bin 4
                # out = 'sele '
                # for chainKey in bin4Resis.keys():
                #     if out == 'sele ':
                #         out += '(resi '
                #     else:
                #         out += ' or (resi '
                #     for r in bin4Resis[chainKey]:
                #         if out[-1] != ' ':
                #             out += '+'
                #         out += str(r)
                #     out += ' and chain ' + chainKey + ')'
                # print(out)
                # print('set_name sele, bin4')
                    
                
                bin1out = [ str(bin1stats[c]/bsResiCnts[0]) if bsResiCnts[0] != 0 else '0' for c in dsspCodes + aaGroups + resis] # Percentages of amino acids with each ss and res characteristic in the order provided in the column names
                bin2out = [ str(bin2stats[c]/bsResiCnts[1]) if bsResiCnts[1] != 0 else '0' for c in dsspCodes + aaGroups + resis] # if/else to avoid divide by 0 error
                bin3out = [ str(bin3stats[c]/bsResiCnts[2]) if bsResiCnts[2] != 0 else '0' for c in dsspCodes + aaGroups + resis]
                bin4out = [ str(bin4stats[c]/bsResiCnts[3]) if bsResiCnts[3] != 0 else '0' for c in dsspCodes + aaGroups + resis]
                bsResiCnts = [str(c) for c in bsResiCnts]
                
                hetAtmOut = [ str(hetAtmTracker[pdb][lig][hetAtmCnt]) for hetAtmCnt in hetAtmResiCnts]
                        
                outFH.write('\t'.join([pdb + '_' + lig] + [pdb] + LigCheck.unilecData[pdb] + [str(len(plipDat.bsites))] + [bs.longname] + [bsSeqs] + bsResiCnts + [str(bs.counts[intType]) for intType in intNames] + bin1out + bin2out + bin3out + bin4out + hetAtmOut) + '\n')
                
    

# Investigating residue at position 155 in different HA types

# HA_pdbs = {}
# HA_pdbs['H1'] = {'a23' : ["1RV0","1RVX","3HTT","3UBJ","3UBQ"],
#                  'a26' : ["1RVT","1RVZ","3HTQ","3UBE","3UBN"]}
# HA_pdbs['H3'] = {'a23' : ["1HGG","1MQM"],
#                  'a26' : ["1MQN"]}
# HA_pdbs['H5'] = {'a23' : ["3ZNK","3ZNL","3ZNM","3ZPB","4BGX","4BGY","4BH1","4CQQ","4CQW","4CQY","4CQZ"],
#                  'a26' : ["4BH0","4BH3","4BH4","4CQR","4CQU","4CQX"]}
# HA_pdbs['H7'] = {'a23' : ["3M5H","4BSI"],
#                  'a26' : ["3M5I","4BSB","4BSC","4BSD","4BSE","4BSF","4BSH"]}
# HA_pdbs['H10'] = {'a26' : ["4D00"]}
# HA_pdbs['B'] = {'a23' : ["2RFT"],
#                 'a26' : ["2RFU"]}


# resi_positions = {'H1' : [155], # Analogous positions for T155V in H1 for other types/subtypes
#                    'H3' : [155],
#                    'H5' : [150, 151],
#                    'H7' : [144, 155],
#                    'H10': [146],
#                    'B'  : [160]}

# HA_resis = collections.defaultdict(lambda: collections.defaultdict(list)) # Nest default dicts to hold {['HA type'] : ['ligand type'] : [residue list]} same structure as HA_pdbs


# for h in HA_pdbs.keys():
#     print(h)
#     for l in HA_pdbs[h].keys():
#         print(l)
#         for pdb in HA_pdbs[h][l]:
#             print(pdb)
#             plipDat = allPLIP[pdb]
#             for lig in plipDat.bsites:
#                 print(lig)
#                 bs = plipDat.bsites[lig]
#                 out = []
#                 for res in bs.bs_res:
#                     if res['resnr'] in resi_positions[h]:
#                         out.append(res)
#                 if out:
#                     if len(out) == 2:
#                         out = out[np.argmin([o['min_dist'] for o in out])]
#                     HA_resis[h][l].extend(out)
#                 #     print(out)
#                 # else:
#                 #     print('----------------\n\tMISSING ' + pdb + ' ' + lig + '\n----------------\n')
                
# ## H1
# h = 'H1'

# h = 'B'

# print("6' with Val")
# print(','.join([str(r['min_dist']) for r in HA_resis[h]['a26'] if r['aa'] == 'VAL']))
# print('Contacts at distances: ' + ','.join([str(r['min_dist']) for r in HA_resis[h]['a26'] if r['aa'] == 'VAL' and r['contact']]))

# print("3' with Val")
# print(','.join([str(r['min_dist']) for r in HA_resis[h]['a23'] if r['aa'] == 'VAL']))
# print('Contacts at distances: ' + ','.join([str(r['min_dist']) for r in HA_resis[h]['a23'] if r['aa'] == 'VAL' and r['contact']]))
# print()


# print("6' with Thr")
# print(','.join([str(r['min_dist']) for r in HA_resis[h]['a26'] if r['aa'] == 'THR']))
# print('Contacts at distances: ' + ','.join([str(r['min_dist']) for r in HA_resis[h]['a26'] if r['aa'] == 'THR' and r['contact']]))

# print("3' with Thr")
# print(','.join([str(r['min_dist']) for r in HA_resis[h]['a23'] if r['aa'] == 'THR']))
# print('Contacts at distances: ' + ','.join([str(r['min_dist']) for r in HA_resis[h]['a23'] if r['aa'] == 'THR' and r['contact']]))





    