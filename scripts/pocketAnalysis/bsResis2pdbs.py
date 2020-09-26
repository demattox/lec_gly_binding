#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 27 17:42:44 2020

@author: dmattox

Save out all binding site residues indicated in the cleaned PLIP records to a separate PDB file
"""

import os
import glob
import dill
import Bio.PDB

import lec_gly as LecGly

os.chdir(LecGly.homeDir)

############################################

plipFile = './data/unilectin/cleanPLIPrecords.p'
zeroResMapping = './data/unilectin/icodeResi2PLIPmapping.p'

bsPDBDir = './data/structures/bsites/residuePDBs/'

############################################

if __name__ == '__main__':
    
    import checkPLIPligands as LigCheck
    
    pdbFiles = glob.glob(LigCheck.originalPDBsDir + '*.pdb')
    
    with open("./data/unilectin/cleanPLIPrecords.p", "rb") as pickleFH:
        allPLIP = dill.load(pickleFH)
    
    with open("./data/unilectin/icodeResi2PLIPmapping.p", "rb")as pickleFH:
        zeroResis = dill.load(pickleFH)
    
    ############################################
    
    # pdb = '1HGG'

    for k,pdb in enumerate(allPLIP.keys()):
        if k % 100 == 0: print("Processing PBD file number " + str(k) + ' of ' + str(len(allPLIP)))
        
        
        if pdb in LigCheck.cifLst:
            parser = Bio.PDB.MMCIFParser(QUIET=True)
            pdbIN = [f for f in LigCheck.cifFiles if f[-8:-4] == pdb][0]
        else:
            parser = Bio.PDB.PDBParser(QUIET=True)
            pdbIN = [f for f in pdbFiles if f[-8:-4] == pdb][0]

        struct = parser.get_structure(pdb, pdbIN)
        model = struct[0] # grab the first model in PDB file
        
        for bs in allPLIP[pdb].bsites:
            # Identify heavy binding site atoms
            
            bsResIDs = [res['aa'] + ':' + res['reschain'] + ':' + str(res['resnr']) for res in allPLIP[pdb].bsites[bs].bs_res]
            if pdb in zeroResis.keys(): # Check for zero numbered resiudes and correct the res num to match the original number including the i_code
                if bs in zeroResis[pdb]:
                    for i,resID in enumerate(bsResIDs):
                        resID = resID.split(':')
                        rID = resID[2]+resID[1] # resnum from PLIP with chain ID
                        if rID in zeroResis[pdb][bs].keys():
                            resID[2] = zeroResis[pdb][bs][rID] # Replace PLIP resnum with true resnum with icode
                            bsResIDs[i] = ':'.join(resID) 
            proAtms = [] # Holds all heavy atoms in residues within the binding site
            for res in model.get_residues():
                resID = LecGly.res2plipID(res)
                if res.full_id[3][2] != ' ': # If the residue has an i_code, append it to the resID to match the modified PLIP resIDs
                    resID = resID + res.full_id[3][2]
                if resID in bsResIDs:
                    for a in res.get_atoms():
                        if a.element != 'H':
                            proAtms.append(a)
            
            outFH = bsPDBDir + pdb + '_' + bs.replace(':','-') + '_' + str(len(proAtms)) + 'atms.pdb'
            
            class BsSelect(Bio.PDB.Select):
                def accept_atom(self, atom):
                    if atom in proAtms:
                        return 1
                    else:
                        return 0
            
            io = Bio.PDB.PDBIO()
            io.set_structure(struct)
            io.save(outFH, select=BsSelect())
