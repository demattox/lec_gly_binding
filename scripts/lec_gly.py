#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 28 15:16:56 2020

@author: dmattox
"""

#homeDir = '/dartfs-hpc/rc/home/y/f002tsy/cbklab/Mattox/lec_gly_binding/'
homeDir = '/Users/dmattox/cbk/lec_gly_binding'

if homeDir[-1] != '/':
    homeDir += '/'

# import os
import numpy as np
# import sys
from rdkit import Chem
# import pliptool.plip.modules.plipxml as plipxml
# import Bio.PDB


# os.chdir(homeDir)


#################
# Function toolbox
#################

## Amino acid code conversion dictionaries
three2one = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M', 'MSE':'M', 'UNK':'U'}
one2three = {v:k for k,v in three2one.items()} # Flip keys and values of previous dict

def res2plipID(bioPDBresi):
    ''' Generate a unique id from Bio.PDB residue following the PLIP format, automated residue vs hetatm detection '''
    if bioPDBresi.full_id[3][0][:2].strip() == 'H_':  # If hetatm
        out = bioPDBresi.full_id[3][0][2:].strip() + ':' + bioPDBresi.full_id[2] + ':' + str(bioPDBresi.full_id[3][1])
    else:
        out = bioPDBresi.get_resname() + ':' + bioPDBresi.full_id[2] + ':' + str(bioPDBresi.full_id[3][1])
    return out


def resCentroid(bioPDBresi):
    ''' Calculate the centroid of a Bio.PDB residue'''
    out = [0,0,0]
    atmCnt = 0
    for atm in bioPDBresi.get_atoms():
        atmCnt += 1
        out[0] += atm.coord[0]; out[1] += atm.coord[1]; out[2] += atm.coord[2]
    out  = [c/atmCnt for c in out]
    return out


def eucDist(coord1, coord2):
    ''' Calculate the euclidean distance between a pair of 3D coordinates in separate lists '''
    return np.sqrt((coord1[0]-coord2[0])**2 + (coord1[1] - coord2[1])**2 + (coord1[2] - coord2[2])**2)


def findClosestAtmPair(link, biopdbModel, thresh, warn = True):
    ''' For a pair of PLIP IDs for ligands in a link tuple, returns the pair of atoms closest to the other residue. If warn is set to True,
    this functon returns a warning if more than one pair is closer than the specified threshold.'''

    res1 = link[0].split(':') # Split PLIP res ID into [het id, chain, res num]
    res2 = link[1].split(':')
    l1 = biopdbModel[res1[1]][('H_'+res1[0], int(res1[2]), ' ')]
    l2 = biopdbModel[res2[1]][('H_'+res2[0], int(res2[2]), ' ')]

    atmPairs = {} # Store the distances between each pair of atoms from the two
    for a1 in l1.get_atoms():
        for a2 in l2.get_atoms():
            atmPairs[eucDist(a1.coord, a2.coord)] = (a1.get_serial_number(), a2.get_serial_number())

    if len([p for p in atmPairs.keys() if p <= thresh]) > 1:
        print("Warning!!!\n\tMultiple pairs of atoms detected within covalent bond distance in residue pair:\n\t" + str(link) + " in PDB " + biopdbModel.get_parent().id)
        print("\tReturning the closest pair of atoms with PDB atom numbers:", atmPairs[min(atmPairs.keys())])
        
    return atmPairs[min(atmPairs.keys())]


def findClosestResPair(resi, resiLst):
    '''  For a bioPDB residue (resi) and a list of bioPDB residues (resiLst), returns a single resiLst residue with the smallest distance between any of its atoms and 
    the atoms in the lone residue (resi), as well as the distance '''
    minDist = float('inf')
    bestRes = ''
    for r2 in resiLst:
        for a1 in resi.get_atoms():
            for a2 in r2.get_atoms():
                d = a1 - a2
                if d < minDist:
                    minDist = d
                    bestRes = r2
    return bestRes, minDist

def newMinDist(proResi, glyResLst):
    '''  For a bioPDB protein residue (proResi) and a list of bioPDB glycan residues (glyResLst), returns the smallest distance between any of the protein atoms and 
    any HEAVY atoms in the glycan residues (any protein atom to any heavy ligand to match PLIP's min_dist values)'''
                                            
    minDist = float('inf')
    for glyRes in glyResLst:
        for proAtm in proResi.get_atoms():
            for glyAtm in glyRes.get_atoms():
                if glyAtm.element != 'H': # Skip ligand hydrogens
                    d = proAtm - glyAtm
                    if d < minDist:
                        minDist = d
    return minDist

def glyInMol(smilesString):
    ''' Returns boolean if at least one ring in the molecule described by the provided SMILES string meets the criteria of:
        Being non-aromatic, containing 4-9 atoms, with only one oxygen and the remaining atoms are carbons'''

    COset = set(['C','O'])
    out = False # Flag vriable to record if a ring meets the above criteria for a glycan

    m = Chem.MolFromSmiles(smilesString, sanitize = False) # Create rdkit.Chem molecule objects from ligand SMILES, sanitize False to avoid valency errors
    Chem.FastFindRings(m) # Do a non-SSSR ring find
    ri = m.GetRingInfo()

    nRings = ri.NumRings() # Count rings in molecule

    if nRings == 0:
        return out # No rings detected, not likely a glycan
    else:
        for ring in ri.BondRings():

            bondtypes = set([str(m.GetBondWithIdx(idx).GetBondType()) for idx in ring]) # Check the type of bonds the compsing the ring
            atmidxs = set() # Holds atom indices for atoms in each bond
            for idx in ring: # For each bond, get the bond index & the atom indices for the bond
                b = m.GetBondWithIdx(idx)
                atmidxs.add(b.GetBeginAtomIdx())
                atmidxs.add(b.GetEndAtomIdx())
            atms = [m.GetAtomWithIdx(ai).GetSymbol() for ai in atmidxs] # convert set of atom indices to atom symbols

            if ( 4 <= len(atms) <= 9 and # Check atom count, bond types, only carbon/oxygen, & only one oxygen
                'AROMATIC' not in bondtypes and
                set(atms) == COset and
                atms.count('O') == 1 ):
                    out = True
    return out

def glyInMolID(hetIDs, id2SMILES):
    ''' Returns boolean if at least one ring in the molecule described by the provided string of heteroatom IDs meets the criteria of:
        Being non-aromatic, containing 4-9 atoms, with only one oxygen and the remaining atoms are carbons'''

    COset = set(['C','O'])
    out = False # Flag variable to record if a ring meets the above criteria for a glycan
    
    ids = hetIDs.split('-')
    
    for i in ids:
        smi = id2SMILES[i]
    
        m = Chem.MolFromSmiles(smi, sanitize = False) # Create rdkit.Chem molecule objects from ligand SMILES, sanitize False to avoid valency errors
        Chem.FastFindRings(m) # Do a non-SSSR ring find
        ri = m.GetRingInfo()
    
        nRings = ri.NumRings() # Count rings in molecule
    
        if nRings == 0:
            continue # No rings detected, not likely a glycan
        else:
            for ring in ri.BondRings():
    
                bondtypes = set([str(m.GetBondWithIdx(idx).GetBondType()) for idx in ring]) # Check the type of bonds the compsing the ring
                atmidxs = set() # Holds atom indices for atoms in each bond
                for idx in ring: # For each bond, get the bond index & the atom indices for the bond
                    b = m.GetBondWithIdx(idx)
                    atmidxs.add(b.GetBeginAtomIdx())
                    atmidxs.add(b.GetEndAtomIdx())
                atms = [m.GetAtomWithIdx(ai).GetSymbol() for ai in atmidxs] # convert set of atom indices to atom symbols
    
                if ( 4 <= len(atms) <= 9 and # Check atom count, bond types, only carbon/oxygen, & only one oxygen
                    'AROMATIC' not in bondtypes and
                    set(atms) == COset and
                    atms.count('O') == 1 ):
                        out = True
    return out

def reformatDSSP(dsspDict):
    '''  Takes Bio.PDB.DSSP dist dictionary and changes keys to string of resnum and chain (ie '25A') '''
    out = {}
    for k in dsspDict.keys():
        newkey = str(k[1][1]) +k[0]
        out[newkey] = dsspDict[k]
    return out

def manualDSSPpull(dsspFH, residueID, insertionCode = ' '):
    '''
    Adaptation of make_dssp_dict from Bio.PDB.DSSP https://biopython.org/DIST/docs/api/Bio.PDB.DSSP%27-pysrc.html#_make_dssp_dict
    Extract DSSP record for a single residue if it was overwritten by original make_dssp_dict function
    
    Input: filehandle for dssp file, residue ID to get info for with chain ID (ie resi num 6 on chain A: '6A'
    Returns: [struct, acc, phi, psi]
    '''
    start = 0 
    out = []
    resChain = residueID[-1]
    resNum = int(residueID[:-1])
    with open(dsspFH, 'r') as handle:
        for l in handle.readlines(): 
            sl = l.split() 
            if len(sl) < 2: 
                continue 
            if sl[1] == "RESIDUE": 
                # Start parsing from here 
                start = 1 
                continue 
            if not start: 
                continue 
            if l[9] == " ": 
                # Skip -- missing residue 
                continue
            resseq = int(l[5:10]) 
            icode = l[10] 
            chainid = l[11]
            #aa = l[13] 
            ss = l[16]
            if ss == " ": 
                ss = "-" 
            try: 
                acc = int(l[34:38]) 
                phi = float(l[103:109]) 
                psi = float(l[109:115])
            except ValueError as exc:
                # DSSP output breaks its own format when there are >9999 
                # residues, since only 4 digits are allocated to the seq num 
                # field.  See 3kic chain T res 321, 1vsy chain T res 6077. 
                # Here, look for whitespace to figure out the number of extra 
                # digits, and shift parsing the rest of the line by that amount. 
                if l[34] != " ": 
                    shift = l[34:].find(" ") 

                    acc = int((l[34 + shift : 38 + shift])) 
                    phi = float(l[103 + shift : 109 + shift]) 
                    psi = float(l[109 + shift : 115 + shift]) 
                else: 
                    raise ValueError(exc)
            if chainid == resChain and resseq == resNum and icode == insertionCode:
                out = [ss, acc, phi, psi]
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

def getCentroid(atmLst):
    '''Given a list of BioPDB atoms, returns the centroid of their coordinates as an array'''
    out = np.array([0,0,0], dtype = 'float32') # Holds the centroid
    for a in atmLst:
        out += a.coord
    out = out/len(atmLst)
    return out

def getArrayCentroid(atmArr):
    '''Given an array of coordiantes, returns the centroid of their coordinates as an array'''
    out = np.array([0,0,0], dtype = 'float32') # Holds the centroid
    for a in atmArr:
        out += a
    out = out/len(atmArr)
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

def cntLocalMaxima(vals):
    ''' From a 1D array or list of integer/float values, returns the number of local maxima,
    where a local maximum is defined as a point i that is larger than points i-1 & i+1
    '''
    out = 0
    for i in range(1,len(vals)-1):
        if vals[i] > vals[i-1] and vals[i] > vals[i+1]:
            out += 1
    return out
    

def smooth(x,window_len=11,window='hanning'):
    """smooth the data using a window with requested size.
    
    Source: adapted from https://scipy-cookbook.readthedocs.io/items/SignalSmooth.html
    
    original code had following copyright disclaimer:
    __________________________________________________________________________
    Copyright (c) 2001, 2002 Enthought, Inc.
    All rights reserved.
    
    Copyright (c) 2003-2017 SciPy Developers.
    All rights reserved.
    
    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:
    
      a. Redistributions of source code must retain the above copyright notice,
         this list of conditions and the following disclaimer.
      b. Redistributions in binary form must reproduce the above copyright
         notice, this list of conditions and the following disclaimer in the
         documentation and/or other materials provided with the distribution.
      c. Neither the name of Enthought nor the names of the SciPy Developers
         may be used to endorse or promote products derived from this software
         without specific prior written permission.
    
    
    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
    AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
    IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
    ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS
    BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
    OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
    SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
    INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
    CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
    ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
    THE POSSIBILITY OF SUCH DAMAGE.
    
    __________________________________________________________________________
    
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal
        
    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    
    see also: 
    
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
 
    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """

    if x.ndim != 1:
        raise ValueError("smooth only accepts 1 dimension arrays.")

    if x.size < window_len:
        raise ValueError("Input vector needs to be bigger than window size.")


    if window_len<3:
        return x


    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError("Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")


    s=np.r_[x[window_len-1:0:-1],x,x[-2:-window_len-1:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')

    y=np.convolve(w/w.sum(),s,mode='valid')
    trim = window_len // 2
    return y[trim:-trim]


