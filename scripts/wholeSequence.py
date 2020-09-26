#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 22 10:31:33 2020

@author: dmattox

Read in a PDB file [arg 1] and write out the chain seqeunces to a FASTA file [arg2]
"""

import sys

""" Reading and constructing a PDB profile from its lines """
class PDBlines:
    def __init__(self, line):
        self.AtomName = line[12:16].strip()
        self.ResName = line[17:20].strip()
        self.Chain = line[20:22].strip()
        self.ResNum = int(line[22:26])
        #self.X = float(line[30:38])
        #self.Y = float(line[38:46])
        #self.Z = float(line[46:54])
        
    @staticmethod
    def read(File):
        pdblines = [PDBlines(l) for l in open(File).readlines() if l[:4] == "ATOM"]
        return pdblines

d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
     'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M',
     'MSE': 'M', 'UNK': 'X', 'XAA': 'X', 'ASX': 'B', 'GLX': 'Z', 'XLE': 'J'}

pdbFile = sys.argv[1]

if pdbFile[-4:] == '.pdb':
    pro=PDBlines.read(pdbFile)
else:
    sys.exit('ERROR: unrecognized file extenion, must be a .pdb file')

outF = sys.argv[2]

with open(outF, 'w') as outFH:
    chainID = ''
    chain = ''
    for l in pro:
        if l.AtomName == 'CA':
            r = d[l.ResName] if l.ResName in d.keys() else 'X' # Get one letter code if available, else use "X"
            if l.Chain != chainID:              # end of a chain
                if chain != '':                 # Save out old chain
                    outFH.write('>'+chainID+'\n')
                    outFH.write(chain+'\n')
                chain = r                       # Start storing new chain sequence
                chainID = l.Chain
            else:                               # Continue current chain
                chain += r

    outFH.write('>'+chainID+'\n')
    outFH.write(chain+'\n')


exit(0)
