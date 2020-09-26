#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import glob

import lec_gly as LecGly
#from pocketAnalysis.bsResis2pdbs import bsPDBDir

os.chdir(LecGly.homeDir)

######################################################

wrlDir = './data/structures/bsites/wrlSurfaces/'
if not os.path.exists(wrlDir):
    os.mkdir(wrlDir)

bsPDBDir = './data/structures/bsites/residuePDBs/'

# mol2Dir = '/dartfs-hpc/rc/home/y/f002tsy/cbklab/Mattox/glycans/unilec3d/structures/bSites/bsiteSurface/surfacePnts/'

######################################################

if __name__ == '__main__':
    from pymol import cmd#, stored

    print('Processing pdbs to wrl representations')
    # Save surface representations from pymol to .wrl files
    
    files = glob.glob(bsPDBDir + '*.pdb')
    
    for i,f in enumerate(files):
        if i % 100 == 0: print("Processing binding site PBD file number " + str(i) + ' of ' + str(len(files)))
        fOut = f.split('/')[-1]
        fOut = fOut.split('_')
        fOut = '_'.join(fOut[:-1])
        cmd.load(f)
        cmd.set('surface_quality', '0')
        cmd.show_as('surface', 'all')
        cmd.set_view('1,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,300,1')
        cmd.save(wrlDir + fOut+'.wrl')
        cmd.delete('all')
    
    print('DONE!\nVDW surfaces for each binding site residue pdb file saved to ' + LecGly.homeDir[:-1] + wrlDir[1:])
    # ######################################################
    
    # print('\n\nExtracting points from wrl files and saving to mol2 files...\n')
    # # Get surface points from .wrl files and save to .mol2 files
    # files = os.listdir(wrlDir)
    # for i,f in enumerate(files):
    #     if i % 100 == 0: print("Processing wrl file number " + str(i) + ' of ' + str(len(files)))
    #     pnts = []
    #     pnts = getPntsFromWRL(wrlDir + f)
    #     if pnts:
    #         coords2Mol2(mol2Dir + f[:-4] + '.mol2', pnts)
    #     else:
    #         print('\t',f,'has no points to extract!')
        
