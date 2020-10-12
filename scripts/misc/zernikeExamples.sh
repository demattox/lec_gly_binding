#!/usr/bin/bash

"""
Created on Mon Oct 12 13:16:14 2020

@author: dmattox
"""

outDir='/dartfs-hpc/rc/home/y/f002tsy/cbklab/Mattox/lec_gly_binding/analysis/zernikeExamples'

arr=(4Z4Z_GLA-D-1 4LK7_04G-D-202 2JDH_TA5-B-1117 4FMH_GAL-X-1 1BOS_GAL-R-4370 4X07_GLA-D-1 4E52_GMH-E-1)

cnt=0
while [ $cnt -lt ${#arr[@]} ]; do
    pdb=${arr[$cnt]}
    echo $pdb
    cp /dartfs-hpc/rc/home/y/f002tsy/cbklab/Mattox/lec_gly_binding/data/structures/bsites/bsiteVoxels/${pdb}_4AngBin.mol2 ${outDir}
    cp /dartfs-hpc/rc/home/y/f002tsy/cbklab/Mattox/lec_gly_binding/data/structures/bsites/zernVoxels/${pdb}_4.vox ${outDir}
    
    /dartfs-hpc/rc/home/y/f002tsy/cbklab/software/zernikeInvariants/bin/vox ${outDir}/${pdb}_4.vox 10 ${outDir}/${pdb}_10th 1
    /dartfs-hpc/rc/home/y/f002tsy/cbklab/software/zernikeInvariants/bin/vox ${outDir}/${pdb}_4.vox 20 ${outDir}/${pdb}_20th 1
    let cnt=cnt+1
done
