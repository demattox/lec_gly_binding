#!/bin/bash

cd /dartfs-hpc/rc/home/y/f002tsy/cbklab/Mattox/lec_gly_binding/data/plip/out/

for pdb in *; do 
    [ -d $pdb ] || continue
    cp ./${pdb}/report.xml ../PLIPreports/${pdb}.xml
    if [ -f ./$pdb/plipfixed.* ]; then
        cp ./${pdb}/plipfixed.* ../../structures/holo/plipFixed/${pdb}.pdb
    else
        cp ./${pdb}/*.pdb ../../structures/holo/plipFixed/${pdb}.pdb
    fi
done
cd ../
