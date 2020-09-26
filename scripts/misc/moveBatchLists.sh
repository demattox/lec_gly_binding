#!/bin/bash

cd /dartfs-hpc/rc/home/y/f002tsy/cbklab/Mattox/lec_gly_binding/data/structures/bsites/batchLists
mkdir ../batchProcessing
for fn in *; do [ -f $fn ] || continue
        mkdir ../batchProcessing/${fn%.txt}
        cp $fn ../batchProcessing/${fn%.txt}/
done
