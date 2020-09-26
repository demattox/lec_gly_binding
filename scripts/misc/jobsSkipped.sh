#!/bin/bash

cd /dartfs-hpc/rc/home/y/f002tsy/cbklab/Mattox/lec_gly_binding/data/structures/bsites/batchProcessing
for batchDir in *; do [ -d $batchDir ] || continue
        if [ ! -f ./${batchDir}/${1} ]; then
                cd ./${batchDir}
                echo 'resubmitting job for '${batchDir}
                mksub $2
                cd ../
        fi
done
