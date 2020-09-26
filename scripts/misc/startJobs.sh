#!/bin/bash

# Provide a pbs script as the 1st (& only) positional argument, copies that script to each subdirectory and submits the job to the Discovery cluster
# Note: run this script from its original directory (./scripts/misc/)
# Example:
## sh startJobs.sh runVoxGen.pbs

scriptDir='/dartfs-hpc/rc/home/y/f002tsy/cbklab/Mattox/lec_gly_binding/scripts/misc/'
cd /dartfs-hpc/rc/home/y/f002tsy/cbklab/Mattox/lec_gly_binding/data/structures/bsites/batchProcessing

for batchDir in *; do [ -d $batchDir ] || continue
        cp ${scriptDir}${1} ./${batchDir}
        cd ./${batchDir}
        echo 'submitting job for '${batchDir}
        mksub ${scriptDir}${1}
        cd /dartfs-hpc/rc/home/y/f002tsy/cbklab/Mattox/lec_gly_binding/data/structures/bsites/batchProcessing
done
