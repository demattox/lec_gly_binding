#!/bin/bash

cnt=1
while [ $cnt -le 15 ]; do
mkdir ./$cnt
cp /dartfs-hpc/rc/home/y/f002tsy/cbklab/Mattox/lec_gly_binding/scripts/prediction/cluster/runTrainValidate.pbs ./$cnt
cd ./$cnt
echo "Submitting ${cnt}"
mksub runTrainValidate.pbs
cd ..
let cnt=cnt+1
done

