#!/bin/bash

cnt=1
while [ $cnt -le 15 ]; do
mkdir ./$cnt
cp /dartfs-hpc/rc/home/y/f002tsy/cbklab/Mattox/lec_gly_binding/scripts/prediction/cluster/${1} ./$cnt
cd ./$cnt
echo "Submitting ${cnt}"
mksub ${1}
cd ..
let cnt=cnt+1
done

