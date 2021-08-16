#!/bin/bash

subcnt=1
while [ $subcnt -le 10 ]; do
mkdir ./$subcnt
cd ./$subcnt
cp /dartfs-hpc/rc/home/y/f002tsy/cbklab/Mattox/lec_gly_binding/scripts/prediction/cluster/${1} ./
echo "Submitting ${subcnt}"
sbatch ${1}
cd ..
let subcnt=subcnt+1
done
