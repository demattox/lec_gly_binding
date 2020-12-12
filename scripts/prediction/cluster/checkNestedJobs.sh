#!/bin/bash

for d in *; do
for sd in ./$d/*; do 
if [ ! -f ./$sd/train_val.e* ]; then             
echo $sd
cd $sd
mksub *.pbs     
cd /dartfs-hpc/rc/home/y/f002tsy/cbklab/Mattox/lec_gly_binding/analysis/training/train_and_validate/nr_CVLOCO_reps/pred/
fi
done
done
