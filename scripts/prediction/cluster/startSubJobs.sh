#!/bin/bash

cnt=1
while [ $cnt -le 15 ]; do
mkdir ./$cnt
cd ./$cnt

subcnt=1
while [ $subcnt -le 10 ]; do
mkdir ./$subcnt
cd ./$subcnt
cp /dartfs-hpc/rc/home/y/f002tsy/cbklab/Mattox/lec_gly_binding/scripts/prediction/cluster/${1} ./
echo "Submitting ${cnt}:${subcnt}"
mksub ${1}
cd ..
let subcnt=subcnt+1
done

cd ..
let cnt=cnt+1
done

