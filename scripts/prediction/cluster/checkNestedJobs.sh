#!/bin/bash

for d in *; do
for sd in ./$d/*; do 
if [ ! -f ./$sd/train_val.e* ]; then             
echo $sd
cd $sd
mksub *.pbs     
cd ../../
fi
done
done
