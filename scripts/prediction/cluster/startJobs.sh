#!/bin/bash

cnt=1
while [ $cnt -le 15 ]; do
mkdir ./$cnt
cp /Users/dmattox/cbk/lec_gly_binding/scripts/prediction/cluster/runTrainValidate.pbs ./$cnt
cd ./$cnt
echo "Submitting ${cnt}"
mksub runTrainValidate.pbs
cd ..
let cnt=cnt+1
done

