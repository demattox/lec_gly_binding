#!/bin/bash

outF="/Users/dmattox/cbk/lec_gly_binding/analysis/d2_distributionFeats.csv"
touch $outF
cd /Users/dmattox/cbk/lec_gly_binding/data/structures/bsites/batchProcessing
fn="d2DistributionFeatures.csv"
head -n 1 ./pdbList_0/${fn} > $outF
for batchDir in *; do [ -d $batchDir ] || continue
    sed '1'd "${batchDir}/${fn}" >> $outF
done