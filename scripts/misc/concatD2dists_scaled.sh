#!/bin/bash

outF="/Users/dmattox/cbk/lec_gly_binding/analysis/d2_scaled_bins.csv"
touch $outF
cd /Users/dmattox/cbk/lec_gly_binding/data/structures/bsites/batchProcessing
fn="scaledBins_D2Distances.csv"
head -n 1 ./pdbList_0/${fn} > $outF
for batchDir in *; do [ -d $batchDir ] || continue
    sed '1'd "${batchDir}/${fn}" >> $outF
done