#!/bin/bash

outD="/Users/dmattox/cbk/lec_gly_binding/analysis/allD2binnedDists"
cd /Users/dmattox/cbk/lec_gly_binding/data/structures/bsites/batchProcessing
fn="allD2Distances.csv"
for batchDir in *; do [ -d $batchDir ] || continue
    cp ./$batchDir/$fn "${outD}/${batchDir}_binnedD2.csv"
done