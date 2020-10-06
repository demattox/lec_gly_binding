#!/bin/bash


# Created on Mon Sep 28 15:09:01 2020
# @author: dmattox

cd /Users/dmattox/cbk/lec_gly_binding/data/structures/holo/seqs/nrSeqs/

outDir="/Users/dmattox/cbk/lec_gly_binding/analysis/seqClustering/"

/Users/dmattox/cdhit/cd-hit -i allSeqs.fst -o ${outDir}id50/clusts -c 0.5 -G 1 -g 1 -b 20 -l 10 -n 3

/Users/dmattox/cdhit/cd-hit -i allSeqs.fst -o ${outDir}id90 -c 0.5 -G 1 -g 1 -b 20 -l 10