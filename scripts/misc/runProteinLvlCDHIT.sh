#!/bin/bash


# Created on Mon Sep 28 15:09:01 2020
# @author: dmattox

cd /Users/dmattox/cbk/lec_gly_binding/data/structures/holo/seqs/nrSeqs/

outDir="/Users/dmattox/cbk/lec_gly_binding/analysis/seqClustering/"

/Users/dmattox/cdhit/cd-hit -i allSeqs.fst -o ${outDir}id50/out -c 0.5 -G 1 -g 1 -b 20 -l 10 -n 3
/Users/dmattox/cdhit/cd-hit -i allSeqs.fst -o ${outDir}id60/out -c 0.6 -G 1 -g 1 -b 20 -l 10 -n 4
/Users/dmattox/cdhit/cd-hit -i allSeqs.fst -o ${outDir}id70/out -c 0.7 -G 1 -g 1 -b 20 -l 10 -n 5
/Users/dmattox/cdhit/cd-hit -i allSeqs.fst -o ${outDir}id80/out -c 0.8 -G 1 -g 1 -b 20 -l 10 -n 5
/Users/dmattox/cdhit/cd-hit -i allSeqs.fst -o ${outDir}id90/out -c 0.9 -G 1 -g 1 -b 20 -l 10 -n 5
