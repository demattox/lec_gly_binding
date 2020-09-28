#!/bin/bash


# Created on Sat Sep 26 07:42:41 2020
# @author: dmattox

cd /Users/dmattox/cbk/lec_gly_binding/data/structures/holo/seqs/allSeqs

outDir="/Users/dmattox/cbk/lec_gly_binding/data/structures/holo/seqs/allSeqsCDHIT90/"
[ -d $outDir ] || mkdir $outDir

for pdb in *.fst; do [ -f $fn ] || continue
    /Users/dmattox/cdhit/cd-hit -i ${pdb} -o ${outDir}${pdb%.fst}_nr.fst -c 0.9 -n 5 -d 0 -G 1 -g 1 -b 20 -l 10 
done
