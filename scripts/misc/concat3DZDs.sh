#!/bin/bash

outF='/Users/dmattox/cbk/lec_gly_binding/analysis/3DZDs_20thOrder.tsv'
touch ${outF}
cd /Users/dmattox/cbk/lec_gly_binding/data/structures/bsites/zernPolys
for fn in *; do [ -f $fn ] || continue
    echo ${fn%.inv}'\t\c' >> ${outF}
    sed '1'd ${fn} >> ${outF}
done