#!/bin/bash

cd /dartfs-hpc/rc/home/y/f002tsy/cbklab/Mattox/lec_gly_binding/

outF='./analysis/3DZDs_20thOrder.tsv'
touch ${outF}

for fn in ./data/structures/bsites/zernPolys/*; do [ -f $fn ] || continue
    name=$(basename "${fn%.inv}")
    echo -e ${name}'\t\c' >> ${outF}
    sed '1'd ${fn} >> ${outF}
done