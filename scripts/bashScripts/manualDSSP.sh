#!/bin/bash

cd /dartfs-hpc/rc/home/y/f002tsy/cbklab/Mattox/lec_gly_binding/data/dssp/origPDB
wget https://files.rcsb.org/download/1BOS.cif
wget https://files.rcsb.org/download/6TIG.cif
wget https://files.rcsb.org/download/6TID.cif
cd ../
/dartfs-hpc/rc/home/y/f002tsy/cbklab/software/dssp-2.3.0/mkdssp -i ./origPDB/6TIG.cif -o ./dsspOut/6TIG.dssp
/dartfs-hpc/rc/home/y/f002tsy/cbklab/software/dssp-2.3.0/mkdssp -i ./origPDB/6TID.cif -o ./dsspOut/6TID.dssp
/dartfs-hpc/rc/home/y/f002tsy/cbklab/software/dssp-2.3.0/mkdssp -i ./origPDB/1BOS.cif -o ./dsspOut/1BOS.dssp

exit 0
