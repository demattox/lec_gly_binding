conda activate py27

cd /dartfs-hpc/rc/home/y/f002tsy/cbklab/Mattox/lec_gly_binding/

python /dartfs-hpc/rc/home/y/f002tsy/cbklab/software/plip/plip/plipcmd.py -f ./data/structures/holo/origPDB/6TID.pdb -x --maxthreads 16 -o ./data/plip/out/6TID/
python /dartfs-hpc/rc/home/y/f002tsy/cbklab/software/plip/plip/plipcmd.py -f ./data/structures/holo/origPDB/6TIG.pdb -x --maxthreads 16 -o ./data/plip/out/6TIG/
python /dartfs-hpc/rc/home/y/f002tsy/cbklab/software/plip/plip/plipcmd.py -f ./data/structures/holo/origPDB/1BOS.pdb -x --maxthreads 16 -o ./data/plip/out/1BOS/
