#!/bin/bash
#BSUB -q normal
#BSUB -J broccoli
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -R "span[ptile=32]"
#BSUB -n 32

source ~/.bashrc
conda activate broccoli
cd /absolute/path/broccoli_work/Metazoa/
python3 /absolute/path/Broccoli-fix_step4/broccoli.py -threads 32 -dir ./ -nb_hits 10 -phylogenies ml -kmer_size 10000
