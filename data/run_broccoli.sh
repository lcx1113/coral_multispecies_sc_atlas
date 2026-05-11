#!/bin/bash
#BSUB -q normal
#BSUB -J broccoli
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -R "span[ptile=32]"
#BSUB -n 32

source ~/.bashrc
conda activate broccoli
cd /share/home/Hu_lab_JX/workplace/broccoli/Metazoa_genome_fix_ahya_plut/
python3 /share/home/Hu_lab_JX/tools/Broccoli-fix_step4/broccoli.py -threads 32 -dir ./ -nb_hits 10 -phylogenies ml -kmer_size 10000
