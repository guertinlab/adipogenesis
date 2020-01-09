#!/bin/bash
#SBATCH -n 1
#SBATCH -t 96:00:00
#SBATCH -o minC100_191202_pro.out
#SBATCH -p largemem
#SBATCH -A guertinlab

cd ~/pro_seq
module load gcc R
Rscript LRT_largemem_clustering_pro.R
