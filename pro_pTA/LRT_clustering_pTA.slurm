#!/bin/bash
#SBATCH -n 1
#SBATCH -t 96:00:00
#SBATCH -o minC100_191230_pTA_1e5.out
#SBATCH -p largemem
#SBATCH -A guertinlab

cd ~/pro_pTA
module load gcc R
Rscript LRT_largemem_clustering_pTA.R
