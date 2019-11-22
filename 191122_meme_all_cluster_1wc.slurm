#!/bin/bash 
#SBATCH -p largemem
#SBATCH -N 1
#SBATCH --output=191122_meme_cluster_1wc.out
#SBATCH -A guertinlab

module load gcc/7.1.0 meme/5.1.0
module load openmpi
module load bedtools

./191122_meme_cluster_1wc.sh

