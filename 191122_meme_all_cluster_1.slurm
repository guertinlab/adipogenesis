#!/bin/bash 
#SBATCH --nodes=1
#SBATCH --time=96:00:00
#SBATCH --output=191122_meme_cluster_1.out
#SBATCH --partition=largemem
#SBATCH -A guertinlab

module load gcc/7.1.0 meme/5.1.0
module load openmpi
module load bedtools

cd ~/clusters
./191122_meme_cluster_1.sh

