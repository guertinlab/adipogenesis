#!/bin/bash 
#SBATCH --nodes=1
#SBATCH --time=96:00:00
#SBATCH --output=191122_meme_cluster_gradual_up.out
#SBATCH --partition=largemem
#SBATCH -A guertinlab

module load gcc/7.1.0 meme/5.1.0
module load openmpi
module load bedtools

cd ~/clusters/final_clusters
./191122_meme_cluster_gradual_up.sh

