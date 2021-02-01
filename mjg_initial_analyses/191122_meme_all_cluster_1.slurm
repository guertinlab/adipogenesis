#!/bin/bash 
#SBATCH --output=191122_meme_cluster_1.out
#SBATCH -n 60
#SBATCH --mem-per-cpu=2000
#SBATCH -t 72:00:00
#SBATCH -p parallel
#SBATCH -A guertinlab

module purge
module load gmvapich2/7.1.0_2.3.1 meme/5.1.0
module load bedtools

cd ~/clusters
./191122_meme_cluster_1.sh

