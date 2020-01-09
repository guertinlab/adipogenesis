#!/bin/bash 
#SBATCH --output=191201_meme_pro_clusters.out
#SBATCH -n 64
#SBATCH -t 72:00:00
#SBATCH -p parallel
#SBATCH -A guertinlab

module load gcc/7.1.0 openmpi/3.1.4 meme/5.1.0
module load bedtools

cd ~/pro_seq
./191201_meme_pro_clusters.sh

