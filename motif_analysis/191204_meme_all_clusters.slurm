#!/bin/bash                                                                                                                                                                                                                       
#SBATCH --output=191204_meme_cluster_1.out                                                                                                                                                                                        
#SBATCH -p standard                                                                                                                                                                                                               
#SBATCH -N 1                                                                                                                                                                                                                      
#SBATCH --ntasks-per-node 40                                                                                                                                                                                                      
#SBATCH -t 7-00:00:00                                                                                                                                                                                                             
#SBATCH -A guertinlab                                                                                                                                                                                                             

module purge
module load gmvapich2/7.1.0_2.3.1 meme/5.1.0
module load bedtools

cd ~/clusters
./191204_meme_cluster.sh
