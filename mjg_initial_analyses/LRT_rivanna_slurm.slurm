#!/bin/bash                                                                                                                                                                                                                       
#SBATCH -n 1                                                                                                                                                                                                                      
#SBATCH -t 96:00:00                                                                                                                                                                                                               
#SBATCH -o minC100_191119_1.out                                                                                                                                                                                                   
#SBATCH -p largemem                                                                                                                                                                                                               
#SBATCH -A guertinlab                                                                                                                                                                                                             

cd ~
module load gcc R
Rscript LRT_largemem_clustering.R
