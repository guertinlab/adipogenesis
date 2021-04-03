#!/bin/bash
#SBATCH -n 1
#SBATCH -t 96:00:00
#SBATCH -o create.final.atac.all.df.out
#SBATCH -p standard
#SBATCH -A guertinlab

module load gcc/7.1.0  openmpi/3.1.4 R/4.0.0

cd /scratch/bhn9by/ATAC

Rscript create.final.atac.all.df.R
