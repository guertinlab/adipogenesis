#!/bin/bash
#SBATCH -p standard
#SBATCH -N 1
#SBATCH --ntasks-per-node 40
#SBATCH -t 96:00:00
#SBATCH -A janeslab 
#SBATCH --array=1-9

sh ChIP.pt3.sh $SLURM_ARRAY_TASK_ID
