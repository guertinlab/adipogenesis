#!/bin/bash
#SBATCH -n 1
#SBATCH -t 96:00:00
#SBATCH -p standard
#SBATCH -A janeslab
#SBATCH --array=1-8

sh ChIP.pt1.sh $SLURM_ARRAY_TASK_ID
