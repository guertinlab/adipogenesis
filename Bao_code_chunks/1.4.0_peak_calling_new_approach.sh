#!/bin/bash

#archive all files for 6d treatment; not used for analysis
mkdir 6day
mv *6d* 6day/

###merge .bam, call peaks, remove blacklisted regions, set summit windows

#generate a unique slurm file for each replicate and run them in parallel:
#header_1   --> sbatch settings
#temp.txt   --> name of .out file
#header_2   --> modules to load
#temp2.txt  --> name of relevant file
#header_3   --> actual commands
for bam in *rep1_atac_rmdup.bam
do
    name=$(echo $bam | awk -F"_rep1_atac_rmdup.bam" '{print $1}')
    echo $name
    echo '#SBATCH -o' $name'.peak.calling.out' > temp.txt
    echo 'name='$name > temp2.txt
    cat peak_calling_slurm_header_1.txt temp.txt peak_calling_slurm_header_2.txt temp2.txt peak_calling_slurm_header_3.txt > $name.peak.calling.slurm
    sbatch $name.peak.calling.slurm
done
