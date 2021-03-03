#!/bin/bash

#perform seqOutBias and convert .bam to .bigwig

#make a unique slurm file for each replicate and run them in parallel:
#header_1   --> sbatch settings
#temp.txt   --> name of .out file
#header_2   --> more sbatch settings and modules to load
#temp2.txt  --> name of relevant file
#header_3   --> actual commands

for bam in *_pro_plus.bam
do
name=$(echo $bam | awk -F"_pro_plus.bam" '{print $1}')
echo $name
    echo '#SBATCH -o' $name'.bigwig.out' > temp.txt
    echo 'name='$name > temp2.txt
    cat bigwig_slurm_header_1.txt temp.txt bigwig_slurm_header_2.txt temp2.txt bigwig_slurm_header_3.txt > $name.bigwig.slurm
    rm temp.txt
    rm temp2.txt
done

#caution: you can't use the tallymer mappability file from ATAC because it has a different read size!
#run one slurm file to completion to generate requisite tallymer mappability file
sbatch 3T3_20min_rep1.bigwig.slurm

#after this is done, start all others (and repeat the first one)
for slurm in *bigwig*slurm
do
    sbatch $slurm
done
