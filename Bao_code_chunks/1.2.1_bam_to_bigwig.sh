#!/bin/bash

#perform seqOutBias and convert .bam to .bigwig

#generate a unique slurm file for each replicate and run them in parallel:
#header_1   --> sbatch settings
#temp.txt   --> name of .out file
#header_2   --> modules to load
#temp2.txt  --> name of relevant file
#header_3   --> actual commands
for bam in *rmdup.bam
do
    name=$(echo $bam | awk -F"_atac_rmdup.bam" '{print $1}')
    echo $name
    echo '#SBATCH -o' $name'.bigwig.out' > temp.txt
    echo 'bam='$bam > temp2.txt
    cat bigwig_slurm_header_1.txt temp.txt bigwig_slurm_header_2.txt temp2.txt bigwig_slurm_header_3.txt > $name.convert.to.bigwig.slurm
done

#run one slurm file to completion to generate requisite tallymer mappability file
sbatch 3T3_20min_rep1.convert.to.bigwig.slurm

#after this is done, start all others (and repeat the first one)
for slurm in *bigwig*slurm
do
    sbatch $slurm
done
