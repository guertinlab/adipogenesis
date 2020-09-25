#!/bin/bash

#require slurm headers
for i in *_atac_PE1.fastq.gz
do
    name=$(echo $i | awk -F"/" '{print $NF}' | awk -F"_atac_PE1.fastq.gz" '{print $1}')
    echo $name
    echo '#SBATCH -o' $name'.align.out' > temp.txt
    echo 'i='$i > temp2.txt
    cat align_slurm_header_1.txt temp.txt align_slurm_header_2.txt temp2.txt align_slurm_header_3.txt > $name.align.slurm
    sbatch $name.align.slurm
    rm temp.txt
    rm temp2.txt
done
