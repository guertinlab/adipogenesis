#!/bin/bash

#run fimo against PSWM and take top 2 million hits
rm -r fimo_composites
mkdir fimo_composites

for i in PSWM*txt
do
    name=$(echo $i | awk -F".txt" '{print $1}')
    echo $name
    echo '#SBATCH -o' $name'.fimo.out' > temp.txt
    echo 'i='../composite_motifs/$name/${name}_meme.txt > temp2.txt
    cat fimo_slurm_header_1.txt temp.txt fimo_slurm_header_2.txt temp2.txt fimo_slurm_header_3.txt > $name.fimo.slurm
    sbatch $name.fimo.slurm                                                                                                                                                           
    rm temp.txt
    rm temp2.txt
done
