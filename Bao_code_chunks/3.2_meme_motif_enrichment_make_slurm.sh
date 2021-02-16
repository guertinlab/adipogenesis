#!/bin/bash

#run MEME on each cluster to find enriched motifs
wget https://hgdownload-test.gi.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes
mkdir meme_motif_enrichment

#generate a unique slurm file for each replicate and run them in parallel:
#header_1   --> sbatch settings
#temp.txt   --> name of .out file
#header_2   --> more sbatch settings and modules to load
#temp2.txt  --> name of relevant file
#header_3   --> actual commands
for i in *cluster*bed
do
    name=$(echo $i | awk -F"_" '{print $NF}' | awk -F".bed" '{print $1}')
    echo $name
    echo '#SBATCH -o' $name'.meme.out' > temp.txt
    echo 'i='$i > temp2.txt
    cat meme_slurm_header_1.txt temp.txt meme_slurm_header_2.txt temp2.txt meme_slurm_header_3.txt > $name.meme.slurm
    sbatch $name.meme.slurm                                                                                                                                                           
    rm temp.txt
    rm temp2.txt
done
