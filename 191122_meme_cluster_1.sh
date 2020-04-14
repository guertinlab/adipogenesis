#!/bin/bash
cd ~/clusters
for i in cluster_bed_cluster1.bed
do
    name=$(echo $i | awk -F"_bed_" '{print $NF}' | awk -F".bed" '{print $1}')
    echo $name	 
    slopBed -i cluster_bed_${name}.bed \
	-g ../mm10.chrom.sizes -b -50 | \
    fastaFromBed -fi /scratch/mjg7y/mm10.fa -bed stdin \
		 -fo ${name}.fasta
    head ${name}.fasta
    srun meme -p 59 -oc ${name}.meme_output -nmotifs 15 -objfun classic -evt 0.01 -searchsize 0 -minw 6 \
	 -maxw 18 -revcomp -dna -markov_order 3 -maxsize 100000000  \
         ${name}.fasta
done
