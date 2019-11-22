#!/bin/bash
cd ~/clusters
for i in cluster_bed_cluster1.bed
do
    name=$(echo $i | awk -F"_bed" '{print $NF}' | awk -F".bed" '{print $1}')
    echo $name	 
    slopBed -i cluster_bed${name}.bed \
	-g ../mm10.chrom.sizes -b -50 | \
    fastaFromBed -fi /scratch/mjg7y/mm10.fa -bed stdin \
		 -fo ${name}.fasta
    head ${name}.fasta
    meme -oc ${name}.meme_output -objfun classic -nmotifs 5 -searchsize 0 -minw 6 \
	 -maxw 18 -revcomp -dna -markov_order 3 -maxsize 100000000  \
         ${name}.fasta
done
