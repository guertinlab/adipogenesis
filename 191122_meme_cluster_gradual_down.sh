#!/bin/bash
cd ~/clusters/final_clusters
for i in final_cluster_bed_gradual.down.bed
do
    name=$(echo $i | awk -F"_bed_" '{print $NF}' | awk -F".bed" '{print $1}')
    echo $name	 
    slopBed -i $i \
	-g ~/mm10.chrom.sizes -b -50 | \
    fastaFromBed -fi /scratch/mjg7y/mm10.fa -bed stdin \
		 -fo final_cluster_${name}.fasta
    head final_cluster_${name}.fasta
    meme -oc final_cluster_${name}.meme_output -objfun classic -nmotifs 5 -searchsize 0 -minw 6 \
	 -maxw 18 -revcomp -dna -markov_order 3 -maxsize 100000000  \
         final_cluster_${name}.fasta
done
