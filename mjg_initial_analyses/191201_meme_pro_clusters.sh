#!/bin/bash
cd ~/pro_seq
for i in cluster_bed_pro_cluster*.bed
do
    name=$(echo $i | awk -F"_bed_pro_" '{print $NF}' | awk -F".bed" '{print $1}')
    echo $name	 
    slopBed -i $i \
	-g ../mm10.chrom.sizes -b 50 | \
        fastaFromBed -fi /scratch/mjg7y/mm10.fa -bed stdin \
		 -fo ${name}_pro.fasta
    head ${name}_pro.fasta
    srun meme -p 59 -oc meme -p 64 -oc ${name}.dREG_meme_output -objfun classic -evt 0.01 -nmotifs 15 -searchsize 0 -minw 6 \
	 -maxw 18 -revcomp -dna -markov_order 3 -maxsize 100000000  \
         ${name}_pro.fasta
done
