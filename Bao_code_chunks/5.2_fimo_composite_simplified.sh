#!/bin/bash

#run fimo against PSWM and take top 2 million hits
rm -r fimo_composites
mkdir fimo_composites

for i in PSWM*txt
do
    name=$(echo $i | awk -F".txt" '{print $1}')
    echo $name
    
    module load gcc/7.1.0  meme/4.10.2

    i=../composite_motifs/$name/{name}_meme.txt
    cd /scratch/bhn9by/ATAC/fimo_composites

    name=$(echo $i | awk -F"/" '{print $NF}'| awk -F"_meme.txt" '{print $1}')

    fimo --thresh 0.001 --text $i /project/genomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa > ${name}_composite_fimo.txt

    #this takes top 2M
    score=$(tail -n +2 ${name}_composite_fimo.txt | sort -nrk6,6 | awk 'FNR == 2000000 {print $6}')
    echo $score
    tail -n +2 ${name}_composite_fimo.txt | awk -v sco="$score" '{ if ($6 >= sco) { print } }' | awk '{OFS="\t";} {print $2,$3,$4,$7,$6,$5,$8}' > ${name}_2M.txt


    #this was to get the order of conformity to consensus.
    tomtom -no-ssc -oc ${name}_ranks.tomtom_output -verbosity 1 -incomplete-scores -min-overlap 1 -dist ed -evalue -thresh 0.05 ../all_query_factors_meme.txt ${name}_composite_fimo.txt
done
  
