#!/bin/bash
#SBATCH -n 1
#SBATCH -t 96:00:00
#SBATCH -o PSWM.fimo.out
#SBATCH -p standard
#SBATCH -A guertinlab

module load gcc/7.1.0  meme/4.10.2

#run FIMO against PSWM and take top 2 million hits
rm -r fimo_composites
mkdir fimo_composites

for i in PSWM*txt
do
    name=$(echo $i | awk -F".txt" '{print $1}')
    echo $name
    
    #run FIMO
    cd /scratch/bhn9by/ATAC/fimo_composites
    fimo --thresh 0.001 --text ../composite_motifs/$name/{name}_meme.txt /project/genomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa > ${name}_composite_fimo.txt

    #this takes top 2M
    score=$(tail -n +2 ${name}_composite_fimo.txt | sort -nrk6,6 | awk 'FNR == 2000000 {print $6}')
    echo $score
    tail -n +2 ${name}_composite_fimo.txt | awk -v sco="$score" '{ if ($6 >= sco) { print } }' | awk '{OFS="\t";} {print $2,$3,$4,$7,$6,$5,$8}' > ${name}_2M.txt

    #this was to get the order of conformity to consensus.
    tomtom -no-ssc -oc ${name}_ranks.tomtom_output -verbosity 1 -incomplete-scores -min-overlap 1 -dist ed -evalue -thresh 0.05 ../all_query_factors_meme.txt ${name}_composite_fimo.txt
done
  
