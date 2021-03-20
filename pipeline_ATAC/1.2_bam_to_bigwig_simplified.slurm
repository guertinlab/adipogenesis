#!/bin/bash
#SBATCH -n 1
#SBATCH -t 96:00:00
#SBATCH -o big_to_bigwig.out
#SBATCH -p standard
#SBATCH -A guertinlab

module load genometools/1.5.10 wigtobigwig/2.8 gcc/7.1.0 seqoutbias/1.2.0

#perform seqOutBias and convert .bam to .bigwig
cd /scratch/bhn9by/ATAC

for bam in *rmdup.bam
do
    name=$(echo $bam | awk -F"_atac_rmdup.bam" '{print $1}')
    echo $name
    
    #if you are running seqOutBias for the first time, the majority will be spent generating the tallymer mappability file
    seqOutBias /project/genomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa \
         $bam --skip-bed --no-scale --bw=${name}.bigWig --only-paired --shift-counts --read-size=38
done
