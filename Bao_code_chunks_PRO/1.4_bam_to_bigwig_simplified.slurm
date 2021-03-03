#!/bin/bash
#SBATCH -n 1
#SBATCH -t 96:00:00
#SBATCH -o big_to_bigwig.out
#SBATCH -p standard
#SBATCH -A guertinlab

module load genometools/1.5.10 wigtobigwig/2.8 gcc/7.1.0 seqoutbias/1.2.0

#perform seqOutBias and convert .bam to .bigwig
cd /scratch/bhn9by/PRO

for bam in *_pro_plus.bam
do
    name=$(echo $bam | awk -F"_pro_plus.bam" '{print $1}')
    echo $name
    
    #tallymer mappability file should already be available from ATAC pipeline so conversion should be fast
    seqOutBias /project/genomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa ${name}_pro_plus.bam --no-scale --skip-bed \
              --bw=${name}_plus_body_0-mer.bigWig --tail-edge --read-size=30
    seqOutBias /project/genomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa ${name}_pro_minus.bam --no-scale --skip-bed \
              --bw=${name}_minus_body_0-mer.bigWig --tail-edge --read-size=30
done      
        
