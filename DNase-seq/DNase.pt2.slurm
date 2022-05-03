#!/bin/bash
#SBATCH -n 1
#SBATCH -t 96:00:00
#SBATCH -p standard
#SBATCH -A janeslab
#SBATCH -o DNase.pt2.out

echo 'Calling peaks'
module purge
module load gcc/7.1.0 bedtools/2.26.0 macs2/2.2.7.1

macs2 callpeak -t *_rmdup.bam -f BAM -n TERT4 --outdir TERT4_macs --keep-dup 50 -q 0.05
slopBed -i TERT4_macs/TERT4_summits.bed -g hg38.chrom.sizes -b 100 | intersectBed -v -a stdin -b hg38-blacklist.v2.bed > TERT4_DNase_peaks.bed

echo 'Converting .bams to .bigWigs'
module purge
module load genometools/1.5.10 wigtobigwig/2.8 gcc/7.1.0 seqoutbias/1.2.0

for bam in *rmdup.bam
do
    name=$(echo $bam | awk -F"_rmdup.bam" '{print $1}')
    echo $name
    
    #if you are running seqOutBias for the first time, the majority will be spent generating the tallymer mappability file
    #determined read length from .bam w/ samtools view TERT4_4h_rep2_rmdup.bam | awk '{print length($10)}' | head -1000 | sort -u
    seqOutBias /project/genomes/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa $bam --skip-bed --no-scale --bw=${name}.bigWig --shift-counts --read-size=50
done



