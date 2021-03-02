#!/bin/bash

module load bioconda/py3.6 gcc/7.1.0 bowtie2/2.2.9 samtools/1.10
source activate myenv

#align .fastq to mm10 genome
for i in *_atac_PE1.fastq.gz
do
    name=$(echo $i | awk -F"/" '{print $NF}' | awk -F"_atac_PE1.fastq.gz" '{print $1}')
    echo $name
    gunzip $name*.gz
    
    echo 'pairing .fastq files'
    fastq_pair ${name}_atac_PE1.fastq ${name}_atac_PE2.fastq
    rm $name*single.fq
    
    echo 'align to mouse genome'
    bowtie2 --maxins 500 -x /project/genomes/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome \
      -1 ${name}_atac_PE1.fastq.paired.fq \
      -2 ${name}_atac_PE2.fastq.paired.fq -S ${name}_atac_smp.sam
      
    echo 'quality filter and remove duplicate amplicons'	
    samtools view -b -q 10 ${name}_atac_smp.sam | samtools sort -n - | \
        samtools fixmate -m - - | samtools sort - | samtools markdup -r - ${name}_atac_rmdup.bam
    rm ${name}_atac_smp.sam
    gzip ${name}*fastq
done
