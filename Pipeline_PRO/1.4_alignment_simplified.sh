#!/bin/bash

module load bioconda/py3.6 gcc/7.1.0 bowtie2/2.2.9 samtools/1.10
source activate myenv

#align .fastq to mm10 genome
for i in *_pro.fastq.gz
do
    name=$(echo $i | awk -F"/" '{print $NF}' | awk -F"_pro.fastq" '{print $1}')
    echo $name

    echo 'trim adapters'
    cutadapt -m 26 -a TGGAATTCTCGGGTGCCAAGG $i | \
        fqdedup -i - -o - | \
        fastx_trimmer -Q33 -f 9 -l 38 | \
        fastx_reverse_complement -Q33 -z -o ${name}.processed.fastq.gz

#updated version of cutadapt commnad: 'cutadapt -a TGGAATTCTCGGGTGCCAAGG -m 5 -O 1 $i | \'
#fqdedup -i - -o -

    echo 'align to mouse genome'
    bowtie2 -p 3 -x /project/genomes/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome -U ${name}.processed.fastq.gz | \
        samtools view -b - | \
        samtools sort - -o ${name}.sorted.bam

    echo 'sort and separate plus/minus reads'
    samtools view -bh -F 20 ${name}.sorted.bam > ${name}_pro_plus.bam
    samtools view -bh -f 0x10 ${name}.sorted.bam > ${name}_pro_minus.bam

    echo 'done'
done
