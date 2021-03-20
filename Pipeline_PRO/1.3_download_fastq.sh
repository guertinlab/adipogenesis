#!/bin/bash

#upload sra.metadata.csv and SRR_Acc_List.txt to Rivanna
cd /scratch/bhn9by/PRO

#remove DOS \r\n\ artifact from .csv (if applicable)
sed -i 's/\r$//' sra.metadata.csv

#download .fastq
#make a unique slurm file for each replicate and run them in parallel
cat SRR_Acc_List.txt | while read acc
do
    echo $acc
    echo '#SBATCH -o' $acc'.out' > temp.txt
    echo 'fasterq-dump' $acc > temp2.txt
    echo 'gzip' $acc'.fastq' > temp3.txt
    cat sra_slurm_header_1.txt temp.txt sra_slurm_header_2.txt temp2.txt temp3.txt > $acc.slurm
    sbatch $acc.slurm
    rm temp.txt
    rm temp2.txt
    rm temp3.txt
done

#Alternatively, you can run fasterq-dump and gzip in series.
module load sratoolkit
cat SRR_Acc_List.txt | while read acc
do
    echo $acc
    fasterq-dump $acc
    gzip $acc.fastq
done
module purge

#After all jobs are done, rename files to actual sample names
for fq in SRR*.fastq.gz
do
    name=$(echo $fq | awk -F".fastq.gz" '{print $1}')
    echo $name
    line=$(grep $name sra.metadata.csv)
    treat=$(echo $line | awk -F',' '{print $31}')
    echo $treat
    mv $fq $treat
done

