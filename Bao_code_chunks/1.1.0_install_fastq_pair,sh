#!/bin/bash

#Install fastq_pair command for alignment
#Technically for this ATAC pipeline you don't need it but we included it for consistency

mkdir /scratch/bhn9by/ATAC
cd /scratch/bhn9by/ATAC

#Retrieve and build fastq_pair binary
wget https://github.com/linsalrob/fastq-pair/archive/master.zip
unzip master.zip
cd fastq-pair-master
gcc -std=gnu99 main.c robstr.c fastq_pair.c is_gzipped.c -o fastq_pair

#put fastq_pair binary on the $PATH
cp fastq_pair /home/bhn9by/bin
cd..
