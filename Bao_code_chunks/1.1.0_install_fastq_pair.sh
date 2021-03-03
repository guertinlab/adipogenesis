#!/bin/bash

#on Rivanna
#Install fastq_pair and put on $PATH
#fastq_pair is not required for adipogenesis project but we use it redundantly for our ATAC pipeline

wget https://github.com/linsalrob/fastq-pair/archive/master.zip
unzip master.zip
cd fastq-pair-master
gcc -std=gnu99 main.c robstr.c fastq_pair.c is_gzipped.c -o fastq_pair
cp fastq_pair /home/bhn9by/bin
cd..
