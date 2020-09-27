#!/bin/bash
#SBATCH -n 1
#SBATCH -t 96:00:00
#SBATCH -o meme.enrichment.out
#SBATCH -p largemem
#SBATCH -A guertinlab

module load gcc/7.1.0 mvapich2/2.3.1 bedtools/2.26.0 meme/5.1.0

wget https://hgdownload-test.gi.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes
mkdir meme_motif_enrichment
cd meme_motif_enrichment

#run MEME on each cluster to find enriched motifs
for i in ../*cluster*bed
do
    name=$(echo $i | awk -F"_" '{print $NF}' | awk -F".bed" '{print $1}')
    echo $name

    slopBed -i ../$i \
      -g mm10.chrom.sizes -b -50 | \
        fastaFromBed -fi /project/genomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa -bed stdin \
              -fo ${name}.fasta
    head ${name}.fasta
    srun meme -p 59 -oc ${name}.meme_output -nmotifs 15 -objfun classic -evt 0.01 -searchsize 0 -minw 6 \
         -maxw 18 -revcomp -dna -markov_order 3 -maxsize 100000000  \
         ${name}.fasta
done
