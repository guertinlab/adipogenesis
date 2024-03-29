#!/bin/bash

module load bioconda ucsc-tools/3.7.4 gcc/7.1.0 openmpi/3.1.4 R/4.0.0

cd /scratch/bhn9by/PRO

# should not run wget in slurm--download connection is often severed
wget https://hgdownload-test.gi.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz
gunzip gencode.vM25.annotation.gtf.gz

# get the first exons for protein coding genes
grep 'transcript_type "protein_coding"' gencode.vM25.annotation.gtf | \
     awk '{if($3=="exon"){print $0}} ' | \
     grep -w "exon_number 1" | \
     cut -f1,4,5,7,9 | tr ";" "\t" | \
     awk '{for(i=5;i<=NF;i++){if($i~/^gene_name/){a=$(i+1)}} print $1,$2,$3,a,"na",$4}' | \
     tr " " "\t" | tr -d '"' > gencode.mm10.firstExon.bed

# get all transcripts for protein coding genes
grep 'transcript_type "protein_coding"' gencode.vM25.annotation.gtf | \
     awk '{if($3=="transcript"){print $0}} ' | \
     cut -f1,4,5,7,9 | tr ";" "\t" | \
     awk '{for(i=5;i<=NF;i++){if($i~/^gene_name/){a=$(i+1)}} print $1,$2,$3,a,"na",$4}' | \
     tr " " "\t" | tr -d '"' > gencode.mm10.transcript.bed

# chrom sizes file
chromsizes=mm10.chrom.sizes
# merge plus
plusfiles=$(ls *plus*bigWig)
bigWigMerge ${plusfiles} pro_plus_merged.bedGraph
LC_COLLATE=C sort -k1,1 -k2,2n pro_plus_merged.bedGraph > pro_plus_merged_sorted.bedGraph
bedGraphToBigWig pro_plus_merged_sorted.bedGraph ${chromsizes} pro_plus_merged.bigWig
# merge minus
minusfiles=$(ls *minus*bigWig)
bigWigMerge ${minusfiles} pro_minus_merged.bedGraph
LC_COLLATE=C sort -k1,1 -k2,2n pro_minus_merged.bedGraph > pro_minus_merged_sorted.bedGraph
bedGraphToBigWig pro_minus_merged_sorted.bedGraph ${chromsizes} pro_minus_merged.bigWig

# this pTA Rscript takes 1.5hr to compile--run using slurm for convenience
Rscript primary_transcript_annotation.R

# discard 6day timepoint from rest of analysis
rm -r 6day
mkdir 6day
mv *6d* 6day



