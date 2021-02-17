#!/bin/bash

cd /scratch/bhn9by/ATAC/fimo_composites

mkdir supp_figure_beds

#transfer bed files for supplemental factors
#check that family number matches up to corresponding motif
cp PSWM_family_12_fimo.bed supp_figure_beds/NFY_fimo.bed
cp PSWM_family_13_fimo.bed supp_figure_beds/NRF_fimo.bed
cp PSWM_family_15_fimo.bed supp_figure_beds/STAT_fimo.bed
cp PSWM_family_16_fimo.bed supp_figure_beds/TFAP2_fimo.bed
cp PSWM_family_17_fimo.bed supp_figure_beds/TEAD_fimo.bed
cp PSWM_family_2_fimo.bed supp_figure_beds/bHLH_fimo.bed
cp PSWM_family_6_fimo.bed supp_figure_beds/CTCFL_fimo.bed
cp PSWM_family_8_fimo.bed supp_figure_beds/ELF_fimo.bed

cp PSWM_family_12_2M.bed supp_figure_beds/NFY_2M.bed
cp PSWM_family_13_2M.bed supp_figure_beds/NRF_2M.bed
cp PSWM_family_15_2M.bed supp_figure_beds/STAT_2M.bed
cp PSWM_family_16_2M.bed supp_figure_beds/TFAP2_2M.bed
cp PSWM_family_17_2M.bed supp_figure_beds/TEAD_2M.bed
cp PSWM_family_2_2M.bed supp_figure_beds/bHLH_2M.bed
cp PSWM_family_6_2M.bed supp_figure_beds/CTCFL_2M.bed
cp PSWM_family_8_2M.bed supp_figure_beds/ELF_2M.bed

#prepare bigWigs for motif enrichment plot
cd /scratch/bhn9by/ATAC/fimo_composites/supp_figure_beds

module load ucsc-tools/3.7.4 gcc/9.2.0 bedtools/2.29.2

for bed in *2M.bed
do
    name=$(echo $bed | awk -F"/" '{print $NF}' | awk -F"_2M.bed" '{print $1}')
    echo $name
    #summing scores of motifs w/in overlapping genomic intervals
    cat $bed | mergeBed -i stdin -c 4 -o sum > ${name}_merged_2M.bed
    bedGraphToBigWig ${name}_merged_2M.bed ../../mm10.chrom.sizes ${name}_mm10_instances.bigWig
done

#Rscripts for supplement are slightly revised from main figure 1
module load gcc/7.1.0  openmpi/3.1.4 R/4.0.0

Rscript /scratch/bhn9by/ATAC/post.composite.fimo.supp.R
Rscript /scratch/bhn9by/ATAC/plot.motif.enrichment.supp.R
