#!/bin/bash

module load gcc/9.2.0 bedtools/2.29.2

cd /scratch/bhn9by/ATAC/fimo_composites

for i in *_2M.txt
do
    name=$(echo $i | awk -F"/" '{print $NF}' | awk -F"_2M.txt" '{print $1}')
    echo $name
    intersectBed -loj -a ../dynamic_peaks.bed -b $i > ${name}_fimo.bed
    intersectBed -loj -a ../nondynamic_peaks.bed -b $i > ${name}_fimo_nondyn.bed
    intersectBed -loj -a ../all_peaks.bed -b $i > ${name}_fimo_all.bed
    cat $i | cut -f1-3,5 | sort -k1,1 -k2,2n > ${name}_2M.bed 
done

#transfer bed files for AP1, GR, CEBP, and TWIST into main figures directory
#check that family number matches up to corresponding motif

rm -r main_figure_beds
mkdir main_figure_beds

cp PSWM_family_1_fimo.bed main_figure_beds/AP1_fimo.bed
cp PSWM_family_3_fimo.bed main_figure_beds/GR_fimo.bed
cp PSWM_family_5_fimo.bed main_figure_beds/CEBP_fimo.bed
cp PSWM_family_18_fimo.bed main_figure_beds/TWIST_fimo.bed

cp PSWM_family_1_2M.bed main_figure_beds/AP1_2M.bed
cp PSWM_family_3_2M.bed main_figure_beds/GR_2M.bed
cp PSWM_family_5_2M.bed main_figure_beds/CEBP_2M.bed
cp PSWM_family_18_2M.bed main_figure_beds/TWIST_2M.bed

#rename 2M files for generating final ATAC dataframe
#skip SP/KLF (family 7)
cp PSWM_family_1_fimo.bed BHLH_fimo_all.bed
cp PSWM_family_2_fimo.bed TCF21_fimo_all.bed
cp PSWM_family_3_fimo.bed GR_fimo_all.bed
cp PSWM_family_4_fimo.bed BHLHA15_fimo_all.bed
cp PSWM_family_5_fimo.bed CEBP_fimo_all.bed
cp PSWM_family_6_fimo.bed CTFCL_fimo_all.bed
#cp PSWM_family_7_fimo.bed SPKLF_fimo_all.bed
cp PSWM_family_8_fimo.bed ETS_fimo_all.bed
cp PSWM_family_9_fimo.bed ZBTB33_fimo_all.bed
cp PSWM_family_10_fimo.bed MAZ_fimo_all.bed
cp PSWM_family_11_fimo.bed NFIC_fimo_all.bed
cp PSWM_family_12_fimo.bed NFY_fimo_all.bed
cp PSWM_family_13_fimo.bed NRF_fimo_all.bed
cp PSWM_family_14_fimo.bed NUR77_fimo_all.bed
cp PSWM_family_15_fimo.bed STAT_fimo_all.bed
cp PSWM_family_16_fimo.bed TFAP2A_fimo_all.bed
cp PSWM_family_17_fimo.bed TEAD_fimo_all.bed
cp PSWM_family_18_fimo.bed TWIST_fimo_all.bed
cp PSWM_family_19_fimo.bed ZNF263_fimo_all.bed








