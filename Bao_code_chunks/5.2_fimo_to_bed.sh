#!/bin/bash

module load gcc/9.2.0 bedtools/2.29.2

cd /scratch/abd3x/Adipogenesis/ATAC/fimo_composites

for i in *_2M.txt
do
    name=$(echo $i | awk -F"/" '{print $NF}' | awk -F"_2M.txt" '{print $1}')
    echo $name
    intersectBed -loj -a ../dynamic_peaks.bed -b $i > ${name}_fimo.bed
    intersectBed -loj -a ../nondynamic_peaks.bed -b $i > ${name}_fimo_nondyn.bed
    #intersectBed -loj -a ../all_peaks.bed -b $i > ${name}_fimo_all.bed #is this necessary???
    cat $i | cut -f1-3,5 | sort -k1,1 -k2,2n > ${name}_2M.bed 
done

mv PSWM_family_1_fimo_all.bed AP1_fimo_all.bed
mv PSWM_family_3_fimo_all.bed GR_fimo_all.bed
mv PSWM_family_5_fimo_all.bed CEBP_fimo_all.bed
mv PSWM_family_18_fimo_all.bed TWIST_fimo_all.bed
rm PSWM*all.bed

#transfer bed files for AP1, GR, CEBP, and TWIST into main figures directory
#check that family number matches up to corresponding motif
mkdir main_figure_beds
cp PSWM_family_1_fimo.bed main_figure_beds/AP1_fimo.bed
#cp PSWM_family_2_fimo.bed main_figure_beds/bHLH_fimo.bed
cp PSWM_family_3_fimo.bed main_figure_beds/GR_fimo.bed
cp PSWM_family_5_fimo.bed main_figure_beds/CEBP_fimo.bed
#cp PSWM_family_14_fimo.bed main_figure_beds/NFY_fimo.bed
#cp PSWM_family_15_fimo.bed main_figure_beds/NRF_fimo.bed
cp PSWM_family_18_fimo.bed main_figure_beds/TWIST_fimo.bed

cp PSWM_family_1_2M.bed main_figure_beds/AP1_2M.bed
#cp PSWM_family_2_2M.bed main_figure_beds/bHLH_2M.bed
cp PSWM_family_3_2M.bed main_figure_beds/GR_2M.bed
cp PSWM_family_5_2M.bed main_figure_beds/CEBP_2M.bed
#cp PSWM_family_14_2M.bed main_figure_beds/NFY_2M.bed
#cp PSWM_family_15_2M.bed main_figure_beds/NRF_2M.bed
cp PSWM_family_18_2M.bed main_figure_beds/TWIST_2M.bed


