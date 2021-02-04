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

mv PSWM_family_1_fimo_all.bed AP1_fimo_all.bed
mv PSWM_family_3_fimo_all.bed GR_fimo_all.bed
mv PSWM_family_5_fimo_all.bed CEBP_fimo_all.bed
mv PSWM_family_18_fimo_all.bed TWIST_fimo_all.bed
rm PSWM*all.bed

#transfer bed files for AP1, GR, CEBP, and TWIST into main figures directory
#check that family number matches up to corresponding motif
mkdir main_figure_beds
cp PSWM_family_1_fimo.bed main_figure_beds/AP1_fimo.bed
cp PSWM_family_3_fimo.bed main_figure_beds/GR_fimo.bed
cp PSWM_family_5_fimo.bed main_figure_beds/CEBP_fimo.bed
cp PSWM_family_18_fimo.bed main_figure_beds/TWIST_fimo.bed

cp PSWM_family_1_2M.bed main_figure_beds/AP1_2M.bed
cp PSWM_family_3_2M.bed main_figure_beds/GR_2M.bed
cp PSWM_family_5_2M.bed main_figure_beds/CEBP_2M.bed
cp PSWM_family_18_2M.bed main_figure_beds/TWIST_2M.bed


#prepare bigWigs for motif enrichment plot
cd /scratch/bhn9by/ATAC/fimo_composites/main_figure_beds

cat SP_unsorted.bed | sort -k1,1 -k2,2n > SP_2M.bed
cat KLF_unsorted.bed | sort -k1,1 -k2,2n > KLF_2M.bed

module load ucsc-tools/3.7.4 gcc/9.2.0 bedtools/2.29.2

intersectBed -loj -a ../../all_peaks.bed -b SP_2M.bed > ../SP_fimo_all.bed
intersectBed -loj -a ../../all_peaks.bed -b KLF_2M.bed > ../KLF_fimo_all.bed

intersectBed -loj -a ../../dynamic_peaks.bed -b SP_2M.bed > SP_fimo.bed
intersectBed -loj -a ../../dynamic_peaks.bed -b KLF_2M.bed > KLF_fimo.bed

for bed in *2M.bed
do
    name=$(echo $bed | awk -F"/" '{print $NF}' | awk -F"_2M.bed" '{print $1}')
    echo $name
    #summing scores of motifs w/in overlapping genomic interval
    cat $bed | mergeBed -i stdin -c 4 -o sum > ${name}_merged_2M.bed
    bedGraphToBigWig ${name}_merged_2M.bed ../../mm10.chrom.sizes ${name}_mm10_instances.bigWig
    rm ${name}_merged_2M.bed
done

#extract up / down composite motifs from combined SP.KLF
cd /scratch/bhn9by/Adipogenesis/ATAC/SP_KLF_split

module load gcc/7.1.0  openmpi/3.1.4 R/4.0.0

Rscript ../extract.motifs.from.combined.family.R

module load gcc/7.1.0 meme/4.10.2

ceqlogo -i SP.KLF_activated.txt -m SP.KLF_activated -o SP.KLF.activated.eps
ceqlogo -i SP.KLF_activated.txt -m SP.KLF_activated -o SP.KLF.activated.rc.eps -r

ceqlogo -i SP.KLF_repressed.txt -m SP.KLF_repressed -o SP.KLF.repressed.eps
ceqlogo -i SP.KLF_repressed.txt -m SP.KLF_repressed -o SP.KLF.repressed.rc.eps -r
