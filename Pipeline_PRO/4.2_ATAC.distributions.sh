module load gcc/7.1.0 bedtools/2.26.0

#generate promoters.bed (150bp upstream & 50bp downstream TSS, strand specific)
cd /scratch/bhn9by/PRO
rm promoters.bed

awk '{OFS="\t";} {$6 == "+";} {print $1,$2,$2+1,$4,$5,$6}' primary_transcript_annotation/primary_transcript_annotation.bed > temp1.bed
awk '{OFS="\t";} {$6 == "-";} {print $1,$3,$3+1,$4,$5,$6}' primary_transcript_annotation/primary_transcript_annotation.bed > temp2.bed
cat temp1.bed temp2.bed > temp3.bed
slopBed -i temp3.bed -g mm10.chrom.sizes -l 150 -r 50 -s > promoters.bed

rm temp1.bed temp2.bed temp3.bed

#all ATAC peaks
rm all_ATAC_peaks_promoters_temp.bed
rm all_ATAC_peaks_intragenic_temp.bed

intersectBed -wa -a all_peaks.bed -b ../PRO/promoters.bed -u > all_ATAC_peaks_promoters_temp.bed

intersectBed -wa -a all_peaks.bed -b ../PRO/primary_transcript_annotation/primary_transcript_annotation.bed -u > all_ATAC_peaks_intragenic_temp.bed

#finalize promoters and intragenic bed files
Rscript ATAC.distributions.R
