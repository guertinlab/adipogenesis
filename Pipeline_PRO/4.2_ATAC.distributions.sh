module load gcc/7.1.0 bedtools/2.26.0

#generate promoters.bed (150bp upstream & 50bp downstream TSS, strand specific)
cd /scratch/bhn9by/PRO

awk '{OFS="\t";} {$6 == "+";} {print $1,$2,$2+1,$4,$5,$6}' primary_transcript_annotation/primary_transcript_annotation.bed > temp1.bed
awk '{OFS="\t";} {$6 == "-";} {print $1,$3,$3+1,$4,$5,$6}' primary_transcript_annotation/primary_transcript_annotation.bed > temp2.bed
cat temp1.bed temp2.bed > temp3.bed
slopBed -i temp3.bed -g mm10.chrom.sizes -l 150 -r 50 -s > promoters.bed

rm temp1.bed temp2.bed temp3.bed

#dynamic ATAC peaks
cd /scratch/bhn9by/ATAC

rm dynamic_ATAC_peaks_promoters_temp.bed
rm dynamic_ATAC_peaks_intragenic_temp.bed

intersectBed -wa -a dynamic_peaks.bed -b ../PRO/promoters.bed -u > dynamic_ATAC_peaks_promoters_temp.bed

intersectBed -wa -a dynamic_peaks.bed -b ../PRO/primary_transcript_annotation/primary_transcript_annotation.bed -u > dynamic_ATAC_peaks_intragenic_temp.bed

#control expected
tab=$'\t'
rm dynamic.expected.distributions.txt
touch dynamic.expected.distributions.txt
echo "promoters""$tab""intragenic""$tab""intergenic" >> dynamic.expected.distributions.txt

for ((i=1;i<=1000;i++)); 
do 
  echo $i
   
  shuffleBed -i dynamic_peaks.bed -g ../PRO/mm10.chrom.sizes > shuffled.bed

  intersectBed -wa -a shuffled.bed -b ../PRO/promoters.bed -u > shuffled_promoters_temp.bed
  intersectBed -wa -a shuffled.bed -b ../PRO/primary_transcript_annotation/primary_transcript_annotation.bed -u > shuffled_intragenic_temp.bed

  Rscript calc.expected.dynamic.R
  
  rm shuffled.bed
  rm shuffled_promoters_temp.bed
  rm shuffled_intragenic_temp.bed
done

#all ATAC peaks
rm all_ATAC_peaks_promoters_temp.bed
rm all_ATAC_peaks_intragenic_temp.bed

intersectBed -wa -a all_peaks.bed -b ../PRO/promoters.bed -u > all_ATAC_peaks_promoters_temp.bed

intersectBed -wa -a all_peaks.bed -b ../PRO/primary_transcript_annotation/primary_transcript_annotation.bed -u > all_ATAC_peaks_intragenic_temp.bed

#control expected
tab=$'\t'
rm all_expected.distributions.txt
touch all_expected.distributions.txt
echo "promoters""$tab""intragenic""$tab""intergenic" >> all_expected.distributions.txt

for ((i=1;i<=1000;i++)); 
do 
  echo $i
   
  shuffleBed -i all_peaks.bed -g ../PRO/mm10.chrom.sizes > shuffled.bed

  intersectBed -wa -a shuffled.bed -b ../PRO/promoters.bed -u > shuffled_promoters_temp.bed
  intersectBed -wa -a shuffled.bed -b ../PRO/primary_transcript_annotation/primary_transcript_annotation.bed -u > shuffled_intragenic_temp.bed

  Rscript calc.observed.all.R
  
  rm shuffled.bed
  rm shuffled_promoters_temp.bed
  rm shuffled_intragenic_temp.bed
done

Rscript ATAC.distributions.R
