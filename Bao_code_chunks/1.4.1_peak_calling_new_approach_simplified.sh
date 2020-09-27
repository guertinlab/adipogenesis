#!/bin/bash

#archive all files for 6d treatment; not used for analysis
mkdir 6day
mv *6d* 6day/

#merge .bam, call peaks, remove blacklisted regions, set summit windows
for bam in *rep1_atac_rmdup.bam
do
    name=$(echo $bam | awk -F"_rep1_atac_rmdup.bam" '{print $1}')
    echo $name
    
    echo 'Merging'
    samtools merge ${name}_merged.bam ${name}_rep1_atac_rmdup.bam ${name}_rep2_atac_rmdup.bam ${name}_rep3_atac_rmdup.bam

    echo 'Calling Peaks'
    macs2 callpeak -t ${name}_merged.bam -f BAMPE -n ${name} --outdir ${name}_macs -g mm -B --call-summits --keep-dup 50 -q 0.05

    echo 'Removing blacklisted regions'
    bedtools subtract -a ${name}_macs/${name}_summits.bed -b mm10.blacklist.bed > \
       ${name}_macs/${name}_summits_bl_removed.bed

    echo 'Setting summit windows'
    awk '{OFS="\t";} {print $1,$2-99,$3+100,$4,$5}' ${name}_macs/${name}_summits_bl_removed.bed > \
        ${name}_summit_window.bed
done

#merge called peaks from individual time points into one file
cat *_summit_window.bed | sort -k1,1 -k2,2n | awk ' $2 >= 0 ' | mergeBed -i stdin > merged_peaks.bed
