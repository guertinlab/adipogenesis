#!/bin/bash

cd /Volumes/GUERTIN_2/adipogenesis/atac
cp */*_sample_atac.bam ./
mkdir 6day
mv 3T3_6d_*bam 6day/
samtools merge preadipMerged.bam *.bam


# isolate only concordant alignments
mv preadipMerged.bam pre.bam

samtools view -b -f 0x2 pre.bam -o preadipMerged.bam
samtools sort -n preadipMerged.bam -o sorted.bam
samtools fixmate sorted.bam fixed.bam

#convert to bed 
bedtools bamtobed -i fixed.bam -bedpe > preadip0.bed

#take columns that span the PE1 start and PE2 start
awk '{OFS="\t";} {print $1,$2,$6,$7,$8,$9}' preadip0.bed > bed
sort -k 1,1 bed > sorted.bed
#sort -k1,1 -k2,2n 

# from bed to bedgraph
bedtools genomecov -bg -trackline -trackopts name=preadip -i sorted.bed -g mm10.chrom.sizes > preadip.bedGraph

# generate the bigwig
#module load ucsc-tools
bedGraphToBigWig preadip.bedGraph mm10.chrom.sizes preadip.bigWig

#https://data.cyverse.org/dav-anon/iplant/home/guertin/preadip.bigWig
#https://genome.ucsc.edu/s/Mike%20Guertin/mm10_0%2D4hr_merged

for fq in *rep1_atac_sample_atac.bam
do
    name=$(echo $fq | awk -F"_rep1_atac_sample_atac.bam" '{print $1}')
    echo $name
    repfiles=$(ls ${name}*_atac_sample_atac.bam)
    samtools merge ${name}_merged.rmdup.bam $repfiles
    samtools view -b -f 0x2 ${name}_merged.rmdup.bam -o ${name}_merged.rmdup.1.bam
    samtools sort -n ${name}_merged.rmdup.1.bam -o ${name}_merged.rmdup.sorted.bam
    samtools fixmate ${name}_merged.rmdup.sorted.bam ${name}_merged.rmdup.fixed.bam
#convert to bed 
    bedtools bamtobed -i ${name}_merged.rmdup.fixed.bam -bedpe > ${name}_merged.rmdup.fixed.bed
#take columns that span the PE1 start and PE2 start
    awk '{OFS="\t";} {print $1,$2,$6,$7,$8,$9}' ${name}_merged.rmdup.fixed.bed > ${name}_merged.rmdup.final.bed
    sort -k1,1 -k2,2n ${name}_merged.rmdup.final.bed > ${name}_merged.rmdup.final.sorted.bed
    # from bed to bedgraph
    reads=$(wc -l ${name}_merged.rmdup.final.sorted.bed | awk -F" " '{print $1}')
#three hard coded    
    norm=$(echo 10000000/$reads | bc -l | xargs printf "%.*f\n" 3)
    genomeCoverageBed -bg -scale ${norm} -trackline -trackopts 'track type=bedGraph name='"${name}"' alwaysZero=on visibility=full' -i ${name}_merged.rmdup.final.sorted.bed -g mm10.chrom.sizes > ${name}.bedGraph
    bedGraphToBigWig ${name}.bedGraph mm10.chrom.sizes ${name}.merged.bigWig
done



#6 day
cd 6day
samtools merge sixdayMerged.bam *.bam
mv sixdayMerged.bam pre.bam

samtools view -b -f 0x2 pre.bam -o sixdayMerged.bam
samtools sort -n sixdayMerged.bam -o sorted.bam
samtools fixmate sorted.bam fixed.bam

#convert to bed 
bedtools bamtobed -i fixed.bam -bedpe > sixdayMerged.bed

#take columns that span the PE1 start and PE2 start
awk '{OFS="\t";} {print $1,$2,$6,$7,$8,$9}' sixdayMerged.bed > bed
sort -k 1,1 bed > sorted.bed
#sort -k1,1 -k2,2n 

# from bed to bedgraph
bedtools genomecov -bg -trackline -trackopts name=sixDay -i sorted.bed -g mm10.chrom.sizes > sixdayMerged.bedGraph

# generate the bigwig
#module load ucsc-tools
bedGraphToBigWig sixdayMerged.bedGraph mm10.chrom.sizes sixdayMerged.bigWig


#I ran this for the UCSC genome browser
for fq in *rep1_atac_sample_atac.bam
do
    name=$(echo $fq | awk -F"_rep1_atac_sample_atac.bam" '{print $1}')
    echo $name
    repfiles=$(ls ${name}*_atac_sample_atac.bam)
    samtools merge ${name}_merged.rmdup.bam $repfiles
    samtools view -b -f 0x2 ${name}_merged.rmdup.bam -o ${name}_merged.rmdup.1.bam
    rm ${name}_merged.rmdup.bam
    samtools sort -n ${name}_merged.rmdup.1.bam -o ${name}_merged.rmdup.sorted.bam
    rm ${name}_merged.rmdup.1.bam
    samtools fixmate ${name}_merged.rmdup.sorted.bam ${name}_merged.rmdup.fixed.bam
#convert to bed 
    bedtools bamtobed -i ${name}_merged.rmdup.fixed.bam -bedpe > ${name}_merged.rmdup.fixed.bed
#take columns that span the PE1 start and PE2 start
    awk '{OFS="\t";} {print $1,$2,$6,$7,$8,$9}' ${name}_merged.rmdup.fixed.bed > ${name}_merged.rmdup.final.bed
    rm ${name}_merged.rmdup.fixed.bed
    sort -k1,1 -k2,2n ${name}_merged.rmdup.final.bed > ${name}_merged.rmdup.final.sorted.bed
    # from bed to bedgraph
    reads=$(wc -l ${name}_merged.rmdup.final.sorted.bed | awk -F" " '{print $1}')
#three hard coded    
    norm=$(echo 10000000/$reads | bc -l | xargs printf "%.*f\n" 3)
    genomeCoverageBed -bg -scale ${norm} -trackline -trackopts 'track type=bedGraph name='"${name}"' alwaysZero=on visibility=full' -i ${name}_merged.rmdup.final.sorted.bed -g ../mm10.chrom.sizes > ${name}.bedGraph
    bedGraphToBigWig ${name}.bedGraph ../mm10.chrom.sizes ${name}.merged.bigWig
done

#I loaded the following tracks on Cyverse:

#https://data.cyverse.org/dav-anon/iplant/home/guertin/adipo_hub/hub.txt
#https://data.cyverse.org/dav-anon/iplant/home/guertin/adipo_hub/genomes.txt
#https://data.cyverse.org/dav-anon/iplant/home/guertin/adipo_hub/mm10/trackDb.txt
#https://data.cyverse.org/dav-anon/iplant/home/guertin/adipo_hub/mm10/3T3_t0.merged.bigWig
#https://data.cyverse.org/dav-anon/iplant/home/guertin/adipo_hub/mm10/3T3_20min.merged.bigWig
#https://data.cyverse.org/dav-anon/iplant/home/guertin/adipo_hub/mm10/3T3_40min.merged.bigWig
#https://data.cyverse.org/dav-anon/iplant/home/guertin/adipo_hub/mm10/3T3_60min.merged.bigWig
#https://data.cyverse.org/dav-anon/iplant/home/guertin/adipo_hub/mm10/3T3_2hr.merged.bigWig
#https://data.cyverse.org/dav-anon/iplant/home/guertin/adipo_hub/mm10/3T3_3hr.merged.bigWig
#https://data.cyverse.org/dav-anon/iplant/home/guertin/adipo_hub/mm10/3T3_4hr.merged.bigWig
#https://data.cyverse.org/dav-anon/iplant/home/guertin/adipo_hub/mm10/3T3_6d.merged.bigWig


#the UCSC session is saved here:
#https://genome.ucsc.edu/s/Mike%20Guertin/mm10_ATAC_seq_hub
mkdir 4hour
cd 4hour
#peak calling
fdr=0.05
files=$(ls *_atac_sample_atac.bam)
name=3T3_atac
species=mm
ndups=50
macs2 callpeak -t ${files} -f BAMPE -n ${name} --outdir atac4h_${fdr} -g ${species} -B --call-summits --keep-dup ${ndups} -q ${fdr}

cd ../6day

fdr=0.05
files=$(ls *_atac_sample_atac.bam)
name=3T3_atac_6day
species=mm
ndups=50
macs2 callpeak -t ${files} -f BAMPE -n ${name} --outdir atac4h_${fdr} -g ${species} -B --call-summits --keep-dup ${ndups} -q ${fdr}

cd ../
# get blacklisted regions
wget http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/mm10-mouse/mm10.blacklist.bed.gz
gunzip mm10.blacklist.bed.gz

# document summit coordinate
cd 4hour/atac4h_0.05/
bedtools subtract -a 3T3_atac_summits.bed -b ../mm10.blacklist.bed > 3T3_atac_summits_bl_removed.bed
awk '{OFS="\t";} {print $1,$2-99,$3+100,$4,$5}' 3T3_atac_summits_bl_removed.bed > 3T3_atac_summit_200window.bed

cd ../../6day/atac6d_0.05
bedtools subtract -a 3T3_atac_6day_summits.bed -b ../mm10.blacklist.bed > 3T3_atac_6day_summits_bl_removed.bed
awk '{OFS="\t";} {print $1,$2-99,$3+100,$4,$5}' 3T3_atac_6day_summits_bl_removed.bed > 3T3_atac_6day_summit_200window.bed


#for counts bigWig file

for bam in *rep*_atac_sample_atac.bam
do
    name=$(echo $bam | awk -F"_atac_sample_atac.bam" '{print $1}')
    echo $name
    seqOutBias mm10.fa $bam --skip-bed --no-scale --bw=${name}.bigWig --only-paired --shift-counts --read-size=38
done


#start from here with new meme outputs
#get de novo motifs
mkdir meme_minimal_atac
cd meme_minimal_atac
for i in /Volumes/GUERTIN_2/adipogenesis/atac/twenty_one_clusters/cluster*.meme_output/meme.txt
do
    python /Users/guertinlab/pyscripts/MEME_individual_from_db.py -i $i
done

for i in /Volumes/GUERTIN_2/adipogenesis/atac/meme_minimal_atac/*_cluster*_meme.txt
do
    name=$(echo $i | awk -F"/" '{print $NF}' | awk -F".txt" '{print $1}')
    echo $name
    /Users/guertinlab/Desktop/meme/libexec/meme-5.0.5/ceqlogo -i $i -m 1 > $name.eps
done

cd ..
#i need to organize this better

cat custom.motifs_meme.txt JASPAR2018_CORE_vertebrates_non-redundant_meme.txt uniprobe_mouse_meme.txt > homer_uniprobe_jaspar.txt
mkdir tomtom
cd tomtom
for i in /Volumes/GUERTIN_2/adipogenesis/atac/meme_minimal_atac/*_cluster*_meme.txt
do
    name=$(echo $i | awk -F"/" '{print $NF}' | awk -F".txt" '{print $1}')
    echo $name
    tomtom -no-ssc -o $name.tomtom_output -verbosity 1 -min-overlap 5 -mi 1 -dist ed -evalue -thresh 0.05 $i /Volumes/GUERTIN_2/adipogenesis/atac/homer_uniprobe_jaspar.txt
    tomtom -no-ssc -o $name.tomtom_output -verbosity 1 -min-overlap 5 -mi 1 -dist ed -text -evalue -thresh 0.05 $i /Volumes/GUERTIN_2/adipogenesis/atac/homer_uniprobe_jaspar.txt > $name.tomtom_output/tomtom.txt
done
cd ..

mkdir motifid_clusters
cd motifid_clusters
for file in ../tomtom/*cluster*tomtom_output/*.txt
do
    name=$(echo $file | awk -F"/" '{print $(NF-1)}' | awk -F"_meme.tomtom_output" '{print $1}')
    echo $name
    motifid=$name
    #echo $file
    linenum=$(awk 'END {print NR}' $file)
    #echo $linenum
    first=2
    i=$first
    if [[ $linenum != 5 ]]
    then
	mkdir $name
	while [[  $i -le $linenum ]]
	do
	    head -$i $file | tail -1 > lastline
	    mapid[$i]=$(awk 'END {print $2}' lastline | awk -F"(" '{print $1}')
	    #echo mapid[$i]_${mapid[$i]}
	    echo ${mapid[$i]} >> $name/motifidlist_$motifid.txt
	    ((i = i + 1))
	done
	head -n $(( $(wc -l $name/motifidlist_$motifid.txt | awk '{print $1}') - 4 )) $name/motifidlist_$motifid.txt > $name/motifidlist_$motifid.final.txt
	rm $name/motifidlist_$motifid.txt
	mv $name/motifidlist_$motifid.final.txt $name/motifidlist_$motifid.txt
    fi
done

rm lastline

cd ..
mkdir final_merged_fimo
cat homer_significant_motifs.txt uniprobe_significant_motifs.txt jaspar_significant_motifs.txt > significant_motifs_all_clusters_jaspar_uniprobe_homer.txt


cd final_merged_fimo

#I think that this file should be more inclusive, we should use jaspar, uniprobe, and cisbp.

for file in ../motifid_clusters/*cluster*/motifidlist*.txt
do
    name=$(echo $file | awk -F"/" '{print $NF}' | awk -F".txt" '{print $1}')
    echo $name
    awk 'FNR==NR{a[$1];next}($1 in a){print}' $file ../significant_motifs_all_clusters_jaspar_uniprobe_homer.txt > ${name}_fimo_match.txt
    wordcount=$(wc -l ${name}_fimo_match.txt | awk 'END {print $1}')
    if [[ $wordcount == 0 ]]
    then
	rm ${name}_fimo_match.txt
    fi
done













wget http://meme-suite.org/doc/examples/fimo_example_output_files/fimo.tsv
#if you have 5 non-data rows
sort -n -k 7 fimo.tsv | tail +6 > fimo_sorted.tsv
#want: 1000000
if [[ $(wc -l <fimo_sorted.tsv) -ge 1000 ]]
then
    value=$(head -1000 fimo_sorted.tsv | tail +1000 | awk '{ print $7 }' | bc)
    LC_ALL=C
    awk -v value="$value" '{if (+$7 <= value) print $0}' fimo_sorted.tsv > fimo_sorted_1000_test.tsv
else
    cat fimo_sorted.tsv > fimo_sorted_1000_test.tsv
fi


#next step is TOMTOM against HOMER data base


#make UCSC for PRO-seq
#use seqOutBias and a variant of this script.

for i in G*.bigWig
do
    name=$(echo $i | awk -F"_3T3" '{print $2}')
    #echo $name
    #    cp $1 3T3${name}
    cp $i 3T3${name}
done


for i in 3T3*rep1_pro_plus.bigWig
do
    name=$(echo $i | awk -F"_rep1_pro_plus.bigWig" '{print $1}')
    echo $name
#    plus_1=${name}_rep1_pro_plus.bigWig
#    minus_1=${name}_rep1_pro_minus.bigWig
#    plus_2=${name}_rep2_pro_plus.bigWig
#    minus_2=${name}_rep2_pro_minus.bigWig
#    plus_3=${name}_rep3_pro_plus.bigWig
#    minus_3=${name}_rep3_pro_minus.bigWig
#    Rscript /Users/guertinlab/rscripts/normalization_factor.R $plus_1 $minus_1 $plus_2 $minus_2 $plus_3 $minus_3
    reads=`awk '{SUM+=$2}END{print SUM}' ${name}_normalization.txt`
    norm=$(echo 10000000/$reads | bc -l | xargs printf "%.*f\n" 3)
    #    seqOutBias  /Volumes/GUERTIN_2/adipogenesis/atac/mm10.fa $bam --regions=${reg}.plus.bed --no-scale --bw=PRO_plus_body_0-mer.bigWig --tail-edge --read-size=30
    echo $norm
    bigWigMerge $plus_1  $plus_2 $plus_3  ${name}_plus_merged.bedGraph
    sort -k1,1 -k2,2n ${name}_plus_merged.bedGraph > ${name}_plus_merged_sorted.bedGraph
    python /Users/guertinlab/pyscripts/normalize_bedGraph.py -i ${name}_plus_merged_sorted.bedGraph -s ${norm} -o ${name}_plus_merged_scaled.bedGraph
    bedGraphToBigWig ${name}_plus_merged_scaled.bedGraph /Volumes/GUERTIN_2/adipogenesis/atac/mm10.chrom.sizes ${name}_plus_merged_normalized.bigWig
    bigWigMerge $minus_1  $minus_2 $minus_3  ${name}_minus_merged.bedGraph
    sort -k1,1 -k2,2n ${name}_minus_merged.bedGraph > ${name}_minus_merged_sorted.bedGraph
    python /Users/guertinlab/pyscripts/normalize_bedGraph.py -i ${name}_minus_merged_sorted.bedGraph -s ${norm} -o ${name}_minus_merged_scaled.bedGraph
    bedGraphToBigWig ${name}_minus_merged_scaled.bedGraph /Volumes/GUERTIN_2/adipogenesis/atac/mm10.chrom.sizes ${name}_minus_merged_normalized.bigWig
done


#this can go into DE analysis

cd ../../../pro_dREG
ls 
#GSM3900947_3T3_t0_rep1_pro_minus.bigWig
#GSM3900947_3T3_t0_rep1_pro_plus.bigWig
#GSM3900948_3T3_t0_rep2_pro_minus.bigWig
#GSM3900948_3T3_t0_rep2_pro_plus.bigWig
#GSM3900949_3T3_t0_rep3_pro_minus.bigWig
#GSM3900949_3T3_t0_rep3_pro_plus.bigWig
#GSM3900950_3T3_t20min_rep1_pro_minus.bigWig
#GSM3900950_3T3_t20min_rep1_pro_plus.bigWig
#GSM3900951_3T3_t20min_rep2_pro_minus.bigWig
#GSM3900951_3T3_t20min_rep2_pro_plus.bigWig
#GSM3900952_3T3_t20min_rep3_pro_minus.bigWig
#GSM3900952_3T3_t20min_rep3_pro_plus.bigWig
#GSM3900953_3T3_t2h_rep1_pro_minus.bigWig
#GSM3900953_3T3_t2h_rep1_pro_plus.bigWig
#GSM3900954_3T3_t2h_rep2_pro_minus.bigWig
#GSM3900954_3T3_t2h_rep2_pro_plus.bigWig
#GSM3900955_3T3_t2h_rep3_pro_minus.bigWig
#GSM3900955_3T3_t2h_rep3_pro_plus.bigWig
#GSM3900956_3T3_t3h_rep1_pro_minus.bigWig
#GSM3900956_3T3_t3h_rep1_pro_plus.bigWig
#GSM3900957_3T3_t3h_rep2_pro_minus.bigWig
#GSM3900957_3T3_t3h_rep2_pro_plus.bigWig
#GSM3900958_3T3_t3h_rep3_pro_minus.bigWig
#GSM3900958_3T3_t3h_rep3_pro_plus.bigWig
#GSM3900959_3T3_t40min_rep1_pro_minus.bigWig
#GSM3900959_3T3_t40min_rep1_pro_plus.bigWig
#GSM3900960_3T3_t40min_rep2_pro_minus.bigWig
#GSM3900960_3T3_t40min_rep2_pro_plus.bigWig
#GSM3900961_3T3_t40min_rep3_pro_minus.bigWig
#GSM3900961_3T3_t40min_rep3_pro_plus.bigWig
#GSM3900962_3T3_t4h_rep1_pro_minus.bigWig
#GSM3900962_3T3_t4h_rep1_pro_plus.bigWig
#GSM3900963_3T3_t4h_rep2_pro_minus.bigWig
#GSM3900963_3T3_t4h_rep2_pro_plus.bigWig
#GSM3900964_3T3_t4h_rep3_pro_minus.bigWig
#GSM3900964_3T3_t4h_rep3_pro_plus.bigWig
#GSM3900965_3T3_t60min_rep1_pro_minus.bigWig
#GSM3900965_3T3_t60min_rep1_pro_plus.bigWig
#GSM3900966_3T3_t60min_rep2_pro_minus.bigWig
#GSM3900966_3T3_t60min_rep2_pro_plus.bigWig
#GSM3900967_3T3_t60min_rep3_pro_minus.bigWig
#GSM3900967_3T3_t60min_rep3_pro_plus.bigWig
#GSM3900968_3T3_t6d_rep1_pro_minus.bigWig
#GSM3900968_3T3_t6d_rep1_pro_plus.bigWig
#GSM3900969_3T3_t6d_rep2_pro_minus.bigWig
#GSM3900969_3T3_t6d_rep2_pro_plus.bigWig
#GSM3900970_3T3_t6d_rep3_pro_minus.bigWig
#GSM3900970_3T3_t6d_rep3_pro_plus.bigWig

#ignoring the 6 day time point

plusfiles=$(ls *plus* | grep -v t6d)
bigWigMerge ${plusfiles} preadip_plus_merged.bedGraph
sortBed -i preadip_plus_merged.bedGraph > preadip_plus_merged.sorted.bedGraph
bedGraphToBigWig preadip_plus_merged.sorted.bedGraph ../atac/mm10.chrom.sizes preadip_plus_merged.bigWig

minusfiles=$(ls *minus* | grep -v t6d)
bigWigMerge ${minusfiles} preadip_minus_merged.bedGraph
sortBed -i preadip_minus_merged.bedGraph > preadip_minus_merged.sorted.bedGraph
bedGraphToBigWig preadip_minus_merged.sorted.bedGraph ../atac/mm10.chrom.sizes preadip_minus_merged.bigWig




#dREG proceessing
awk '{OFS="\t";} {print $1,$2,$6}' 3T3_adipogenesis.dREG.peak.full.bed > \
    3T3_adipogenesis.dREG.peak.minus.bed
awk '{OFS="\t";} {print $1,$6,$3}' 3T3_adipogenesis.dREG.peak.full.bed > \
    3T3_adipogenesis.dREG.peak.plus.bed

#go to R




for fq in *rep1_atac_sample_atac.bam
do
    name=$(echo $fq | awk -F"_rep1_atac_sample_atac.bam" '{print $1}')
    echo $name
    repfiles=$(ls ${name}*_atac_sample_atac.bam)
    samtools merge ${name}_merged.rmdup.bam $repfiles
    samtools view -b -f 0x2 ${name}_merged.rmdup.bam -o ${name}_merged.rmdup.1.bam
    rm ${name}_merged.rmdup.bam
    samtools sort -n ${name}_merged.rmdup.1.bam -o ${name}_merged.rmdup.sorted.bam
    rm ${name}_merged.rmdup.1.bam
    samtools fixmate ${name}_merged.rmdup.sorted.bam ${name}_merged.rmdup.fixed.bam
#convert to bed 
    bedtools bamtobed -i ${name}_merged.rmdup.fixed.bam -bedpe > ${name}_merged.rmdup.fixed.bed
#take columns that span the PE1 start and PE2 start
    awk '{OFS="\t";} {print $1,$2,$6,$7,$8,$9}' ${name}_merged.rmdup.fixed.bed > ${name}_merged.rmdup.final.bed
    rm ${name}_merged.rmdup.fixed.bed
    sort -k1,1 -k2,2n ${name}_merged.rmdup.final.bed > ${name}_merged.rmdup.final.sorted.bed
    # from bed to bedgraph
    reads=$(wc -l ${name}_merged.rmdup.final.sorted.bed | awk -F" " '{print $1}')
#three hard coded    
    norm=$(echo 10000000/$reads | bc -l | xargs printf "%.*f\n" 3)
    genomeCoverageBed -bg -scale ${norm} -trackline -trackopts 'track type=bedGraph name='"${name}"' alwaysZero=on visibility=full' -i ${name}_merged.rmdup.final.sorted.bed -g ../mm10.chrom.sizes > ${name}.bedGraph
    bedGraphToBigWig ${name}.bedGraph ../mm10.chrom.sizes ${name}.merged.bigWig
done


#dreg
#Guertin_PRO_early_adipogenesis_c98b87bf-c98d-45c1-8793-ee14d53a873a
#LRT_atac-seq
#TMixClust
#count
