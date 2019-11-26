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

https://data.cyverse.org/dav-anon/iplant/home/guertin/preadip.bigWig
https://genome.ucsc.edu/s/Mike%20Guertin/mm10_0%2D4hr_merged


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
    bedGraphToBigWig ${name}.bedGraph mm10.chrom.sizes ${name}.merged.bigWig
done

#6 day
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


# get blacklisted regions
cd /Volumes/GUERTIN_2/adipogenesis/atac
wget http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/mm10-mouse/mm10.blacklist.bed.gz
gunzip mm10.blacklist.bed.gz

# document summit coordinate
cd 4hour/atac4h_0.05/
bedtools subtract -a 3T3_atac_summits.bed -b ../../mm10.blacklist.bed > 3T3_atac_summits_bl_removed.bed
awk '{OFS="\t";} {print $1,$2-99,$3+100,$4,$5}' 3T3_atac_summits_bl_removed.bed > 3T3_atac_summit_200window.bed

cd ../../6day/atac6d_0.05
bedtools subtract -a 3T3_atac_6day_summits.bed -b ../../mm10.blacklist.bed > 3T3_atac_6day_summits_bl_removed.bed
awk '{OFS="\t";} {print $1,$2-99,$3+100,$4,$5}' 3T3_atac_6day_summits_bl_removed.bed > 3T3_atac_6day_summit_200window.bed

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

#dreg
#Guertin_PRO_early_adipogenesis_c98b87bf-c98d-45c1-8793-ee14d53a873a
#LRT_atac-seq
#TMixClust
#count


#for counts bigWig file

for bam in *rep*_atac_sample_atac.bam
do
    name=$(echo $bam | awk -F"_atac_sample_atac.bam" '{print $1}')
    echo $name
    seqOutBias mm10.fa $bam --skip-bed --no-scale --bw=${name}.bigWig --only-paired --shift-counts --read-size=38
done
