
mkdir 4hour


fdr=0.05
files=$(ls ../*_atac_sample_atac.bam)
cd 4hour
name=3T3_atac
species=mm
ndups=50
macs2 callpeak -t ${files} -f BAMPE -n ${name} --outdir atac4h_${fdr} -g ${species} \
     -B --call-summits --keep-dup ${ndups} -q ${fdr}

cd ../6day

fdr=0.05
files=$(ls *_atac_sample_atac.bam)
name=3T3_atac_6day
species=mm
ndups=50
macs2 callpeak -t ${files} -f BAMPE -n ${name} --outdir atac4h_${fdr} -g ${species} \
     -B --call-summits --keep-dup ${ndups} -q ${fdr}

cd ../
       
# get blacklisted regions
wget http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/mm10-mouse/mm10.blacklist.bed.gz
gunzip mm10.blacklist.bed.gz

# document summit coordinate
cd 4hour/atac4h_0.05/
bedtools subtract -a 3T3_atac_summits.bed -b ../mm10.blacklist.bed > \
     3T3_atac_summits_bl_removed.bed
awk '{OFS="\t";} {print $1,$2-99,$3+100,$4,$5}' 3T3_atac_summits_bl_removed.bed > \
     3T3_atac_summit_200window.bed

cd ../../6day/atac6d_0.05
bedtools subtract -a 3T3_atac_6day_summits.bed -b ../mm10.blacklist.bed > \
     3T3_atac_6day_summits_bl_removed.bed
awk '{OFS="\t";} {print $1,$2-99,$3+100,$4,$5}' 3T3_atac_6day_summits_bl_removed.bed > \
     3T3_atac_6day_summit_200window.bed
