cd /scratch/abd3x/Adipogenesis/ChIP

i=$SLURM_ARRAY_TASK_ID
dir=`sed "${i}q;d" peak_files_list.txt`
name=`sed "${i}q;d" peak_files_list.txt | awk -F"_macs" '{print $1}'`

echo $name

module purge
module load gcc/9.2.0 bedtools/2.29.2 ucsc-tools/3.7.4

echo 'Converting to .fasta'
slopBed -i ${dir}/${name}_ChIP_summits.bed -g mm10.chrom.sizes -b 50 | intersectBed -v -a stdin -b mm10-blacklist.v2.bed > ${name}_ChIP_peaks.bed
shuf -n 10000 ${name}_ChIP_peaks.bed > ${name}_10K.bed
fastaFromBed -fi /project/genomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa -bed ${name}_10K.bed -fo ${name}_ChIP_peaks.fasta

module purge
#module load gmvapich2/7.1.0_2.3.1 meme/5.1.0
module load gmvapich2/9.2.0_2.3.3 meme/5.1.0

echo 'Starting Meme'

srun meme -p 40 -oc ${name}.meme_output -nmotifs 15 -objfun classic -evt 0.01 -searchsize 0 -minw 6 \
         -maxw 18 -revcomp -dna -markov_order 3 -maxsize 100000000  \
         ${name}_ChIP_peaks.fasta
echo 'DONE'
