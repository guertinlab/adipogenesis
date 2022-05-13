module purge
module load gmvapich2/9.2.0_2.3.3 meme/5.1.0

i=$SLURM_ARRAY_TASK_ID
name=`sed "${i}q;d" motif.names.txt | awk -F"_macs" '{print $1}'`

fimo --thresh 0.001 --skip-matched-sequence $name.round.1.meme /project/genomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa > $name.round.1.fimo.txt

    #this takes top 1M
score=$(tail -n +2 $name.round.1.fimo.txt | sort -nrk6,6 | awk 'FNR == 1000000 {print $6}')
tail -n +2 $name.round.1.fimo.txt | awk -v sco="$score" '{ if ($6 >= sco) { print } }' | awk '{OFS="\t";} {print $2,$3,$4,$7,$6,$5}' > $name.round.1.fimo.1M.txt

module purge
module load gcc/9.2.0 bedtools/2.29.2

intersectBed -c -a ${name}_ChIP_peaks.bed -b $name.round.1.fimo.1M.txt > $name.round.1.peaks.bed
#intersectBed -wa -a ${name}_ChIP_peaks.bed -b $name.round.1.fimo.1M.txt -u > $name.round.1.peaks.bed
intersectBed -v -a ${name}_ChIP_peaks.bed -b $name.round.1.fimo.1M.txt > $name.round.1.nomotif.peaks.bed

fastaFromBed -fi /project/genomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa -bed \
	     $name.round.1.nomotif.peaks.bed -fo $name.round.1.nomotif.peaks.fasta

module purge
module load gmvapich2/9.2.0_2.3.3 meme/5.1.0

echo 'Starting Meme'

srun meme -p 40 -oc ${name}.round2.meme_output -nmotifs 15 -objfun classic -evt 0.01 -searchsize 0 -minw 6 \
         -maxw 15 -revcomp -dna -markov_order 3 -minsites 2000 -csites 10000 \
         $name.round.1.nomotif.peaks.fasta
