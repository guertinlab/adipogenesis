cd /scratch/abd3x/Adipogenesis/ChIP

i=$SLURM_ARRAY_TASK_ID
sra=`sed "${i}q;d" ChIP_sra_table.txt | awk -F"\t" '{print $1}'`
name=`sed "${i}q;d" ChIP_sra_table.txt | awk -F"\t" '{print $2}'`

echo $name

module purge
module load sratoolkit/2.10.5

echo 'Downloading .fastq'
fasterq-dump -3 $sra -o $name.fastq

module purge
module load gcc/7.1.0 bowtie2/2.2.9

echo 'Starting alignment'
bowtie2 -x /project/genomes/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome -U $name.fastq -S $name.sam
