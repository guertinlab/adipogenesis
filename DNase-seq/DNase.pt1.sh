cd /scratch/abd3x/Adipogenesis/DNase

i=$SLURM_ARRAY_TASK_ID
sra=`sed "${i}q;d" DNase_acc_list.txt | awk -F"\t" '{print $1}'`
name=`sed "${i}q;d" DNase_acc_list.txt | awk -F"\t" '{print $2}'`

echo $name

module purge
module load sratoolkit/2.10.5

echo 'Downloading .fastq'
fasterq-dump -3 $sra -o $name.fastq

module purge
module load gcc/7.1.0 bowtie2/2.2.9 samtools/1.10

echo 'Aligning'
bowtie2 -x /project/genomes/Homo_sapiens/UCSC/hg38/Sequence/Bowtie2Index/genome -U $name.fastq -S $name.sam

echo 'Removing duplicates'
samtools view -b -q 10 ${name}.sam | samtools sort -n - | samtools fixmate -m - - | samtools sort - | samtools markdup -r - ${name}_rmdup.bam
rm ${name}.sam
gzip ${name}*fastq

echo 'DONE'
