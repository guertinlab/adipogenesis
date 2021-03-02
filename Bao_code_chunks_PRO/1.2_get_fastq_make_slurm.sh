#On Rivanna
mkdir /scratch/bhn9by/PRO
cd /scratch/bhn9by/PRO

#Upload sra.metadata.csv and SRR_Acc_List.txt to PRO directory

#Set up conda env to get access to sratoolkit
module load bioconda
conda create -n myenv2 -c bioconda ucsc-bedgraphtobigwig ucsc-bigwigmerge

#get mm10 tallymer files for seqOutBias
cp ../ATAC/genome* $PWD

#download .fastq from SRA
#make a unique slurm file for each replicate and run them in parallel
cat SRR_Acc_List.txt | while read acc
do
    echo $acc
    echo '#SBATCH -o' $acc'.out' > temp.txt
    echo 'fasterq-dump' $acc > temp2.txt
    echo 'gzip' $acc'.fastq' > temp3.txt
    cat sra_slurm_header_1.txt temp.txt sra_slurm_header_2.txt temp2.txt temp3.txt > $acc.slurm
    sbatch $acc.slurm
    rm temp.txt
    rm temp2.txt
    rm temp3.txt
done

#Alternatively, you can run fasterq-dump and gzip in series. Approximately 7-10 minutes per accession.
cat SRR_Acc_List.txt | while read acc
do
    echo $acc
    fasterq-dump $acc
    gzip $acc.fastq
done


#After all jobs are done, rename files to actual sample names
for fq in SRR*.fastq.gz
do
    name=$(echo $fq | awk -F".fastq.gz" '{print $1}')
    echo $name
    line=$(grep $name sra.metadata.csv)
    treat=$(echo $line | awk -F',' '{print $31}')
    echo $treat
    mv $fq $treat
done
