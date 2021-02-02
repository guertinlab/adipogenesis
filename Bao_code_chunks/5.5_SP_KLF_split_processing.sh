#module load gcc/9.2.0  mvapich2/2.3.3 meme/5.1.0 bedtools/2.29.2

cd /scratch/bhn9by/ATAC/SP_KLF_split

#bedtools getfasta -fo sp.klf.fasta -fi /project/genomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa -bed /scratch/bhn9by/ATAC/fimo_composites/PSWM_family_8_fimo.bed
module load gcc/7.1.0  openmpi/3.1.4 R/4.0.0

echo Starting prep R script
Rscript ../prep.SP.KLF.fimo.R

module load gcc/9.2.0  mvapich2/2.3.3 meme/5.1.0

echo Starting SP FIMO
fimo --thresh 0.01 --text SP/SP_composite_meme.txt sp_fimo.txt > output_sp1.txt
echo Starting KLF FIMO
fimo --thresh 0.01 --text KLF/KLF_composite_meme.txt sp_fimo.txt > output_klf.txt

module load gcc/7.1.0 bedtools/2.26.0
cp /scratch/bhn9by/ATAC/fimo_composites/PSWM_family_8_2M.bed $PWD/sp_klf_2M.bed
bedtools getfasta -fi /project/genomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa -bed sp_klf_2M.bed > sp_klf_2M.fasta

module load gcc/9.2.0  mvapich2/2.3.3 meme/5.1.0

fimo --thresh 0.01 --text SP/SP_composite_meme.txt sp_klf_2M.txt > output_sp1_2M.txt
fimo --thresh 0.01 --text KLF/KLF_composite_meme.txt sp_klf_2M.txt > output_klf_2M.txt

echo Starting split R script
module load gcc/7.1.0  openmpi/3.1.4 R/4.0.0
Rscript ../SP_KLF_split.R

echo DONE
