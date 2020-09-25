mkdir 6day
mv *6d* 6day/

for bam in *rep1_atac_rmdup.bam
do
    name=$(echo $bam | awk -F"_rep1_atac_rmdup.bam" '{print $1}')
    echo $name
    echo '#SBATCH -o' $name'.peak.calling.out' > temp.txt
    echo 'name='$name > temp2.txt
    cat peak_calling_slurm_header_1.txt temp.txt peak_calling_slurm_header_2.txt temp2.txt peak_calling_slurm_header_3.txt > $name.peak.calling.slurm
    sbatch $name.peak.calling.slurm
    
done

