for i in *cluster*bed
do
    name=$(echo $i | awk -F"_" '{print $NF}' | awk -F".bed" '{print $1}')
    echo $name
    echo '#SBATCH -o' $name'.meme.out' > temp.txt
    echo 'i='$i > temp2.txt
    cat meme_slurm_header_1.txt temp.txt meme_slurm_header_2.txt temp2.txt meme_slurm_header_3.txt > $name.meme.slurm
    sbatch $name.meme.slurm                                                                                                                                                           
    rm temp.txt
    rm temp2.txt
done
