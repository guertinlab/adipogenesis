#!/bin/bash
#SBATCH -n 1
#SBATCH -t 96:00:00
#SBATCH -o motif.families.out
#SBATCH -p standard
#SBATCH -A guertinlab

module load gcc/9.2.0  mvapich2/2.3.3 meme/5.1.0

#define motif families
cd /scratch/bhn9by/ATAC/tomtom_all_query_factors

for meme in *meme.txt
do
    name=$(echo $meme | awk -F".txt_" '{print $NF}' | awk -F"_meme.txt" '{print $1}')
    #echo $name
    tomtom -no-ssc -o $name.tomtom_output -verbosity 1  -incomplete-scores -min-overlap 1 -dist ed -evalue -thresh 0.0005 $meme ../all_query_factors_meme.txt
    cd $name.tomtom_output
    cut -f1,2,5 tomtom.tsv | tail -n +2 | sed '$d' | sed '$d' | sed '$d' | sed '$d' >> ../3_col_combined_motif_db_pre.txt
    cd ..
done

grep -v '#' 3_col_combined_motif_db_pre.txt > 3_col_combined_motif_db.txt
rm 3_col_combined_motif_db_pre.txt
