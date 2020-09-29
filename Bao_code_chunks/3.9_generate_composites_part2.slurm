#!/bin/bash
#SBATCH -n 1
#SBATCH -t 96:00:00
#SBATCH -o motif.composites.out
#SBATCH -p standard
#SBATCH -A guertinlab

#generate final composite motif PSWM and corresponding ceqlogo
module load gcc/9.2.0  mvapich2/2.3.3 meme/5.1.0

cd /scratch/bhn9by/ATAC/composite_motifs

for family in PSWM*
do
    cd $family
    num=$(ls *txt | wc -l)

    if [[ $num -ge 2 ]]
    then
      module load gcc/7.1.0  openmpi/3.1.4 R/4.0.0
      Rscript ../../generate_composite_motif.R $family
      cat ../../meme_header.txt ${family}_composite_PSWM.txt > ${family}_meme.txt	
    else
      line=`grep MOTIF *meme.txt`
      cp *meme.txt ${family}_meme.txt
      sed -i "s;${line};MOTIF   Composite;g" ${family}_meme.txt
    fi
      module load gcc/7.1.0 meme/4.10.2
      ceqlogo -i ${family}_meme.txt -m Composite -o ${family}.eps
      ceqlogo -i ${family}_meme.txt -m Composite -o ${family}.rc.eps -r
    cd ..
    
done

cd ..
