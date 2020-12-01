#!/bin/bash
#SBATCH -n 1
#SBATCH -t 96:00:00
#SBATCH -o fimo.enrichment.out
#SBATCH -p standard
#SBATCH -A guertinlab

module load gcc/7.1.0 bedtools/2.26.0 openmpi/3.1.4 R/4.0.0

#change the directory to wherever you are keeping the ATAC cluster .bed files
cd /scratch/bhn9by/ATAC

tab=$'\t'

mkdir fimo_motif_enrichment
cd fimo_motif_enrichment

#loop through each motif
for motif in /scratch/bhn9by/ATAC/Top_motif/*bed
do
    motif_name=$(echo $motif | awk -F"/" '{print $NF}' | awk -F".bed" '{print $1}')
    echo $motif_name
    touch ${motif_name}.txt
    echo "name""$tab""with.motif""$tab""without.motif" >> ${motif_name}.txt
    
    #find how many nondynamic peaks contain the motif 
    nondyn_with_motif=$(intersectBed -wa -a ../nondynamic_peaks.bed -b $motif | sort -u | wc -l)
    nondyn_without_motif=$(intersectBed -v -a ../nondynamic_peaks.bed -b $motif | sort -u | wc -l)

    echo "nondynamic$tab$nondyn_with_motif$tab$nondyn_without_motif" >> ${motif_name}.txt
    
    #loop though each cluster
    for bed in ../cluster*bed
    do
	name=$(echo $bed | awk -F"cluster_bed_" '{print $2}' | awk -F".bed" '{print $1}')
	#find how many peaks in each cluster contain the motif

	#calculate percentage of cluster peaks that have the motif
	cluster_with_motif=$(intersectBed -wa -a $bed -b $motif | sort -u | wc -l)
	cluster_without_motif=$(intersectBed -v -a $bed -b $motif | sort -u | wc -l)
	
	echo "$name$tab$cluster_with_motif$tab$cluster_without_motif" >> ${motif_name}.txt

    done

    Rscript ../fimo_motif_enrichment.R ${motif_name}.txt

done
cd ..
