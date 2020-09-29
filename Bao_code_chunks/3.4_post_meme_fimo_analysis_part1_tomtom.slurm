#!/bin/bash
#SBATCH -n 1
#SBATCH -t 96:00:00
#SBATCH -o post.meme.fimo.out
#SBATCH -p largemem
#SBATCH -A guertinlab

module load gcc/9.2.0  mvapich2/2.3.3 meme/5.1.0

#pulls each de novo motif from all the clusters and saves each as its own file
wget https://github.com/guertinlab/adipogenesis/blob/2fdf0bbab4fe6f368f5a60e42f7899b6570ff71c/motif_analysis/MEME_individual_from_db.py

#run tomtom on meme output
mkdir tomtom
cd tomtom

echo 'Running TOMTOM'
for i in ../meme_motif_enrichment/*output
do
    name=$(echo $i | awk -F"/" '{print $NF}' | awk -F".meme" '{print $1}')
    echo $name
    tomtom -no-ssc -o $name.tomtom_output -verbosity 1 -min-overlap 5 -mi 1 -dist ed -evalue -thresh 0.05 $i/meme.txt ../homer_uniprobe_jaspar_edited.txt
    tomtom -no-ssc -o $name.tomtom_output -verbosity 1 -min-overlap 5 -mi 1 -dist ed -text -evalue -thresh 0.05 $i/meme.txt ../homer_uniprobe_jaspar_edited.txt > $name.tomtom_output/tomtom.txt
done

cd ..

mkdir motifid_clusters
cd motifid_clusters

#extract motif names from tomtom output
echo 'Extracting Motif Names'
for file in ../tomtom/*cluster*tomtom_output/*.txt
do
    name=$(echo $file | awk -F"/" '{print $(NF-1)}' | awk -F".tomtom_output" '{print $1}')
    echo $name
    motifid=$name
    #echo $file
    linenum=$(awk 'END {print NR}' $file)
    #echo $linenum
    first=2
    i=$first
    if [[ $linenum != 5 ]]
    then
	mkdir $name
	while [[  $i -le $linenum ]]
	do
	    head -$i $file | tail -1 > lastline
	    mapid[$i]=$(awk 'END {print $2}' lastline | awk -F"(" '{print $1}')
	    #echo mapid[$i]_${mapid[$i]}
	    echo ${mapid[$i]} >> $name/motifidlist_$motifid.txt
	    ((i = i + 1))
	done
	head -n $(( $(wc -l $name/motifidlist_$motifid.txt | awk '{print $1}') - 4 )) $name/motifidlist_$motifid.txt > $name/motifidlist_$motifid.final.txt
	rm $name/motifidlist_$motifid.txt
	mv $name/motifidlist_$motifid.final.txt $name/motifidlist_$motifid.txt
    fi
done

rm lastline
cd ..
