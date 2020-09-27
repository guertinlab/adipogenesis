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


#organize de novo motifs by supercluster
echo 'Sort MEME motifs by supercluster'
mkdir supercluster_meme_motifs
cd supercluster_meme_motifs

mkdir up.down
mkdir up.flat
mkdir grad.up
mkdir down.up
mkdir grad.down

#up.flat - 9,23,17,11
#grad.up - 5,8,10,1
#up.down - 6,2,12
#grad.down - 7,3,4,13,14
#down.up - 24

cp -r ../motifid_clusters/cluster6 up.down
cp -r ../motifid_clusters/cluster2 up.down
cp -r ../motifid_clusters/cluster12 up.down

cp -r ../motifid_clusters/cluster9 up.flat
cp -r ../motifid_clusters/cluster23 up.flat
cp -r ../motifid_clusters/cluster17 up.flat
cp -r ../motifid_clusters/cluster11 up.flat

cp -r ../motifid_clusters/cluster5 grad.up
cp -r ../motifid_clusters/cluster8 grad.up
cp -r ../motifid_clusters/cluster10 grad.up
cp -r ../motifid_clusters/cluster1 grad.up

#no motifs tomtom'd out from cluster24 meme result
#cp -r ../motifid_clusters/cluster24 down.up
rm -r down.up

cp -r ../motifid_clusters/cluster7 grad.down
cp -r ../motifid_clusters/cluster3 grad.down
cp -r ../motifid_clusters/cluster4 grad.down
cp -r ../motifid_clusters/cluster13 grad.down
cp -r ../motifid_clusters/cluster14 grad.down

for dir in *.*
do
    echo $dir
    cd $dir
    for int_dir in cluster*
    do
	cd $int_dir
	mv motifidlist* ..
	cd ..
    done
    cat *txt > $dir.denovo.motifs.txt
    mv $dir.denovo.motifs.txt ..
    cd ..
done

cd ..

module load gcc/7.1.0  openmpi/3.1.4 R/4.0.0

#organize fimo motifs by supercluster
echo 'Sort FIMO motifs by supercluster'
mkdir supercluster_fimo_motifs

Rscript sort_fimo_motifs_supercluster.R

#Match the de novo with the fimo identified motifs
echo 'Match MEME and FIMO motifs by supercluster'
mkdir fimo_denovo_match
cd fimo_denovo_match

for file in ../supercluster_meme_motifs/*.txt
do
    name=$(echo $file | awk -F"/" '{print $NF}' | awk -F".denovo" '{print $1}')
    echo $name
    awk 'FNR==NR{a[$1];next}($1 in a){print}' $file ../supercluster_fimo_motifs/$name.fimo.motifs.txt > ${name}_meme_fimo_match.txt
    wordcount=$(wc -l ${name}_meme_fimo_match.txt | awk 'END {print $1}')
    if [[ $wordcount == 0 ]]
    then
	rm ${name}_meme_fimo_match.txt
    fi
done

cat * > all_matched_motifs.txt
cp all_matched_motifs.txt ..

cd ..
