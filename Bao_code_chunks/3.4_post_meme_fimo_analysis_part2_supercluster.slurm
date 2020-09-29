#!/bin/bash
#SBATCH -n 1
#SBATCH -t 96:00:00
#SBATCH -o post.meme.fimo.out
#SBATCH -p largemem
#SBATCH -A guertinlab

module load gcc/9.2.0  mvapich2/2.3.3 meme/5.1.0

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
