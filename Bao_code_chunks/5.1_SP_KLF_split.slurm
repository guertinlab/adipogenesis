#!/bin/bash
#SBATCH -n 1
#SBATCH -t 96:00:00
#SBATCH -o sp.klf.split.out
#SBATCH -p standard
#SBATCH -A guertinlab

cd /scratch/bhn9by/ATAC/

rm -r SP_KLF_split
mkdir SP_KLF_split
cd SP_KLF_split

echo Generate composite SP and KLF motifs
mkdir SP
cd SP

module load gcc/9.2.0  mvapich2/2.3.3 meme/5.1.0

motifs=`grep -E 'SP|Sp' ../../PSWM_family_7.txt`

for motif in $motifs
do
    cp  /scratch/bhn9by/ATAC/individual_memes/$motif*meme.txt $PWD
done

max_motif=0
final_query=''

for line in *meme.txt
do
    name=$(echo $line | awk -F"_meme.txt" '{print $1}')
    echo $name
    query_factor=$line
    mv ${query_factor} ..
    rm ref_factors.txt
    cat *_meme.txt > ref_factors.txt
    mv ../${query_factor} $PWD
    #increased e-value threshold to 0.05!
    tomtom -no-ssc -oc $name.tomtom_output -verbosity 1 -min-overlap 5 -mi 1 -dist pearson -evalue -thresh 0.05 ${query_factor} ref_factors.txt
    if [[ $(wc -l < $name.tomtom_output/tomtom.tsv) -ge $max_motif ]]
    then
	max_motif=$(wc -l < $name.tomtom_output/tomtom.tsv)
	final_query=$name
    fi
done

echo FINAL_QUERY IS $final_query
wc -l $final_query.tomtom_output/tomtom.tsv
cd $final_query.tomtom_output
python ../../../tomtom_output_to_composite.py -i tomtom.xml
mv tomtom.xml_test_index_pswm.txt ../composite.values.txt
mv tomtom.xml_test_index_rc_offset.txt ../composite.index.txt
cd ..

module load gcc/7.1.0  openmpi/3.1.4 R/4.0.0

Rscript ../../generate_composite_motif.SP.R

cat ../../meme_header.txt SP_composite_PSWM.txt > SP_composite_meme.txt	

module load gcc/7.1.0 meme/4.10.2

ceqlogo -i SP_composite_meme.txt -m Composite -o SP_composite.eps
ceqlogo -i SP_composite_meme.txt -m Composite -o SP_composite.rc.eps -r

cd ..

mkdir KLF
cd KLF

module load gcc/9.2.0  mvapich2/2.3.3 meme/5.1.0

motifs=`grep -E 'KLF|Klf' ../../PSWM_family_7.txt`

for motif in $motifs
do
    cp  /scratch/bhn9by/ATAC/individual_memes/$motif*meme.txt $PWD
done

max_motif=0
final_query=''

for line in *meme.txt
do
    name=$(echo $line | awk -F"_meme.txt" '{print $1}')
    echo $name
    query_factor=$line
    mv ${query_factor} ..
    rm ref_factors.txt
    cat *_meme.txt > ref_factors.txt
    mv ../${query_factor} $PWD
    #increased e-value threshold to 0.05!
    tomtom -no-ssc -oc $name.tomtom_output -verbosity 1 -min-overlap 5 -mi 1 -dist pearson -evalue -thresh 0.05 ${query_factor} ref_factors.txt
    if [[ $(wc -l < $name.tomtom_output/tomtom.tsv) -ge $max_motif ]]
    then
	max_motif=$(wc -l < $name.tomtom_output/tomtom.tsv)
	final_query=$name
    fi
done

echo FINAL_QUERY IS $final_query
wc -l $final_query.tomtom_output/tomtom.tsv
cd $final_query.tomtom_output
python ../../../tomtom_output_to_composite.py -i tomtom.xml
mv tomtom.xml_test_index_pswm.txt ../composite.values.txt
mv tomtom.xml_test_index_rc_offset.txt ../composite.index.txt
cd ..

module load gcc/7.1.0  openmpi/3.1.4 R/4.0.0

Rscript ../../generate_composite_motif.KLF.R

cat ../../meme_header.txt KLF_composite_PSWM.txt > KLF_composite_meme.txt	

module load gcc/7.1.0 meme/4.10.2

ceqlogo -i KLF_composite_meme.txt -m Composite -o KLF_composite.eps
ceqlogo -i KLF_composite_meme.txt -m Composite -o KLF_composite.rc.eps -r

cd ..
