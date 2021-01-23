#!/bin/bash

cd /scratch/bhn9by/ATAC

wget https://raw.githubusercontent.com/guertinlab/adipogenesis/master/Bao_code_chunks/misc_scripts/tomtom_output_to_composite.py

mkdir composite_motifs
cd composite_motifs

#query tomtom for each factor against all others
module purge
module load gcc/9.2.0  mvapich2/2.3.3 meme/5.1.0
for txt in ../PSWM_family*.txt
do

    dir_name=$(echo $txt | awk -F'../' '{print $2}' | awk -F'.txt' '{print $1}')
    echo $dir_name
    mkdir $dir_name
    cd $dir_name
    
    cat ../$txt | while read line
    do
	echo $line
	cp ../../individual_memes/${line}_meme.txt $PWD
    done
    echo ''
    
    query_factor=`head -1 ../$txt`
    
    if [[ $(wc -l < ../$txt) -ge 2 ]]
    then   	
	cat ../$txt | { while read line
	do
	    query_factor=$line
	    rm ref_factors_meme.txt
	    mv ${query_factor}_meme.txt ..
	    cat *_meme.txt > ref_factors_meme.txt
	    mv ../${query_factor}_meme.txt $PWD
	    tomtom -no-ssc -o ${query_factor}.tomtom_output -verbosity 1 -incomplete-scores -min-overlap 1 -dist ed -evalue -thresh 0.0005 ${query_factor}_meme.txt ref_factors_meme.txt
	    
	    if [[ $(wc -l < ${query_factor}.tomtom_output/tomtom.tsv) -ge $max_motif ]]
	    then
		max_motif=$(wc -l < ${query_factor}.tomtom_output/tomtom.tsv)
		final_query=$query_factor
	    fi
	    
	done
	echo FINAL_QUERY IS $final_query
	wc -l ${final_query}.tomtom_output/tomtom.tsv
	cd ${final_query}.tomtom_output
	python2.7 ../../../tomtom_output_to_composite.py -i tomtom.xml
	mv tomtom.xml_test_index_pswm.txt ../composite.values.txt
	mv tomtom.xml_test_index_rc_offset.txt ../composite.index.txt
	cd ../..
	}
    fi
    
    cd ..
    
done
module purge
