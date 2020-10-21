#!/bin/bash

cd /scratch/bhn9by/ATAC/

#split SP/KLF
rm sp.txt
rm -r SP_composite

#Change family number to match SP_KLF
grep Sp PSWM_family_4.txt > sp.txt
grep SP PSWM_family_4.txt >> sp.txt

#generate composite for SP only
mkdir SP_composite
cd SP_composite

cat ../sp.txt | while read line
do
    echo $line
    cp ../../individual_memes/${line}_meme.txt $PWD
done

#query tomtom for each factor against all others
module load gcc/9.2.0  mvapich2/2.3.3 meme/5.1.0
cat ../sp.txt | { while read line
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

#generate composite PSWM
module load gcc/7.1.0  openmpi/3.1.4 R/4.0.0
Rscript ../generate_composite_motif.R SP
cat ../meme_header.txt SP_composite_PSWM.txt > SP_meme.txt

#generate logo
module load gcc/7.1.0 meme/4.10.2
ceqlogo -i SP_meme.txt -m Composite -o SP.eps
ceqlogo -i SP_meme.txt -m Composite -o SP.rc.eps -r

cd ..
