#!/bin/bash
#SBATCH -n 1
#SBATCH -t 96:00:00
#SBATCH -o motif.clustering.spklf.out
#SBATCH -p largemem
#SBATCH -A guertinlab

module load gcc/9.2.0  mvapich2/2.3.3 meme/5.1.0
cd /scratch/bhn9by/ATAC/

#split SP/KLF
rm klf.txt
rm sp.txt
rm -r SP_composite
rm -r KLF_composite

#Change family number to match SP_KLF
grep Klf PSWM_family_3.txt > klf.txt
grep KLF PSWM_family_3.txt >> klf.txt

grep Sp PSWM_family_3.txt > sp.txt
grep SP PSWM_family_3.txt >> sp.txt

mkdir SP_composite
cd SP_composite

cat ../sp.txt | while read line
do
    echo $line
    cp ../../individual_memes/${line}_meme.txt $PWD
done

echo ''
    
cat ../sp.txt | { while read line
	do
	    query_factor=$line
	    rm ref_factors_meme.txt
	    mv ${query_factor}_meme.txt ..
	    cat *_meme.txt > ref_factors_meme.txt
	    mv ../${query_factor}_meme.txt $PWD
	    tomtom -no-ssc -o ${query_factor}.tomtom_output -verbosity 1 -incomplete-scores -min-overlap 0 -dist ed -evalue -thresh 0.0005 ${query_factor}_meme.txt ref_factors_meme.txt
	    
	    if [[ $(wc -l < ${query_factor}.tomtom_output/tomtom.tsv) -ge $max_motif ]]
	    then
		max_motif=$(wc -l < ${query_factor}.tomtom_output/tomtom.tsv)
		final_query=$query_factor
	    fi
	    
	done
	echo FINAL_QUERY IS $final_query
	wc -l ${final_query}.tomtom_output/tomtom.tsv
	cd ${final_query}.tomtom_output
	python ../../../tomtom_output_to_composite.py -i tomtom.xml
	mv tomtom.xml_test_index_pswm.txt ../composite.values.txt
	mv tomtom.xml_test_index_rc_offset.txt ../composite.index.txt
	cd ../..
}

Rscript ../generate_composite_motif.R SP
cat ../meme_header.txt SP_composite_PSWM.txt > SP_meme.txt	
ceqlogo -i SP_meme.txt -m Composite -o SP.eps
ceqlogo -i SP_meme.txt -m Composite -o SP.rc.eps -r

cd ..

mkdir KLF_composite
cd KLF_composite

cat ../klf.txt | while read line
do
    echo $line
    cp ../../individual_memes/${line}_meme.txt $PWD
done

echo ''
    
cat ../klf.txt | { while read line
	do
	    query_factor=$line
	    rm ref_factors_meme.txt
	    mv ${query_factor}_meme.txt ..
	    cat *_meme.txt > ref_factors_meme.txt
	    mv ../${query_factor}_meme.txt $PWD
	    tomtom -no-ssc -o ${query_factor}.tomtom_output -verbosity 1 -incomplete-scores -min-overlap 0 -dist ed -evalue -thresh 0.0005 ${query_factor}_meme.txt ref_factors_meme.txt
	    
	    if [[ $(wc -l < ${query_factor}.tomtom_output/tomtom.tsv) -ge $max_motif ]]
	    then
		max_motif=$(wc -l < ${query_factor}.tomtom_output/tomtom.tsv)
		final_query=$query_factor
	    fi
	    
	done
	echo FINAL_QUERY IS $final_query
	wc -l ${final_query}.tomtom_output/tomtom.tsv
	cd ${final_query}.tomtom_output
	python ../../../tomtom_output_to_composite.py -i tomtom.xml
	mv tomtom.xml_test_index_pswm.txt ../composite.values.txt
	mv tomtom.xml_test_index_rc_offset.txt ../composite.index.txt
	cd ../..
}

Rscript ../generate_composite_motif.R KLF
cat ../meme_header.txt KLF_composite_PSWM.txt > KLF_meme.txt	
ceqlogo -i KLF_meme.txt -m Composite -o KLF.eps
ceqlogo -i KLF_meme.txt -m Composite -o KLF.rc.eps -r

cd ..
