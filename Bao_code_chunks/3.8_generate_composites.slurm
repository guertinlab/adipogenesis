#!/bin/bash
#SBATCH -n 1
#SBATCH -t 96:00:00
#SBATCH -o motif.composites.out
#SBATCH -p standard
#SBATCH -A guertinlab

#make composite motifs
module load gcc/9.2.0  mvapich2/2.3.3 meme/5.1.0

mkdir composite_motifs
cd composite_motifs

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
	    #increased e-value threshold to 0.05!
	    tomtom -no-ssc -o ${query_factor}.tomtom_output -verbosity 1 -min-overlap 5 -mi 1 -dist pearson -evalue -thresh 0.05 ${query_factor}_meme.txt ref_factors_meme.txt
	    
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
    fi
    
    cd ..
    
done

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
