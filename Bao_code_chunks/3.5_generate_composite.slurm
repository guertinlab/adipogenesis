#!/bin/bash
#SBATCH -n 1
#SBATCH -t 96:00:00
#SBATCH -o motif.clustering.out
#SBATCH -p largemem
#SBATCH -A guertinlab

module load gcc/9.2.0  mvapich2/2.3.3 meme/5.1.0

#extract individual meme files from combined database

mkdir individual_memes
cd individual_memes

python ../MEME_individual_from_db.py -i ../homer_uniprobe_jaspar_edited.txt

for file in *meme.txt 
do
    name=$(echo $file | awk -F"homer_uniprobe_jaspar_edited.txt_" '{print $2}')
    mv $file $name
    
done

cd ..

#tomtom all query factors against all others
mkdir tomtom_all_query_factors
cd tomtom_all_query_factors

cat ../all_matched_motifs.txt | while read factor
do
    echo $factor
    cp ../individual_memes/${factor}_meme.txt $PWD
done

#remove Ptf1a motif b/c it's a forced palindrome
rm Ptf1a_homer_meme.txt
#remove Tal1:Tcf3 motif b/c motif doesn't make biological sense - looks like neither Tal1 nor Tcf3
rm TAL1::TCF3_jaspar_meme.txt
#remove ZNF740 motif
ZNF740_jaspar.bed

cat *meme.txt > ../all_query_factors_meme.txt

#define motif families

for meme in *meme.txt
do
    name=$(echo $meme | awk -F".txt_" '{print $NF}' | awk -F"_meme.txt" '{print $1}')
    #echo $name
    tomtom -no-ssc -o $name.tomtom_output -verbosity 1  -incomplete-scores -min-overlap 0 -dist ed -evalue -thresh 0.0005 $meme ../all_query_factors_meme.txt
    cd $name.tomtom_output
    cut -f1,2,5 tomtom.tsv | tail -n +2 | sed '$d' | sed '$d' | sed '$d' | sed '$d' >> ../3_col_combined_motif_db_pre.txt
    cd ..
done

grep -v '#' 3_col_combined_motif_db_pre.txt > 3_col_combined_motif_db.txt
rm 3_col_combined_motif_db_pre.txt

Rscript visualize_pswms.R

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
