mkdir individual_memes
cd individual_memes

python ../MEME_individual_from_db.py -i ../homer_uniprobe_jaspar_edited.txt

for file in *meme.txt 
do
    name=$(echo $file | awk -F"homer_uniprobe_jaspar_edited.txt_" '{print $2}')
    mv $file $name
    
done

cd ..

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
	python2.7 ../../../tomtom_output_to_composite.py -i tomtom.xml
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
	Rscript ../../generate_composite_motif.R $family
	cat ../../meme_header.txt ${family}_composite_PSWM.txt > ${family}_meme.txt	
    else
	line=`grep MOTIF *meme.txt`
	cp *meme.txt ${family}_meme.txt
	#sed command written for OS-X, not LINUX
	sed -i '' "s;${line};MOTIF   Composite;g" ${family}_meme.txt
    fi
    ceqlogo -i ${family}_meme.txt -m Composite -o ${family}.eps
    ceqlogo -i ${family}_meme.txt -m Composite -o ${family}.rc.eps -r
    cd ..
    
done

cd ..


fimo --thresh 0.001 --text /Users/guertinlab/Desktop/mike_generate_composites/composite_motifs/PSWM_family_1/PSWM_family_1_meme.txt /Users/guertinlab/genomes/hg38/hg38.fa > AP1_family_composite_fimo.txt
fimo --thresh 0.001 --text /Users/guertinlab/Desktop/mike_generate_composites/composite_motifs/PSWM_family_2/PSWM_family_2_meme.txt /Users/guertinlab/genomes/hg38/hg38.fa > GR_family_composite_fimo.txt
fimo --thresh 0.001 --text /Users/guertinlab/Desktop/mike_generate_composites/composite_motifs/PSWM_family_3/PSWM_family_3_meme.txt /Users/guertinlab/genomes/hg38/hg38.fa > SP_family_composite_fimo.txt
fimo --thresh 0.001 --text /Users/guertinlab/Desktop/mike_generate_composites/composite_motifs/PSWM_family_4/PSWM_family_4_meme.txt /Users/guertinlab/genomes/hg38/hg38.fa > bHLH_family_composite_fimo.txt
fimo --thresh 0.001 --text /Users/guertinlab/Desktop/mike_generate_composites/composite_motifs/PSWM_family_5/PSWM_family_5_meme.txt /Users/guertinlab/genomes/hg38/hg38.fa > TWIST_family_composite_fimo.txt

#this was to get the order of conformity to consensus.
for i in composite_motifs/PSWM_family_*/PSWM_family_*_meme.txt
do
    name=$(echo $i | awk -F"/" '{print $NF}' | awk -F"_meme." '{print $1}')
    echo $i
    tomtom -no-ssc -oc ${name}_ranks.tomtom_output -verbosity 1 -incomplete-scores -min-overlap 0 -dist ed -evalue -thresh 0.05 all_query_factors_meme.txt $i
done
