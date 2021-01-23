#!/bin/bash

for family in PSWM*
do
    cd $family
    num=$(ls *txt | wc -l)

    if [[ $num -ge 2 ]]
    then
    
    #generate composite PSWM
    module load gcc/7.1.0  openmpi/3.1.4 R/4.0.0
    Rscript ../../generate_composite_motif.R $family
    cat ../../meme_header.txt ${family}_composite_PSWM.txt > ${family}_meme.txt	
    else
    line=`grep MOTIF *meme.txt`
    cp *meme.txt ${family}_meme.txt
    sed -i "s;${line};MOTIF   Composite;g" ${family}_meme.txt
    fi
    
    #generate logo
    module load gcc/7.1.0 meme/4.10.2
    ceqlogo -i ${family}_meme.txt -m Composite -o ${family}.eps
    ceqlogo -i ${family}_meme.txt -m Composite -o ${family}.rc.eps -r
    cd ..
    
done

cd ..
module purge
