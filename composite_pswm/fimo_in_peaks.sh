for i in *_composite_fimo_1M.txt
do
    name=$(echo $i | awk -F"/" '{print $NF}' | awk -F"_composite_fimo_1M.txt" '{print $1}')
    echo $name
    intersectBed -loj -a dynamic_peaks.bed -b $i > ${name}_fimo.bed
done



#this type of output is great for input into a seqLogo, but it is a challenge to carry over motif info into the massive data frame.
#for the massive data frame, two columns: presence/absence; 0/score
#
