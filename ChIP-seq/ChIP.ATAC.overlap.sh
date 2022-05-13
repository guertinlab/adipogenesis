rm *.atac.peaks.overlap.bed
for bed in *.round.1.peaks.bed
do
	name=`echo $bed | awk -F".round.1.peaks.bed" '{print $1}'`
	echo $name
	intersectBed -loj -wb -a $bed -b ~/Adipogenesis/ATAC_analysis_redo/all_peaks.bed > $name.atac.peaks.overlap.bed
done



#rm *.motif.count.atac.overlap.bed
#for bed in *.motif.count.bed
#do
#	name=`echo $bed | awk -F".motif.count.bed" '{print $1}'`
#	echo $name
#	intersectBed -loj -wb -a $bed -b ~/Adipogenesis/ATAC_analysis_redo/all_peaks.bed > $name.motif.count.atac.overlap.bed
#done
