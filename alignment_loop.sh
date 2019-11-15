#this is complicated in order to deal with Drosophila reads. However, since they Drosophila reads that erroneously align to human are so minimal, this complicated alignemtn schme is not necessary. Note that drosophila reads were spiked in to experiment, but not used doe to low abundance.

#I need to download data directly from SRA/GEO to see if 20min rep2 and 60 min rep 2 are switched

for i in *_PE1.fastq.gz
do
	name=$(echo $i | awk -F"/" '{print $NF}' | awk -F"_atac_PE1.fastq.gz" '{print $1}')
	echo $name
	mkdir $name
	cp ${name}_atac_PE1.fastq.gz $name
	cp ${name}_atac_PE2.fastq.gz $name
	cd $name
	gunzip *.gz
	fastq_pair ${name}_atac_PE1.fastq ${name}_atac_PE2.fastq
	rm *single.fq

	printf "\ninitial alignment to sample species \n"
	bowtie2 --maxins 500 -x /Volumes/GUERTIN_2/adipogenesis/atac/mm10 -1 ${name}_atac_PE1.fastq.paired.fq -2 ${name}_atac_PE2.fastq.paired.fq -S ${name}_atac_smp.sam

	# sort the sam file, filter based on MAPQ=10, remove duplicates, and generate a bam file
	samtools view -b -q 10 ${name}_atac_smp.sam | samtools sort -n - | samtools fixmate -m - - | samtools sort - | samtools markdup -r - ${name}_atac_smp_rmdup.bam
	printf "\nnumber of sample MAPQ-filtered and deduplicated reads \n"
	samtools view -c ${name}_atac_smp_rmdup.bam
	rm ${name}_atac_smp.sam
	
	# convert the bam file into a fastq file for re-alignment to the 'spike-in' species
	# spike-in reads will be removed from the generated fastq files
	samtools sort -n ${name}_atac_smp_rmdup.bam -o ${name}_atac_smp_sort.bam
	bedtools bamtofastq -i ${name}_atac_smp_sort.bam -fq ${name}_atac_smp_r1.fastq -fq2 ${name}_atac_smp_r2.fastq
	rm ${name}_atac_smp_rmdup.bam ${name}_atac_smp_sort.bam

	# align processed read data to 'spike-in' species
	# here we are aligning de-duplicated and 'uniquely' mapped reads from alignment to the sample species
	# the reads aligned here are the 'intersect' reads  
	# these reads will be removed from the sample and re-aligned with no mismatches to seed
	fastq_pair ${name}_atac_smp_r1.fastq ${name}_atac_smp_r2.fastq
	printf "\nalign sample-aligned and MAPQ-filtered data to 'spike-in' species \n"
	bowtie2 --maxins 500 -x /Volumes/GUERTIN_2/adipogenesis/atac/dm6 -1 ${name}_atac_smp_r1.fastq.paired.fq -2 ${name}_atac_smp_r2.fastq.paired.fq -S ${name}_atac_spk.sam
	samtools view -b -q 10 ${name}_atac_spk.sam | samtools sort -n - | samtools fixmate -m - - | samtools sort - | samtools markdup -r - ${name}_atac_spk_rmdup.bam
	printf "\nnumber of sample MAPQ-filtered and deduplicated reads \n"
	samtools view -c ${name}_atac_spk_rmdup.bam
	samtools sort -n ${name}_atac_spk_rmdup.bam -o ${name}_atac_spk_sort.bam
	bedtools bamtofastq -i ${name}_atac_spk_sort.bam -fq ${name}_atac_spk_r1.fastq -fq2 ${name}_atac_spk_r2.fastq
	rm ${name}_atac_spk.sam ${name}_atac_spk_rmdup.bam ${name}_atac_spk_sort.bam *single*

	
	# remove doubly unique (intersect) reads from sample
	# input 1 to python script - reads will be removed from this file
	# input 2 to python script - these reads will be removed
	# input 3 to python script - output file
	fastq_pair ${name}_atac_spk_r1.fastq ${name}_atac_spk_r2.fastq
	printf "\nremove intersect from sample \n"
	python removeOverlap2.py ${name}_atac_smp_r1.fastq.paired.fq ${name}_atac_spk_r1.fastq.paired.fq ${name}_atac_smp_r1_spkrem.fastq
	python removeOverlap2.py ${name}_atac_smp_r2.fastq.paired.fq ${name}_atac_spk_r2.fastq.paired.fq ${name}_atac_smp_r2_spkrem.fastq
	printf "\nline count after removing overlap reads \n"
	python pyShell.py ${name}_atac_smp_r1_spkrem.fastq
	python pyShell.py ${name}_atac_smp_r2_spkrem.fastq
	rm ${name}_atac_smp_r1.fastq ${name}_atac_smp_r2.fastq *single*

	# realign thes doubly unique intersect reads to the sample species 
	# do not allow mismatches
	printf "\nalign these doubly overlap to the sample species without mismatches \n" 
	bowtie2 --maxins 500 --no-1mm-upfront -x /Volumes/GUERTIN_2/adipogenesis/atac/mm10 -1 ${name}_atac_spk_r1.fastq.paired.fq -2 ${name}_atac_spk_r2.fastq.paired.fq -S ${name}_atac_smp_spkrem_smp.sam

	# realign these doubly unique intersect reads to the spike-in species without mismatches
	# uniquely aligned reads should be removed from the sample reads 
	printf "\nalign these doubly overlap to the spike-in species without mismatches \n"
	bowtie2 --maxins 500 --no-1mm-upfront -x /Volumes/GUERTIN_2/adipogenesis/atac/dm6 -1 ${name}_atac_spk_r1.fastq.paired.fq -2 ${name}_atac_spk_r2.fastq.paired.fq -S ${name}_atac_smp_spkrem_spk.sam

	# extract unique alignments in the sample no mismatch file
	samtools view -b -q 10 ${name}_atac_smp_spkrem_smp.sam | samtools sort -n - | samtools fixmate -m - - | samtools sort - | samtools markdup -r - ${name}_atac_smp_spkrem_smp_rmdup.bam
	printf "\nnumber of sample no-mm reads \n"
	samtools view -c ${name}_atac_smp_spkrem_smp_rmdup.bam
	samtools sort -n ${name}_atac_smp_spkrem_smp_rmdup.bam -o ${name}_atac_smp_spkrem_smp_sort.bam
	bedtools bamtofastq -i ${name}_atac_smp_spkrem_smp_sort.bam -fq ${name}_atac_smp_spkrem_smp_r1.fastq -fq2 ${name}_atac_smp_spkrem_smp_r2.fastq
	rm ${name}_atac_smp_spkrem_smp.sam ${name}_atac_smp_spkrem_smp_rmdup.bam
	rm ${name}_atac_smp_spkrem_smp_sort.bam
	rm ${name}_atac_spk_r1.fastq ${name}_atac_spk_r2.fastq
	
	# extract unique alignments in the spike-in no mismatch file
	samtools view -b -q 10 ${name}_atac_smp_spkrem_spk.sam | samtools sort -n - | samtools fixmate -m - - | samtools sort - | samtools markdup -r - ${name}_atac_smp_spkrem_spk_rmdup.bam
	printf "\nnumber of spike-in no-mm reads \n"
	samtools view -c ${name}_atac_smp_spkrem_spk_rmdup.bam
	samtools sort -n ${name}_atac_smp_spkrem_spk_rmdup.bam -o ${name}_atac_smp_spkrem_spk_sort.bam
	bedtools bamtofastq -i ${name}_atac_smp_spkrem_spk_sort.bam -fq ${name}_atac_smp_spkrem_spk_r1.fastq -fq2 ${name}_atac_smp_spkrem_spk_r2.fastq
	rm ${name}_atac_smp_spkrem_spk.sam ${name}_atac_smp_spkrem_spk_rmdup.bam
	rm ${name}_atac_smp_spkrem_spk_sort.bam
	
	# remove reads from the no-mm sample that align without mismatches to the spike-in
	# input 1 to python script - reads will be removed from this file
	# input 2 to python script - these reads will be removed 
	# input 3 to python script - output file
	fastq_pair ${name}_atac_smp_spkrem_spk_r1.fastq ${name}_atac_smp_spkrem_spk_r2.fastq
	fastq_pair ${name}_atac_smp_spkrem_smp_r1.fastq ${name}_atac_smp_spkrem_smp_r2.fastq
	printf "\nremove spike in no-mm from sample no-mm reads to get reads that should go back to the sample \n"
	python backToSample2.py ${name}_atac_smp_spkrem_smp_r1.fastq.paired.fq ${name}_atac_smp_spkrem_spk_r1.fastq.paired.fq ${name}_atac_smp_spkrem2_r1.fastq
	python backToSample2.py ${name}_atac_smp_spkrem_smp_r2.fastq.paired.fq ${name}_atac_smp_spkrem_spk_r2.fastq.paired.fq ${name}_atac_smp_spkrem2_r2.fastq
	rm *single*
	
	# add the sample-specific no-mm reads back to the sample
	cat ${name}_atac_smp_r1_spkrem.fastq ${name}_atac_smp_spkrem2_r1.fastq > ${name}_atac_smp_r1.fastq
	cat ${name}_atac_smp_r2_spkrem.fastq ${name}_atac_smp_spkrem2_r2.fastq > ${name}_atac_smp_r2.fastq
	rm ${name}_atac_smp_r1_spkrem.fastq ${name}_atac_smp_spkrem2_r1.fastq
	rm ${name}_atac_smp_r2_spkrem.fastq ${name}_atac_smp_spkrem2_r2.fastq
	
	# re-align final sample reads and generate the output bam file
	fastq_pair ${name}_atac_smp_r1.fastq ${name}_atac_smp_r2.fastq
	printf "\nfinal alignment for sample data \n"
	bowtie2 --maxins 500 -x /Volumes/GUERTIN_2/adipogenesis/atac/mm10 -1 ${name}_atac_smp_r1.fastq.paired.fq -2 ${name}_atac_smp_r2.fastq.paired.fq -S ${name}_atac_sample1.sam
	samtools view -b -q 10 ${name}_atac_sample1.sam | samtools sort -n - | samtools fixmate -m - - | samtools sort - | samtools markdup -r - ${name}_atac_sample_atac.bam
	printf "\nfinal number of sample reads \n"
	samtools view -c ${name}_atac_sample_atac.bam
	rm ${name}_atac_sample1.sam *single*
	
	
	####################################################################
	## isolate spike-in-specific reads 
	####################################################################
	
	# align to the spike genome
	printf "\ninitial alignment to spike-in species \n"
	bowtie2 --maxins 500 -x /Volumes/GUERTIN_2/adipogenesis/atac/dm6 -1 ${name}_atac_PE1.fastq.paired.fq -2 ${name}_atac_PE2.fastq.paired.fq -S ${name}_atac_spk.sam
	
	# sort the sam file, filter based on MAPQ=10, remove duplicates, and generate a bam file
	samtools view -b -q 10 ${name}_atac_spk.sam | samtools sort -n - | samtools fixmate -m - - | samtools sort - | samtools markdup -r - ${name}_atac_spk_rmdup.bam
	printf "\nnumber of spike-in MAPQ-filtered and deduplicated reads \n"
	samtools view -c ${name}_atac_spk_rmdup.bam
	rm ${name}_atac_spk.sam

	# convert the bam file into a fastq file for further processing
	samtools sort -n ${name}_atac_spk_rmdup.bam -o ${name}_atac_spk_sort.bam
	bedtools bamtofastq -i ${name}_atac_spk_sort.bam -fq ${name}_atac_spk_r1.fastq -fq2 ${name}_atac_spk_r2.fastq
	rm ${name}_atac_spk_rmdup.bam ${name}_atac_spk_sort.bam
	
	# align processed read data to 'sample' species
	# here we are aligning de-duplicated and 'uniquely' mapped reads from alignment to the spike-in species
	# the reads aligned here are the 'intersect' reads  
	# these reads will be removed from the spike-in and re-aligned with no mismatches to seed
	fastq_pair ${name}_atac_spk_r1.fastq ${name}_atac_spk_r2.fastq
	printf "\nalign spike-aligned and MAPQ-filtered data to 'sample' species \n"
	bowtie2 --maxins 500 -x /Volumes/GUERTIN_2/adipogenesis/atac/mm10 -1 ${name}_atac_spk_r1.fastq.paired.fq -2 ${name}_atac_spk_r2.fastq.paired.fq -S ${name}_atac_spk.sam
	
	samtools view -b -q 10 ${name}_atac_spk.sam | samtools sort -n - | samtools fixmate -m - - | samtools sort - | samtools markdup -r - ${name}_atac_spk_rmdup.bam
	printf "\nnumber of spike MAPQ-filtered and deduplicated reads \n"
	samtools view -c ${name}_atac_spk_rmdup.bam
	samtools sort -n ${name}_atac_spk_rmdup.bam -o ${name}_atac_spk_sort.bam
	bedtools bamtofastq -i ${name}_atac_spk_sort.bam -fq ${name}_atac_spkrm_r1.fastq -fq2 ${name}_atac_spkrm_r2.fastq
	rm ${name}_atac_spk.sam ${name}_atac_spk_rmdup.bam ${name}_atac_spk_sort.bam *single*

	# remove doubly unique (intersect) reads from spike-in
	# input 1 to python script - reads will be removed from this file
	# input 2 to python script - these reads will be removed
	# input 3 to python script - output file
	fastq_pair ${name}_atac_spkrm_r1.fastq ${name}_atac_spkrm_r2.fastq

	printf "\nremove intersect from sample \n"
	python removeOverlap2.py ${name}_atac_spk_r1.fastq.paired.fq ${name}_atac_spkrm_r1.fastq.paired.fq ${name}_atac_spk_r1_smprem.fastq
	python removeOverlap2.py ${name}_atac_spk_r2.fastq.paired.fq ${name}_atac_spkrm_r2.fastq.paired.fq ${name}_atac_spk_r2_smprem.fastq
	printf "\nline count after removing overlap reads \n"
	python pyShell.py ${name}_atac_spk_r1_smprem.fastq
	python pyShell.py ${name}_atac_spk_r2_smprem.fastq
	rm ${name}_atac_spkrm_r1.fastq ${name}_atac_spkrm_r2.fastq
	rm ${name}_atac_spk_r1.fastq ${name}_atac_spk_r2.fastq
	
	# remove reads from the no-mm spike-in that align without mismatches to the sample
	# input 1 to python script - reads will be removed from this file
	# input 2 to python script - these reads will be removed 
	# input 3 to python script - output file
	printf "\nremove sample in no-mm from spike no-mm reads to get reads that should go back to the spike-in \n"
	python backToSample2.py ${name}_atac_smp_spkrem_spk_r1.fastq.paired.fq ${name}_atac_smp_spkrem_smp_r1.fastq.paired.fq backtospike_r1.fastq
	python backToSample2.py ${name}_atac_smp_spkrem_spk_r2.fastq.paired.fq ${name}_atac_smp_spkrem_smp_r2.fastq.paired.fq backtospike_r2.fastq
	cat ${name}_atac_spk_r1_smprem.fastq backtospike_r1.fastq > finalspike_r1.fastq
	cat ${name}_atac_spk_r2_smprem.fastq backtospike_r2.fastq > finalspike_r2.fastq
	rm ${name}_atac_smp_spkrem_spk_r1.fastq ${name}_atac_smp_spkrem_smp_r1.fastq
	rm ${name}_atac_smp_spkrem_spk_r2.fastq ${name}_atac_smp_spkrem_smp_r2.fastq
	
	# re-align final spike-in reads and generate the output bam file
	fastq_pair finalspike_r1.fastq finalspike_r2.fastq
	printf "\nfinal alignment for spike-in data \n"
	bowtie2 --maxins 500 -x /Volumes/GUERTIN_2/adipogenesis/atac/dm6 -1 finalspike_r1.fastq.paired.fq -2 finalspike_r2.fastq.paired.fq -S ${name}_atac_spike1.sam
	samtools view -b -q 10 ${name}_atac_spike1.sam | samtools sort -n - | samtools fixmate -m - - | samtools sort - | samtools markdup -r - ${name}_atac_spike_atac.bam
	printf "\nfinal number of spike reads \n"
	samtools view -c ${name}_atac_spike_atac.bam
	rm *_r1.fastq *_r2.fastq *_smprem.fastq ${name}_atac_spike1.sam nomm_intersects.fastq
	cd ..
done
