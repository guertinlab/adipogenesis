#!/bin/bash
#this is complicated in order to deal with Drosophila reads.
#However, since they Drosophila reads that erroneously align to human are so minimal,
#this complicated alignment scheme is not necessary.
#Note that Drosophila reads were spiked in to experiment,
#but not used in the analysis due to low abundance.
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
#align to mouse genome
	bowtie2 --maxins 500 -x mm10 -1 ${name}_atac_PE1.fastq.paired.fq \
		-2 ${name}_atac_PE2.fastq.paired.fq -S ${name}_atac_smp.sam
#quality filter and remove duplicate amplicons
	samtools view -b -q 10 ${name}_atac_smp.sam | samtools sort -n - | \
	    samtools fixmate -m - - | samtools sort - | samtools markdup -r - ${name}_atac_smp_rmdup.bam
	rm ${name}_atac_smp.sam
#align to Drosophila genome (spike in cells) and remove ambiguous dm6 and mm10 alignments
	samtools sort -n ${name}_atac_smp_rmdup.bam -o ${name}_atac_smp_sort.bam
	bedtools bamtofastq -i ${name}_atac_smp_sort.bam -fq ${name}_atac_smp_r1.fastq \
		 -fq2 ${name}_atac_smp_r2.fastq
	rm ${name}_atac_smp_rmdup.bam ${name}_atac_smp_sort.bam
	fastq_pair ${name}_atac_smp_r1.fastq ${name}_atac_smp_r2.fastq
	bowtie2 --maxins 500 -x dm6 -1 ${name}_atac_smp_r1.fastq.paired.fq \
		-2 ${name}_atac_smp_r2.fastq.paired.fq -S ${name}_atac_spk.sam
	samtools view -b -q 10 ${name}_atac_spk.sam | samtools sort -n - | \
	    samtools fixmate -m - - | samtools sort - | samtools markdup -r - ${name}_atac_spk_rmdup.bam
	samtools sort -n ${name}_atac_spk_rmdup.bam -o ${name}_atac_spk_sort.bam
	bedtools bamtofastq -i ${name}_atac_spk_sort.bam \
		 -fq ${name}_atac_spk_r1.fastq -fq2 ${name}_atac_spk_r2.fastq
	rm ${name}_atac_spk.sam ${name}_atac_spk_rmdup.bam ${name}_atac_spk_sort.bam *single*

	fastq_pair ${name}_atac_spk_r1.fastq ${name}_atac_spk_r2.fastq
	python removeOverlap2.py ${name}_atac_smp_r1.fastq.paired.fq \
	       ${name}_atac_spk_r1.fastq.paired.fq ${name}_atac_smp_r1_spkrem.fastq
	python removeOverlap2.py ${name}_atac_smp_r2.fastq.paired.fq \
	       ${name}_atac_spk_r2.fastq.paired.fq ${name}_atac_smp_r2_spkrem.fastq
	python pyShell.py ${name}_atac_smp_r1_spkrem.fastq
	python pyShell.py ${name}_atac_smp_r2_spkrem.fastq
	rm ${name}_atac_smp_r1.fastq ${name}_atac_smp_r2.fastq *single*

#align ambigous reads with no mismatches to assign to the appropriate genome 
	bowtie2 --maxins 500 --no-1mm-upfront -x mm10 -1 ${name}_atac_spk_r1.fastq.paired.fq \
		-2 ${name}_atac_spk_r2.fastq.paired.fq -S ${name}_atac_smp_spkrem_smp.sam

	bowtie2 --maxins 500 --no-1mm-upfront -x dm6 -1 ${name}_atac_spk_r1.fastq.paired.fq \
		-2 ${name}_atac_spk_r2.fastq.paired.fq -S ${name}_atac_smp_spkrem_spk.sam

	samtools view -b -q 10 ${name}_atac_smp_spkrem_smp.sam | samtools sort -n - | \
	    samtools fixmate -m - - | samtools sort - | \
	    samtools markdup -r - ${name}_atac_smp_spkrem_smp_rmdup.bam
	samtools sort -n ${name}_atac_smp_spkrem_smp_rmdup.bam -o ${name}_atac_smp_spkrem_smp_sort.bam
	bedtools bamtofastq -i ${name}_atac_smp_spkrem_smp_sort.bam \
		 -fq ${name}_atac_smp_spkrem_smp_r1.fastq -fq2 ${name}_atac_smp_spkrem_smp_r2.fastq
	rm ${name}_atac_smp_spkrem_smp.sam ${name}_atac_smp_spkrem_smp_rmdup.bam
	rm ${name}_atac_smp_spkrem_smp_sort.bam
	rm ${name}_atac_spk_r1.fastq ${name}_atac_spk_r2.fastq
	
	samtools view -b -q 10 ${name}_atac_smp_spkrem_spk.sam | samtools sort -n - | \
	    samtools fixmate -m - - | samtools sort - | \
	    samtools markdup -r - ${name}_atac_smp_spkrem_spk_rmdup.bam
	samtools sort -n ${name}_atac_smp_spkrem_spk_rmdup.bam -o ${name}_atac_smp_spkrem_spk_sort.bam
	bedtools bamtofastq -i ${name}_atac_smp_spkrem_spk_sort.bam \
		 -fq ${name}_atac_smp_spkrem_spk_r1.fastq -fq2 ${name}_atac_smp_spkrem_spk_r2.fastq
	rm ${name}_atac_smp_spkrem_spk.sam ${name}_atac_smp_spkrem_spk_rmdup.bam
	rm ${name}_atac_smp_spkrem_spk_sort.bam

	fastq_pair ${name}_atac_smp_spkrem_spk_r1.fastq ${name}_atac_smp_spkrem_spk_r2.fastq
	fastq_pair ${name}_atac_smp_spkrem_smp_r1.fastq ${name}_atac_smp_spkrem_smp_r2.fastq
	python backToSample2.py ${name}_atac_smp_spkrem_smp_r1.fastq.paired.fq \
	       ${name}_atac_smp_spkrem_spk_r1.fastq.paired.fq ${name}_atac_smp_spkrem2_r1.fastq
	python backToSample2.py ${name}_atac_smp_spkrem_smp_r2.fastq.paired.fq \
	       ${name}_atac_smp_spkrem_spk_r2.fastq.paired.fq ${name}_atac_smp_spkrem2_r2.fastq
	rm *single*	
	cat ${name}_atac_smp_r1_spkrem.fastq ${name}_atac_smp_spkrem2_r1.fastq > ${name}_atac_smp_r1.fastq
	cat ${name}_atac_smp_r2_spkrem.fastq ${name}_atac_smp_spkrem2_r2.fastq > ${name}_atac_smp_r2.fastq
	rm ${name}_atac_smp_r1_spkrem.fastq ${name}_atac_smp_spkrem2_r1.fastq
	rm ${name}_atac_smp_r2_spkrem.fastq ${name}_atac_smp_spkrem2_r2.fastq
	
#align processed sample reads and convert to bam
	fastq_pair ${name}_atac_smp_r1.fastq ${name}_atac_smp_r2.fastq
	bowtie2 --maxins 500 -x mm10 -1 ${name}_atac_smp_r1.fastq.paired.fq \
		-2 ${name}_atac_smp_r2.fastq.paired.fq -S ${name}_atac_sample1.sam
	samtools view -b -q 10 ${name}_atac_sample1.sam | samtools sort -n - | \
	    samtools fixmate -m - - | samtools sort - | samtools markdup -r - ${name}_atac_sample_atac.bam

	samtools view -c ${name}_atac_sample_atac.bam
	rm ${name}_atac_sample1.sam *single*
#the same process can be applied to spike in reads, but since the proportion
#of spike in reads are so low, we do not prceed with them for this analysis
	cd ..
done
