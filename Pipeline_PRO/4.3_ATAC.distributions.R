setwd('/scratch/bhn9by/ATAC')

#sort all ATACs peaks into promoters, intragenic, and intergenic

x = read.table('all_ATAC_peaks_promoters_temp.bed')
x$location = paste0(x$V1,':',x$V2,'-',x$V3)

y = read.table('all_ATAC_peaks_intragenic_temp.bed')
y$location = paste0(y$V1,':',y$V2,'-',y$V3)

z = read.table('all_peaks.bed')
z$location = paste0(z$V1,':',z$V2,'-',z$V3)

promoters = x$location
intragenic = setdiff(y$location,x$location)
intergenic = setdiff(z$location,y$location)
intergenic = setdiff(intergenic,x$location)

#write promoter all ATAC peaks to bed file
chr = sapply(strsplit(promoters, ':'), '[', 1)
x = sapply(strsplit(promoters, ':'), '[', 2)
start = sapply(strsplit(x, '-'), '[', 1)
end = sapply(strsplit(x, '-'), '[', 2)
bed = data.frame(chr,start,end)
write.table(bed,file='all_ATAC_peaks_promoters.bed',sep='\t',quote=F,col.names = F,row.names = F)

#write intragenic all ATAC peaks to bed file
chr = sapply(strsplit(intragenic, ':'), '[', 1)
x = sapply(strsplit(intragenic, ':'), '[', 2)
start = sapply(strsplit(x, '-'), '[', 1)
end = sapply(strsplit(x, '-'), '[', 2)
bed = data.frame(chr,start,end)
write.table(bed,file='all_ATAC_peaks_intragenic.bed',sep='\t',quote=F,col.names = F,row.names = F)

#write intergenic all ATAC peaks to bed file
chr = sapply(strsplit(intergenic, ':'), '[', 1)
x = sapply(strsplit(intergenic, ':'), '[', 2)
start = sapply(strsplit(x, '-'), '[', 1)
end = sapply(strsplit(x, '-'), '[', 2)
bed = data.frame(chr,start,end)
write.table(bed,file='all_ATAC_peaks_intergenic.bed',sep='\t',quote=F,col.names = F,row.names = F)
