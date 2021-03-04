#Manually download SraRunTable.txt from https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA550096&o=acc_s%3Aa - 'download Metadata'
#Manually download SRR_Acc_List.txt from same site - 'download Accession List'

#Fix formatting of SraRunTable in R
setwd('C:/School/UVA/Research/Adipogenesis/code_PRO')

df = read.csv('SraRunTable.txt',header=T)

df$names = c(paste0('3T3_t0_rep',1:3,'_pro.fastq.gz'),
             paste0('3T3_20min_rep',1:3,'_pro.fastq.gz'),
             paste0('3T3_2hr_rep',1:3,'_pro.fastq.gz'),
             paste0('3T3_3hr_rep',1:3,'_pro.fastq.gz'),
             paste0('3T3_40min_rep',1:3,'_pro.fastq.gz'),
             paste0('3T3_4hr_rep',1:3,'_pro.fastq.gz'),
             paste0('3T3_60min_rep',1:3,'_pro.fastq.gz'),
             paste0('3T3_6d_rep',1:3,'_pro.fastq.gz')
             )

write.csv(df,'sra.metadata.csv',row.names=F,quote=F)
