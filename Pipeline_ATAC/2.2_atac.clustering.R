#Rscript should be run on Rivanna due to high memory requirement; see .slurm

library(DESeq2)
library(DEGreport)
library(tibble)
library(lattice)

setwd('/scratch/bhn9by/ATAC')

load('cluster_rlog_pval_1e8.Rdata')
     
clusters.all.test.1e8 <- degPatterns(cluster_rlog, metadata = meta, minc = 100, time = "sample.conditions", col=NULL, eachStep = TRUE)

save(clusters.all.test.1e8, file = 'clusters.all.minc100.1e8.Rdata')
