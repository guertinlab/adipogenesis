library(DESeq2)
library(DEGreport)
library(tibble)
library(lattice)

setwd('/scratch/bhn9by/PRO')

load('cluster_rlog_pval_1e40.Rdata')
     
clusters.all.test.1e40 <- degPatterns(cluster_rlog, metadata = meta, minc = 100, time = "sample.conditions", col=NULL, eachStep = TRUE)

save(clusters.all.test.1e40, file = 'clusters.all.minc100.1e40.Rdata')
