library(DESeq2)
library(DEGreport)
library(tibble)
library(lattice)
load('pro_cluster_rlog_pval_0.00001.Rdata')
     
clusters.all.test.0.00001 <- degPatterns(cluster_rlog, metadata = meta, minc = 100, time = "sample.conditions", col=NULL, eachStep = TRUE)

save(clusters.all.test.0.00001, file = '191130_clusters.all.pro.minc100.pval0.00001.Rdata')
