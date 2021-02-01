library(DESeq2)
library(DEGreport)
library(tibble)
library(lattice)
load('cluster_rlog_pval_0.00000001.Rdata')
     
clusters.all.test.0.00000001 <- degPatterns(cluster_rlog, metadata = meta, minc = 100, time = "sample.conditions", col=NULL, eachStep = TRUE)

save(clusters.all.test.0.00000001, file = '191119_clusters.all.minc100.pval0.00000001.Rdata')
