library(DESeq2)
library(DEGreport)
library(tibble)
library(lattice)
load('pTA_pro_cluster_rlog_pval_1e5.Rdata')
setwd('~/pro_pTA')     
clusters.all.1e5 <- degPatterns(cluster_rlog.pTA, metadata = meta.pTA, minc = 10, time = "sample.conditions", col=NULL, eachStep = TRUE)

save(clusters.all.1e5, file = '191230_clusters.pTA.pro.minc10.pval1e5.Rdata')
