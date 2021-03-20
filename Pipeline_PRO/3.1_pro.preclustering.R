library(bigWig)
library(DESeq2)
library(DEGreport)
library(tibble)
library(lattice)
library(tidyr)
source('https://raw.githubusercontent.com/mjg54/znf143_pro_seq_analysis/master/docs/ZNF143_functions.R')

dir = '/scratch/bhn9by/PRO/'

setwd(dir)

gene.file = read.table(paste0(dir,'primary_transcript_annotation/primary_transcript_annotation.bed'), header=FALSE)

get.raw.counts.interval.pro <- function(df, path.to.bigWig, file.prefix = 'H', file.suffix = '.bigWig') {
    vec.names = c()
    inten.df=data.frame(matrix(ncol = 0, nrow = nrow(df)))
    
    for (mod.bigWig in Sys.glob(file.path(path.to.bigWig, paste0(file.prefix, "*plus", file.suffix)))) {
        print(mod.bigWig)
        factor.name = strsplit(strsplit(mod.bigWig, "/")[[1]][length(strsplit(mod.bigWig, "/")[[1]])], '_plus')[[1]][1]
        print(factor.name)
        vec.names = c(vec.names, factor.name)
        loaded.bw.plus = load.bigWig(mod.bigWig)
        loaded.bw.minus = load.bigWig(paste0(path.to.bigWig,'/',factor.name,'_minus',file.suffix))
        mod.inten = bed6.region.bpQuery.bigWig(loaded.bw.plus, loaded.bw.minus, df)
        inten.df = cbind(inten.df, mod.inten)
    }
    colnames(inten.df) = vec.names
    r.names = paste(df[,1], ':', df[,2], '-', df[,3],'_', df[,4], sep='')
    row.names(inten.df) = r.names
    return(inten.df)
}

df.3t3 = get.raw.counts.interval.pro(gene.file, dir, file.prefix = '3T3',file.suffix = '_body_0-mer.bigWig')

colnames(df.3t3) = sapply(strsplit(colnames(df.3t3), '3T3_'), '[', 2)

save(df.3t3, file = 'df.3t3.Rdata')

#for normalizing bigWigs for browser
write.table(estimateSizeFactorsForMatrix(df.3t3),file = 'norm.bedGraph.sizeFactor.txt', quote =F, col.names=F)

sample.conditions = factor(sapply(strsplit(as.character(colnames(df.3t3)), '_'), '[', 1))

sample.conditions = factor(sample.conditions, levels=c("t0","20min","40min","60min","2hr","3hr","4hr"))

deseq.counts.table = DESeqDataSetFromMatrix(df.3t3, as.data.frame(sample.conditions), ~ sample.conditions)

dds = DESeq(deseq.counts.table)

#counts table
normalized.counts.pro = counts(dds, normalized=TRUE)
save(normalized.counts.pro,file='normalized.counts.pro.Rdata')

#PCA
rld = rlog(dds, blind=TRUE)

#I think two ATAC-seq samples were labeled incorrectly. Contacted GEO, confirmed with an old email from Warren.
x = plotPCA(rld, intgroup="sample.conditions", returnData=TRUE)
plotPCAlattice(x, file = 'PCA_pro.pdf')

#clustering
dds.lrt = DESeq(dds, test="LRT", reduced = ~ 1)

res.lrt = results(dds.lrt)
save(res.lrt, file = 'res.lrt.Rdata')

padj.cutoff = 1e-40

siglrt.re = res.lrt[res.lrt$padj < padj.cutoff & !is.na(res.lrt$padj),]

rld_mat <- assay(rld)
cluster_rlog = rld_mat[rownames(siglrt.re),]
meta = as.data.frame(sample.conditions)
rownames(meta) = colnames(cluster_rlog)
save(cluster_rlog, meta, sample.conditions, file = 'cluster_rlog_pval_1e40.Rdata')

#define dynamic genes here rather than after clustering b/c ~600 genes not sorted into clusters
a = strsplit(rownames(cluster_rlog),'_')
gene = unlist(a)[2*(1:nrow(cluster_rlog))]
a = unlist(a)[2*(1:nrow(cluster_rlog))-1]
a = strsplit(a,':')
chr = unlist(a)[2*(1:nrow(cluster_rlog))-1]
a = unlist(a)[2*(1:nrow(cluster_rlog))]
a = strsplit(a,'-')
start = unlist(a)[2*(1:nrow(cluster_rlog))-1]
end = unlist(a)[2*(1:nrow(cluster_rlog))]

bed = data.frame(chr=chr,start=start,end=end,gene=gene)

write.table(bed, file = 'dynamic_genes.bed', quote = FALSE, sep = '\t', col.names=FALSE, row.names=FALSE)
