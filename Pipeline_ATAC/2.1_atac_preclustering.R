library(bigWig)
library(DESeq2)
library(lattice)
library(DEGreport)
source('https://raw.githubusercontent.com/mjg54/znf143_pro_seq_analysis/master/docs/ZNF143_functions.R')

directory = '/scratch/bhn9by/ATAC'
setwd(directory)
preadipo.file = read.table("/scratch/bhn9by/ATAC/old_peak_calling_macs/old_peak_calling_summit_window.bed", header = F, sep = "\t")

get.raw.counts.interval <- function(df, path.to.bigWig, file.prefix = 'H') {
    df = df[,1:5]
    vec.names = c()
    inten.df=data.frame(matrix(ncol = 0, nrow = nrow(df)))
    for (mod.bigWig in Sys.glob(file.path(path.to.bigWig,
                                          paste(file.prefix, "*.bigWig", sep ='')))) {
        factor.name = strsplit(strsplit(mod.bigWig, "/")[[1]][length(strsplit(mod.bigWig, "/")[[1]])], '\\.')[[1]][1]
        print(factor.name)
        vec.names = c(vec.names, factor.name)
        loaded.bw = load.bigWig(mod.bigWig)
        mod.inten = bed.region.bpQuery.bigWig(loaded.bw, df[,1:3])
        inten.df = cbind(inten.df, mod.inten)
    }
    colnames(inten.df) = vec.names
    r.names = paste(df[,1], ':', df[,2], '-', df[,3], sep='')
    row.names(inten.df) = r.names
    return(inten.df)
}

df.preadipo = get.raw.counts.interval(preadipo.file, directory, file.prefix = '3')
save(df.preadipo, file= 'df.preadipo.Rdata')

sample.conditions = factor(sapply(strsplit(as.character(colnames(df.preadipo)), '_'), '[', 2))

sample.conditions = factor(sample.conditions, levels=c("t0","20min","40min","60min","2hr","3hr","4hr"))

deseq.counts.table = DESeqDataSetFromMatrix(df.preadipo, as.data.frame(sample.conditions), ~ sample.conditions)

dds = DESeq(deseq.counts.table)

#counts table
normalized.counts.atac = counts(dds, normalized=TRUE)
save(normalized.counts.atac,file='normalized.counts.atac.Rdata')

#PCA
rld = rlog(dds, blind=TRUE)

x = plotPCA(rld, intgroup="sample.conditions", returnData=TRUE)
plotPCAlattice(x, file = 'PCA_atac.pdf')

#clustering
dds.lrt = DESeq(dds, test="LRT", reduced = ~ 1)

res.lrt = results(dds.lrt)
save(res.lrt, file = 'res.lrt.Rdata')

padj.cutoff = 0.00000001 #1e-8

siglrt.re = res.lrt[res.lrt$padj < padj.cutoff & !is.na(res.lrt$padj),]

rld_mat <- assay(rld)
cluster_rlog = rld_mat[rownames(siglrt.re),]
meta = as.data.frame(sample.conditions)
rownames(meta) = colnames(cluster_rlog)
save(cluster_rlog, meta, sample.conditions, file = 'cluster_rlog_pval_1e8.Rdata')

#generate nondynamic peaks set for FIMO
not.different = rownames(res.lrt[res.lrt$padj > 0.5 & !is.na(res.lrt$padj) & !is.na(res.lrt$log2FoldChange) & abs(res.lrt$log2FoldChange) < 0.25,])
#not.different = rownames(res.lrt[res.lrt$padj > 0.1 & !is.na(res.lrt$padj) & !is.na(res.lrt$baseMean) & res.lrt$baseMean > 10,])
chr = sapply(strsplit(not.different, ':'), '[', 1)
x = sapply(strsplit(not.different, ':'), '[', 2)
start = sapply(strsplit(x, '-'), '[', 1)
end = sapply(strsplit(x, '-'), '[', 2)

curated.not.different = data.frame(chr,start,end)
write.table(curated.not.different,file='nondynamic_peaks.bed',sep='\t',col.names=F,row.names=F,quote=F)

#generate all peaks set
chr = sapply(strsplit(rownames(res.lrt), ':'), '[', 1)
z = sapply(strsplit(rownames(res.lrt), ':'), '[', 2)
start = sapply(strsplit(z, '-'), '[', 1)
end = sapply(strsplit(z, '-'), '[', 2)

bed = data.frame(chr,start,end,res.lrt$baseMean)
write.table(bed[,1:3], file = 'all_peaks.bed',sep='\t',col.names=F,row.names=F,quote=F)

