library(bigWig)
library(DESeq2)
library(lattice)
library(DEGreport)
source('https://raw.githubusercontent.com/mjg54/znf143_pro_seq_analysis/master/docs/ZNF143_functions.R')

directory = '~/Adipogenesis/DNase/'
setwd(directory)
dnase.file = read.table("TERT4_DNase_peaks.bed", header = F, sep = "\t")

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

df.dnase = get.raw.counts.interval(dnase.file, directory, file.prefix = '')
save(df.dnase, file= 'df.dnase.Rdata')

sample.conditions = factor(sapply(strsplit(as.character(colnames(df.dnase)), '_'), '[', 2))

sample.conditions = factor(sample.conditions, levels=c("0h","4h"))

deseq.counts.table = DESeqDataSetFromMatrix(df.dnase, as.data.frame(sample.conditions), ~ sample.conditions)

dds = DESeq(deseq.counts.table)

                                        #counts table
normalized.counts.dnase = counts(dds, normalized=TRUE)
save(normalized.counts.dnase,file='normalized.counts.dnase.Rdata')

                                        #PCA
rld = rlog(dds, blind=TRUE)

x = plotPCA(rld, intgroup="sample.conditions", returnData=TRUE)
plotPCAlattice(x, file = 'PCA_dnase.pdf')

                                        #DESeq
unt = 2
trt = 2

hg.dnase = DESeqDataSetFromMatrix(df.dnase, DataFrame(sample.conditions), ~ sample.conditions)

sizeFactors(hg.dnase) = estimateSizeFactorsForMatrix(df.dnase)
hg.dnase = estimateDispersions(hg.dnase)
hg.dnase = nbinomWaldTest(hg.dnase)
res.hg.dnase = as.data.frame(results(hg.dnase))

lattice = categorize.deseq.df(res.hg.dnase, fdr = 0.01, log2fold = 0.0, treat = '')

# very few differentially accessible peaks
