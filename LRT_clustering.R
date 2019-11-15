library(bigWig)
library(DESeq2)
library(DEGreport)
library(tibble)
library(lattice)
source('https://raw.githubusercontent.com/mjg54/znf143_pro_seq_analysis/master/docs/ZNF143_functions.R')

directory = '/Volumes/GUERTIN_2/adipogenesis/atac/'
setwd(directory)
preadipo.file = read.table("4hour/atac4h_0.05/3T3_atac_summit_200window.bed")

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
                                        #follow this: https://hbctraining.github.io/DGE_workshop/lessons/08_DGE_LRT.html
sample.conditions = factor(sapply(strsplit(as.character(colnames(df.preadipo)), '_'), '[', 2))

sample.conditions = factor(sample.conditions, levels=c("t0","20min","40min","60min","2hr","3hr","4hr"))

deseq.counts.table = DESeqDataSetFromMatrix(df.preadipo, as.data.frame(sample.conditions), ~ sample.conditions)

dds = DESeq(deseq.counts.table)

                                        #exploratory
rld = rlog(dds, blind=TRUE)

                                        #I think two ATAC-seq samples were labeled incorrectly. Contacted GEO, confirmed bwith an old email from Warren.
#
x = plotPCA(rld, intgroup="sample.conditions", returnData=TRUE)
plotPCAlattice(x, file = 'PCA_adipogenesis_lattice_guertin.pdf')

                                        #lrt

dds.lrt = DESeq(dds, test="LRT", reduced = ~ 1)

#use for outside R
normalized.counts = counts(dds, normalized=TRUE)

res.lrt = results(dds_lrt)

padj.cutoff = 0.000001

siglrt.re = res.lrt[res.lrt$padj < padj.cutoff & !is.na(res.lrt$padj),]

dim(siglrt.re)                                    

rld_mat <- assay(rld)
cluster_rlog = rld_mat[rownames(siglrt.re),]
meta = as.data.frame(sample.conditions)
rownames(meta) = colnames(cluster_rlog)

save.image('191115_adipogenesis.Rdata')

#need to play around with this:
clusters <- degPatterns(cluster_rlog[1:1000,], metadata = meta, minc = 200, time = "sample.conditions", col=NULL, eachStep = TRUE)


class(clusters)
head(clusters$df)
cluster_groups <- clusters$df
group1 <- clusters$df %>%
          filter(cluster == 1)








#useful for network inference section 5.6 two wave network
                                        #https://kateto.net/netscix2016.html
#layered graph bipartite
#made a cyclical https://rdrr.io/cran/igraph/man/layout_with_sugiyama.html
#term sugiyama 
#https://stackoverflow.com/questions/51194653/r-is-not-taking-the-parameter-hgap-in-layout-with-sugiyama
