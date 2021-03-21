library(DESeq2)

categorize.deseq.df <- function(df, fdr = 0.05, log2fold = 0.0, treat
= 'Auxin') {

    df.effects.lattice = df
    df.effects.lattice$response = 'All Other Genes'
    
    if(nrow(df.effects.lattice[df.effects.lattice$padj < fdr & !is.na(df.effects.lattice$padj) & df.effects.lattice$log2FoldChange > log2fold,]) > 0) {
        df.effects.lattice[df.effects.lattice$padj < fdr & !is.na(df.effects.lattice$padj) & df.effects.lattice$log2FoldChange > log2fold,]$response = 'Activated'
    }
    if(nrow(df.effects.lattice[df.effects.lattice$padj < fdr & !is.na(df.effects.lattice$padj) & df.effects.lattice$log2FoldChange < log2fold,]) > 0) {
        df.effects.lattice[df.effects.lattice$padj < fdr & !is.na(df.effects.lattice$padj) & df.effects.lattice$log2FoldChange < log2fold,]$response = 'Repressed'
    }
    if(nrow(df.effects.lattice[df.effects.lattice$padj > 0.5 & !is.na(df.effects.lattice$padj) & abs(df.effects.lattice$log2FoldChange) < 0.25,]) > 0) {
        df.effects.lattice[df.effects.lattice$padj > 0.5 & !is.na(df.effects.lattice$padj) & abs(df.effects.lattice$log2FoldChange) < 0.25,]$response = 'Unchanged'
    }
    
    return(df.effects.lattice)
}

setwd("/scratch/bhn9by/PRO/")

#dynamic.genes = read.table('dynamic_genes.bed',sep='\t')

pTA = read.table('primary_transcript_annotation/primary_transcript_annotation.bed')

df = data.frame(matrix(nrow=nrow(pTA),ncol=0))

rownames(df) = paste0(pTA[,1],':',pTA[,2],'-',pTA[,3],'_',pTA[,4])

df$chr = pTA[,1]
df$start = pTA[,2]
df$end = pTA[,3]
df$gene = pTA[,4]
df$strand = pTA[,6]

df$location = paste0(df$chr,':',df$start,'-',df$end)

#add TSS
func <- function(gene) {
    if(df[df$gene == gene,]$strand == '+') {
        return(df[df$gene == gene,]$start)
    } else {
        return(df[df$gene == gene,]$end)
    }
}

df$TSS = sapply(df$gene,func)

#add size
df$size = abs(df$start-df$end)

#add counts
load('normalized.counts.pro.Rdata')

zero.min = c()
twenty.min = c()
forty.min = c()
sixty.min = c()
onetwenty.min = c()
oneeighty.min = c()
twoforty.min = c()
genes = c()

print('Getting PRO means')
for (i in 1:nrow(normalized.counts.pro)) {
    zero.min = append(zero.min, mean(normalized.counts.pro[i,19:21]))
    twenty.min = append(twenty.min, mean(normalized.counts.pro[i,1:3]))
    forty.min = append(forty.min, mean(normalized.counts.pro[i,10:12]))
    sixty.min = append(sixty.min, mean(normalized.counts.pro[i,16:18]))
    onetwenty.min = append(onetwenty.min,mean(normalized.counts.pro[i,4:6]))
    oneeighty.min = append(oneeighty.min, mean(normalized.counts.pro[i,7:9]))
    twoforty.min = append(twoforty.min, mean(normalized.counts.pro[i,13:15]))
    genes = append(genes,rownames(normalized.counts.pro)[i])
}

pro.means = data.frame(row.names = genes, min.0 = zero.min, min.20 = twenty.min,
                   min.40 = forty.min, min.60 = sixty.min, min.120 = onetwenty.min,
                   min.180 = oneeighty.min, min.240 = twoforty.min)

save(pro.means,file='pro.means.Rdata')

df = merge(df,pro.means,by='row.names',all=TRUE)
rownames(df) = df[,1]
df = df[,-1]

#add pairwise comparisons
print('Adding Pairwise Comparisons')

load('df.3t3.Rdata')

time.pts.key = data.frame(time.pts=c(rep(20,3),rep(120,3),rep(180,3),rep(40,3),rep(240,3),rep(60,3),rep(0,3)),
                          cols = c(1:21))

time.pts = c(0,20,40,60,120,180,240)

comparisons.df = data.frame(row.names=rownames(df))
for (i in 1:(length(time.pts)-1)) {
    for (j in (i+1):length(time.pts)) {

        print(paste0(time.pts[j],'.v.',time.pts[i]))
        
        a = time.pts.key[time.pts.key$time.pts == time.pts[i],]$cols
        b = time.pts.key[time.pts.key$time.pts == time.pts[j],]$cols
        merged.counts.small = df.3t3[,c(a,b)]

                                        # number of replicates per condition
        unt = 3
        trt = 3

        sample.conditions = factor(c(rep("untreated",unt), rep("treated",trt)), levels=c("untreated","treated"))
        mm.deseq.counts.table = DESeqDataSetFromMatrix(merged.counts.small, DataFrame(sample.conditions), ~ sample.conditions)

        mm.atac = mm.deseq.counts.table
        atac.size.factors = estimateSizeFactorsForMatrix(merged.counts.small)
        
        sizeFactors(mm.atac) = atac.size.factors
        mm.atac = estimateDispersions(mm.atac)
        mm.atac = nbinomWaldTest(mm.atac)
        res.mm.atac = results(mm.atac)
        
        lattice = categorize.deseq.df(res.mm.atac, fdr = 0.001, log2fold = 0.0, treat = '')
        lattice = as.data.frame(lattice[,c(2,6,7)])
        colnames(lattice) = paste0(colnames(lattice),'.',time.pts[j],'.v.',time.pts[i])
        
        comparisons.df = merge(comparisons.df,lattice,by='row.names',all=TRUE)
        rownames(comparisons.df) = comparisons.df$Row.names
        comparisons.df = comparisons.df[,-1]
    }
}

save(comparisons.df,file='comparisons.df.Rdata')

df = merge(df,comparisons.df,by = 'row.names',all = TRUE)
rownames(df) = df$Row.names
df = df[,-1]

#add baseMean and overall time course info
print('Add baseMean and overall time course info')

load('res.lrt.Rdata')
res.lrt = as.data.frame(res.lrt[,c(1,2,6)])
res.lrt$response = 'Nondynamic'
res.lrt[res.lrt$padj < 1e-40 & !is.na(res.lrt$padj),]$response = 'Dynamic'
colnames(res.lrt) = paste0(colnames(res.lrt),'.time.course')

df = merge(df,res.lrt,by = 'row.names',all = TRUE)
rownames(df) = df$Row.names
df = df[,-1]

final.pro.all.df=df
save(final.pro.all.df,file='final.pro.all.df.Rdata')
