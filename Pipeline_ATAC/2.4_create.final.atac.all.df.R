library(DESeq2)

categorize.deseq.df <- function(df, fdr = 0.05, log2fold = 0.0, treat
                                = 'Auxin') {
    
    df.effects.lattice = df
    df.effects.lattice$response = 'All Other Peaks'
    
    if(nrow(df.effects.lattice[df.effects.lattice$padj < fdr & !is.na(df.effects.lattice$padj) & df.effects.lattice$log2FoldChange > log2fold,]) > 0) {
        df.effects.lattice[df.effects.lattice$padj < fdr & !is.na(df.effects.lattice$padj) & df.effects.lattice$log2FoldChange > log2fold,]$response = 'Increased'
    }
    if(nrow(df.effects.lattice[df.effects.lattice$padj < fdr & !is.na(df.effects.lattice$padj) & df.effects.lattice$log2FoldChange < log2fold,]) > 0) {
        df.effects.lattice[df.effects.lattice$padj < fdr & !is.na(df.effects.lattice$padj) & df.effects.lattice$log2FoldChange < log2fold,]$response = 'Decreased'
    }
    if(nrow(df.effects.lattice[df.effects.lattice$padj > 0.5 & !is.na(df.effects.lattice$padj) & abs(df.effects.lattice$log2FoldChange) < 0.25,]) > 0) {
        df.effects.lattice[df.effects.lattice$padj > 0.5 & !is.na(df.effects.lattice$padj) & abs(df.effects.lattice$log2FoldChange) < 0.25,]$response = 'Unchanged'
    }
    
    return(df.effects.lattice)
}

setwd("/scratch/bhn9by/ATAC")

peaks = unique(read.table('all_peaks.bed',sep='\t'))

peaks$location = paste0(peaks$V1,':',peaks$V2,'-',peaks$V3)

df = data.frame(row.names = peaks$location,chr=peaks[,1],start=peaks[,2],end=peaks[,3])

#add counts
load('normalized.counts.atac.Rdata')

zero.min = c()
twenty.min = c()
forty.min = c()
sixty.min = c()
onetwenty.min = c()
oneeighty.min = c()
twoforty.min = c()
genes = c()

print('Getting ATAC means')
for (i in 1:nrow(normalized.counts.atac)) {
    zero.min = append(zero.min, mean(normalized.counts.atac[i,19:21]))
    twenty.min = append(twenty.min, mean(normalized.counts.atac[i,1:3]))
    forty.min = append(forty.min, mean(normalized.counts.atac[i,10:12]))
    sixty.min = append(sixty.min, mean(normalized.counts.atac[i,16:18]))
    onetwenty.min = append(onetwenty.min,mean(normalized.counts.atac[i,4:6]))
    oneeighty.min = append(oneeighty.min, mean(normalized.counts.atac[i,7:9]))
    twoforty.min = append(twoforty.min, mean(normalized.counts.atac[i,13:15]))
    genes = append(genes,rownames(normalized.counts.atac)[i])
}

atac.means = data.frame(row.names = genes, min.0 = zero.min, min.20 = twenty.min,
                   min.40 = forty.min, min.60 = sixty.min, min.120 = onetwenty.min,
                   min.180 = oneeighty.min, min.240 = twoforty.min)

save(atac.means,file='atac.means.Rdata')

df = merge(df,atac.means,by='row.names',all=TRUE)
rownames(df) = df[,1]
df = df[,-1]

#add cluster information
print('Adding Cluster Info')
print(head(df))

#up.flat - 9,23,17,11
#grad.up - 5,8,10,1
#up.down - 6,2,12
#grad.down - 7,3,4,13,14
#down.up - 24

supercluster.df = data.frame(cluster=c(9,23,17,11,5,8,10,1,6,2,12,7,3,4,13,14,24),
                             supercluster=c(rep('up.flat',4),rep('grad.up',4),rep('up.down',3),rep('grad.down',5),rep('down.up',1)))

cluster.df = data.frame()
for (file in Sys.glob(file.path(paste0('cluster_bed_cluster*.bed')))) {
    cluster = as.numeric(strsplit(strsplit(file,'cluster_bed_cluster')[[1]][2],'.bed')[[1]][1])
    print(cluster)
    sc = supercluster.df[supercluster.df$cluster == cluster,]$supercluster
    x = read.table(file)
    x$location = paste0(x$V1,':',x$V2,'-',x$V3)
    cluster.df = rbind(cluster.df,data.frame(peak=x$location,cluster=cluster,supercluster=sc))
}

df = merge(df,cluster.df,by.x='row.names',by.y=1,all.x=TRUE)
rownames(df) = df[,1]
df = df[,-1]

#add tf family scores
print('Adding TF Family Scores')
print(head(df))

scores.df = data.frame()
for (file in Sys.glob('fimo_composites/*_fimo_all.bed')) {
    factor = strsplit(strsplit(file,'/')[[1]][2],'_fimo_all.bed')[[1]][1]
    print(factor)
    x = read.table(file,sep='\t')
    x = x[x$V5 != -1,]
    if(factor %in% c('SP','KLF')) {
        y = aggregate(as.numeric(V7)~V1+V2+V3, data=x, FUN=sum)
    } else {
        y = aggregate(as.numeric(V8)~V1+V2+V3, data=x, FUN=sum)
    }
    colnames(y) = c('chr', 'start', 'end', factor)
    rownames(y) = paste0(y[,1], ':', y[,2], '-', y[,3])
    
    scores.df = rbind(scores.df,data.frame(peak = rownames(y),score=y[,4],factor = factor))
}

fimo.scores.all.atac = data.frame(peak = unique(scores.df$peak))

for(factor in unique(scores.df$factor)) {
    temp = scores.df[scores.df$factor == factor,]
    temp = temp[,c(1,2)]
    fimo.scores.all.atac = merge(fimo.scores.all.atac,temp,by='peak',all.x=TRUE)
    colnames(fimo.scores.all.atac)[ncol(fimo.scores.all.atac)] = factor
}

rownames(fimo.scores.all.atac) = fimo.scores.all.atac$peak
fimo.scores.all.atac = fimo.scores.all.atac[,-1]
save(fimo.scores.all.atac,file='fimo.scores.all.atac.Rdata')

df = merge(df,fimo.scores.all.atac,by = 'row.names',all = TRUE)
rownames(df) = df$Row.names
df = df[,-1]

#add pairwise comparisons
print('Adding Pairwise Comparisons')
print(head(df))

load('df.preadipo.Rdata')

time.pts.key = data.frame(time.pts=c(rep(20,3),rep(120,3),rep(180,3),rep(40,3),rep(240,3),rep(60,3),rep(0,3)),
                          cols = c(1:21))

time.pts = c(0,20,40,60,120,180,240)

comparisons.df = data.frame()
for (i in 1:(length(time.pts)-1)) {
    for (j in (i+1):length(time.pts)) {

        print(paste0(time.pts[j],'.v.',time.pts[i]))
        
        a = time.pts.key[time.pts.key$time.pts == time.pts[i],]$cols
        b = time.pts.key[time.pts.key$time.pts == time.pts[j],]$cols
        merged.counts.small = df.preadipo[,c(a,b)]

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

print('Add baseMean and overall time course info')

load('res.lrt.Rdata')
res.lrt = as.data.frame(res.lrt[,c(1,2,6)])
res.lrt$response = 'Nondynamic'
res.lrt[res.lrt$padj < 0.00000001 & !is.na(res.lrt$padj),]$response = 'Dynamic'
colnames(res.lrt) = paste0(colnames(res.lrt),'.time.course')

df = merge(df,res.lrt,by = 'row.names',all = TRUE)
rownames(df) = paste0(df$chr,':',df$start,'-',df$end)
df = df[,-1]

print('Add distribution')

df$distribution = 'Intergenic'

x = read.table('all_ATAC_peaks_intragenic.bed')
x$location = paste0(x$V1,':',x$V2,'-',x$V3)
df[rownames(df) %in% x$location,]$distribution = 'Intragenic'

x = read.table('all_ATAC_peaks_promoters.bed')
x$location = paste0(x$V1,':',x$V2,'-',x$V3)
df[rownames(df) %in% x$location,]$distribution = 'Promoter'

final.atac.all.df=df
save(final.atac.all.df,file='final.atac.all.df.Rdata')
