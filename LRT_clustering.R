library(bigWig)
library(DESeq2)
library(DEGreport)
library(tibble)
library(lattice)
library(tidyr)
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
save(df.preadipo, file= 'df.preadipo.Rdata')

load('df.preadipo.Rdata')
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

res.lrt = results(dds.lrt)

padj.cutoff = 0.00000001 #1e-8

siglrt.re = res.lrt[res.lrt$padj < padj.cutoff & !is.na(res.lrt$padj),]

dim(siglrt.re)                                    

rld_mat <- assay(rld)
cluster_rlog = rld_mat[rownames(siglrt.re),]
meta = as.data.frame(sample.conditions)
rownames(meta) = colnames(cluster_rlog)

save.image('191119_adipogenesis.Rdata')
save(cluster_rlog, meta, sample.conditions, file = 'cluster_rlog_pval_0.00000001.Rdata')


####run on Rivanna: LRT_rivanna_slurm.slurm

#!/bin/bash
#SBATCH -n 1
#SBATCH -t 96:00:00
#SBATCH -o minC100_191119_1.out
#SBATCH -p largemem
#SBATCH -A guertinlab

cd ~
module load gcc R
Rscript LRT_largemem_clustering.R


###### LRT_largemem_clustering.R
library(DESeq2)
library(DEGreport)
library(tibble)
library(lattice)
load('cluster_rlog_pval_0.00000001.Rdata')
     
clusters.all.test.0.00000001 <- degPatterns(cluster_rlog, metadata = meta, minc = 100, time = "sample.conditions", col=NULL, eachStep = TRUE)

save(clusters.all.test.0.00000001, file = '191119_clusters.all.minc100.pval0.00000001.Rdata')

######
#https://stackoverflow.com/questions/51295402/r-on-macos-error-vector-memory-exhausted-limit-reached
#https://stackoverflow.com/questions/51248293/error-vector-memory-exhausted-limit-reached-r-3-5-0-macos
                                        #this output can be used to make custom graphics:
                                        #see https://rdrr.io/bioc/DEGreport/man/degPatterns.html for the list of outputs:



#load data from Rivanna
load('191119_clusters.all.minc100.pval0.00000001.Rdata')


#I am using plotiing_clusters_LRT.R to optimize plotting functions with a small data frame

plot.df = clusters.all.test.0.00000001$normalized

plot.df$sample.conditions = as.character(plot.df$sample.conditions)
plot.df$sample.conditions[plot.df$sample.conditions == 't0'] = 0
plot.df$sample.conditions[plot.df$sample.conditions == '20min'] = 20
plot.df$sample.conditions[plot.df$sample.conditions == '40min'] = 40
plot.df$sample.conditions[plot.df$sample.conditions == '60min'] = 60
plot.df$sample.conditions[plot.df$sample.conditions == '2hr'] = 120
plot.df$sample.conditions[plot.df$sample.conditions == '3hr'] = 180
plot.df$sample.conditions[plot.df$sample.conditions == '4hr'] = 240
plot.df$sample.conditions = as.numeric(plot.df$sample.conditions)
plot.df = plot.df[order(plot.df$genes),]
plot.df = plot.df[order(plot.df$sample.conditions),]

plot.df$cluster = paste('cluster', as.character(plot.df$cluster), sep = '')



#this has a bwplot, loess curve, points, and lines, which can be removed from the panel function to customize
pdf(paste('Clusters_minc100_stringent_pvalue','.pdf', sep=''), width=11, height=15)

trellis.par.set(box.umbrella = list(lty = 1, col="black", lwd=1),
                box.rectangle = list( lwd=1.0, col="black", alpha = 1.0),
                plot.symbol = list(col="black", lwd=1.0, pch ='.'))
print(
xyplot(value ~  sample.conditions | cluster, group = genes, data = plot.df, type = c('l'),#type = c('l','p'),
       scales=list(x=list(cex=1.0,relation = "free", rot = 45), y =list(cex=1.0, relation="free")),
       aspect=1.0,
       between=list(y=0.5, x=0.5),
       ylab = list(label = 'Normalized ATAC signal', cex =1.0),
       xlab = list(label = 'Time (minutes)', cex =1.0),
       par.settings = list(superpose.symbol = list(pch = c(16),
                                                   col=c('grey20'), cex =0.5),
                           strip.background=list(col="grey80"),
                           superpose.line = list(col = c('#99999980'), lwd=c(1),
                                                 lty = c(1))),
       panel = function(x, y, ...) {
           panel.xyplot(x, y, ...)
           panel.bwplot(x, y, pch = '|', horizontal = FALSE, box.width = 15, do.out = FALSE)
           panel.loess(x, y, ..., col = "blue", lwd =2.0,  span = 1/2, degree = 1, family = c("gaussian"))
           
})

      )
dev.off()

                                        #just plot loess



                                        #Propose using these data and plotting on a single plot for clustering
loess.df = data.frame(matrix(ncol = 3, nrow = 0),  stringsAsFactors = FALSE)
colnames(loess.df) = c('ATAC','time','cluster')
for (i in unique(plot.df$cluster)){
    ATAC = loess.smooth(plot.df[plot.df$cluster == i,]$sample.conditions, plot.df[plot.df$cluster == i,]$value,  span = 1/2, degree = 1, family = c("gaussian"))$y
    time = loess.smooth(plot.df[plot.df$cluster == i,]$sample.conditions, plot.df[plot.df$cluster == i,]$value, span = 1/2, degree = 1, family = c("gaussian"))$x
    new = as.data.frame(cbind(ATAC, time, i),  stringsAsFactors = FALSE)
    colnames(new) = c('ATAC','time','cluster')
    new$ATAC = as.numeric(new$ATAC)
    new$time = as.numeric(new$time)
    loess.df = rbind(loess.df, new)
}

colors.clusters = rep('white', 21)
colors.clusters[c(6,13,14)] = 'red'
colors.clusters[c(19, 15,16,20,3)] = 'blue'
colors.clusters[c(10)] = 'green'
colors.clusters[c(21,18,1,4)] = 'black'
colors.clusters[c(12,8,9,5,7)] = 'purple'
colors.clusters[c(17, 2, 11)] = 'orange'



pdf(paste('Clusters_individual_loess','.pdf', sep=''), width=4, height=4)

trellis.par.set(box.umbrella = list(lty = 1, col="black", lwd=1),
                box.rectangle = list( lwd=1.0, col="black", alpha = 1.0),
                plot.symbol = list(col="black", lwd=1.0, pch ='.'))
print(
xyplot(ATAC ~  time, groups = cluster, data = loess.df, type = c('l'),#type = c('l','p'),
       scales=list(x=list(cex=1.0,relation = "free", rot = 45), y =list(cex=1.0, relation="free")),
       aspect=1.0,
       between=list(y=0.5, x=0.5),
       ylab = list(label = 'cluster loess of ATAC signal', cex =1.0),
       xlab = list(label = 'Time (minutes)', cex =1.0),
       par.settings = list(superpose.symbol = list(pch = c(16),
                                                   col=colors.clusters, cex =0.5),
                           strip.background=list(col="grey80"),
                           superpose.line = list(col = colors.clusters, lwd=c(3),
                                                 lty = c(1))))

      )
dev.off()
    




plot.df$chr = sapply(strsplit(plot.df$genes, ':'), '[', 1)
range.chr = sapply(strsplit(plot.df$genes, ':'), '[', 2)
plot.df$start = sapply(strsplit(range.chr, '-'), '[', 1)
plot.df$end = sapply(strsplit(range.chr, '-'), '[', 2)


#I just noticed that I screwed this up, but no big deal because MEME skips the duplicate sequence names
#I corrected the code here, but this is not what is running in MEME, although the outcome should be identical
for (i in unique(plot.df$cluster)) {
    print(i)
    write.table(plot.df[plot.df$cluster == i,
                        c('chr','start','end', 'value', 'cluster')][!duplicated(plot.df[plot.df$cluster == i,]$genes),],
                file = paste0('cluster_bed_',
                              gsub(" ", "", i, fixed = TRUE),'.bed'),
                quote = FALSE, row.names = FALSE, col.names = FALSE, sep = '\t')
}


#use these for MEME in Rivanna
#191122_meme_all_cluster_1.slurm
#191122_meme_all_cluster_not1.slurm
#191122_meme_cluster_1wc.sh
#191122_meme_all_cluster_1wc.slurm
#191122_meme_cluster_1.sh
#191122_meme_cluster_not1.sh



#long vs wide conversion http://www.cookbook-r.com/Manipulating_data/Converting_data_between_wide_and_long_format/
library(data.table)

plot.df.cluster = dcast(plot.df, genes + cluster ~ sample.conditions, value.var="value")

avg.clusters = as.data.frame(matrix(nrow = 0, ncol = 7))
colnames(y) = colnames(plot.df.cluster[3:9])
for (i in unique(plot.df.cluster$cluster)) {
    z = data.frame(matrix(colMeans(plot.df.cluster[plot.df.cluster$cluster == i,3:9]), ncol = 7, nrow = 1))
    rownames(z) = c(i)
    colnames(z) = as.character(colnames(plot.df.cluster)[3:9])
    avg.clusters = rbind(avg.clusters, z)
}


dd = dist(avg.clusters)
hc = hclust(dd, method = "complete")

                                        #8

pdf(paste('dendrogram_HC_eight','.pdf', sep=''), width=8, height=5)
plot(hc, xlab = "Clusters", main = ' ', hang = -1)
abline(h = 1.38, lty = 2)
dev.off()

pdf(paste('Clusters_minc100_stringent_pvalue_post_HC_eight','.pdf', sep=''), width=11, height=20)

trellis.par.set(box.umbrella = list(lty = 1, col="black", lwd=1),
                box.rectangle = list( lwd=1.0, col="black", alpha = 1.0),
                plot.symbol = list(col="black", lwd=1.0, pch ='.'))
print(
    xyplot(value ~  sample.conditions | cluster, group = genes, data = plot.df, type = c('l'),#type = c('l','p'),
       scales=list(x=list(cex=1.0,relation = "free", rot = 45), y =list(cex=1.0, relation="free")),
       aspect=1.0,
       layout = c(5,8),
       between=list(y=0.5, x=0.5),
       index.cond=list(
           c(6,13,14,
             19, 15,16,20,3,
             10,
             21,18,
             1,4,
             12,8,9,5,7,
             17,
             2, 11)),
       skip = c(F, F, F, T, T,
                F, F, F, F, F,
                F, T, T, T, T,
                F, F, T, T, T,
                F, F, T, T, T,
                F, F, F, F, F,
                F, T, T, T, T,
                F, F, T, T, T) ,
       ylab = list(label = 'Normalized ATAC signal', cex =1.0),
       xlab = list(label = 'Time (minutes)', cex =1.0),
       par.settings = list(superpose.symbol = list(pch = c(16),
                                                   col=c('grey20'), cex =0.5),
                           strip.background=list(col="grey80"),
                           superpose.line = list(col = c('#99999980'), lwd=c(1),
                                                 lty = c(1))),
       panel = function(x, y, ...) {
           panel.xyplot(x, y, ...)
           panel.bwplot(x, y, pch = '|', horizontal = FALSE, box.width = 15, do.out = FALSE)
           panel.loess(x, y, ..., col = "blue", lwd =2.0,  span = 1/2, degree = 1, family = c("gaussian"))
           
})

      )
dev.off()

                                        #6

pdf(paste('dendrogram_HC_six','.pdf', sep=''), width=8, height=5)
plot(hc, xlab = "Clusters", main = ' ', hang = -1)
abline(h = 2, lty = 2)
dev.off()

pdf(paste('Clusters_minc100_stringent_pvalue_post_HC_six','.pdf', sep=''), width=11, height=16)

trellis.par.set(box.umbrella = list(lty = 1, col="black", lwd=1),
                box.rectangle = list( lwd=1.0, col="black", alpha = 1.0),
                plot.symbol = list(col="black", lwd=1.0, pch ='.'))
print(
    xyplot(value ~  sample.conditions | cluster, group = genes, data = plot.df, type = c('l'),#type = c('l','p'),
       scales=list(x=list(cex=1.0,relation = "free", rot = 45), y =list(cex=1.0, relation="free")),
       aspect=1.0,
       layout = c(5,6),
       between=list(y=0.5, x=0.5),
       index.cond=list(
           c(6,13,14,
             19, 15,16,20,3,
             10,
             21,18,1,4,
             12,8,9,5,7,
             17, 2, 11)),
       skip = c(F, F, F, T, T,
                F, F, F, F, F,
                F, T, T, T, T,
                F, F, F, F, T,
                F, F, F, F, F,
                F, F, F, T, T) ,
       ylab = list(label = 'Normalized ATAC signal', cex =1.0),
       xlab = list(label = 'Time (minutes)', cex =1.0),
       par.settings = list(superpose.symbol = list(pch = c(16),
                                                   col=c('grey20'), cex =0.5),
                           strip.background=list(col="grey80"),
                           superpose.line = list(col = c('#99999980'), lwd=c(1),
                                                 lty = c(1))),
       panel = function(x, y, ...) {
           panel.xyplot(x, y, ...)
           panel.bwplot(x, y, pch = '|', horizontal = FALSE, box.width = 15, do.out = FALSE)
           panel.loess(x, y, ..., col = "blue", lwd =2.0,  span = 1/2, degree = 1, family = c("gaussian"))
           
})

      )
dev.off()


                                        #choose 6 or 8
down.flat = plot.df[plot.df$cluster == 'cluster14' |
                    plot.df$cluster == 'cluster22' |
                    plot.df$cluster == 'cluster23',]

gradual.down = plot.df[plot.df$cluster == 'cluster7' |
                    plot.df$cluster == 'cluster3' |
                    plot.df$cluster == 'cluster4' |
                    plot.df$cluster == 'cluster8' |
                    plot.df$cluster == 'cluster11',]

down.up = plot.df[plot.df$cluster == 'cluster18',]

gradual.up = plot.df[plot.df$cluster == 'cluster9' |
                    plot.df$cluster == 'cluster6' |
                    plot.df$cluster == 'cluster1' |
                    plot.df$cluster == 'cluster12',]

up.flat = plot.df[plot.df$cluster == 'cluster21' |
                    plot.df$cluster == 'cluster16' |
                    plot.df$cluster == 'cluster17' |
                    plot.df$cluster == 'cluster13' |
                    plot.df$cluster == 'cluster15',]

up.down = plot.df[plot.df$cluster == 'cluster5' |
                    plot.df$cluster == 'cluster10' |
                    plot.df$cluster == 'cluster2',]

down.flat$cluster.final = 'down.flat'
gradual.down$cluster.final = 'gradual.down'
down.up$cluster.final = 'down.up'
gradual.up$cluster.final = 'gradual.up'
up.flat$cluster.final = 'up.flat'
up.down$cluster.final = 'up.down'

nrow(down.flat)/7
nrow(gradual.down)/7
nrow(down.up)/7
nrow(gradual.up)/7                       #
nrow(up.flat)/7
nrow(up.down)/7

final.plot.df = rbind(down.flat,
      gradual.down,
      down.up,
      gradual.up,
      up.flat,
      up.down)

pdf(paste('Clusters_SIX_stringent_pvalue','.pdf', sep=''), width=6.83, height=5)

trellis.par.set(box.umbrella = list(lty = 1, col="black", lwd=1),
                box.rectangle = list( lwd=1.0, col="black", alpha = 1.0),
                plot.symbol = list(col="black", lwd=1.0, pch ='.'))
print(
xyplot(value ~  sample.conditions | cluster.final, group = genes, data = final.plot.df, type = c('l'),#type = c('l','p'),
       scales=list(x=list(cex=1.0,relation = "free", rot = 45), y =list(cex=1.0, relation="free")),
       aspect=1.0,
       between=list(y=0.5, x=0.5),
       layout = c(3,2),
       index.cond=list(
           c(1,2,3,6,5,4)),
       ylab = list(label = 'Normalized ATAC signal', cex =1.0),
       xlab = list(label = 'Time (minutes)', cex =1.0),
       par.settings = list(superpose.symbol = list(pch = c(16),
                                                   col=c('grey20'), cex =0.5),
                           strip.background=list(col="grey80"),
                           superpose.line = list(col = c('#99999980'), lwd=c(1),
                                                 lty = c(1))),
       panel = function(x, y, ...) {
           panel.xyplot(x, y, ...)
           panel.bwplot(x, y, pch = '|', horizontal = FALSE, box.width = 15, do.out = FALSE)
           panel.loess(x, y, ..., col = "blue", lwd =2.0,  span = 1/2, degree = 1, family = c("gaussian"))
           
})

      )
dev.off()



for (i in unique(final.plot.df$cluster.final)) {
    print(i)
    write.table(final.plot.df[final.plot.df$cluster.final == i,
                              c('chr','start','end', 'value', 'cluster.final')][!duplicated(final.plot.df[final.plot.df$cluster.final == i,]$genes),],
                file = paste0('final_cluster_bed_',
                              gsub(" ", "", i, fixed = TRUE),'.bed'),
                quote = FALSE, row.names = FALSE, col.names = FALSE, sep = '\t')
}



#next thing I want to do is cluster motifs identified de novo using the method from our 2017 PLOS Gen paper
python /Users/guertinlab/pyscripts/MEME_individual_from_db.py -i /Volumes/GUERTIN_2/adipogenesis/atac/twenty_one_clusters/cluster1.meme_output/meme.txt

setwd('/Volumes/GUERTIN_2/adipogenesis/atac/meme_minimal_atac/mapping_list')



                                        #DEseq2 pTA output

library(bigWig)
library(DESeq2)
library(DEGreport)
library(tibble)
library(lattice)
library(tidyr)
source('https://raw.githubusercontent.com/mjg54/znf143_pro_seq_analysis/master/docs/ZNF143_functions.R')



setwd('/Volumes/GUERTIN_2/adipogenesis/pro_dREG')

direc = '/Volumes/GUERTIN_2/adipogenesis/pro_dREG'

gene.file = read.table('/Volumes/GUERTIN_2/adipogenesis/primary_transcript_annotation.bed', header=FALSE)

get.raw.counts.interval.pro <- function(df, path.to.bigWig, file.prefix = 'H') {
    vec.names = c()
    inten.df=data.frame(matrix(ncol = 0, nrow = nrow(df)))
    
    for (mod.bigWig in Sys.glob(file.path(path.to.bigWig, paste(file.prefix, "*plus.bigWig", sep ='')))) {
        factor.name = strsplit(strsplit(mod.bigWig, "/")[[1]][length(strsplit(mod.bigWig, "/")[[1]])], '_plus')[[1]][1]
        print(factor.name)
        vec.names = c(vec.names, factor.name)
        loaded.bw.plus = load.bigWig(mod.bigWig)
        print(mod.bigWig)
        print(paste(path.to.bigWig,'/',factor.name, '_minus.bigWig', sep=''))
        loaded.bw.minus = load.bigWig(paste(path.to.bigWig,'/',factor.name, '_minus.bigWig', sep=''))
        mod.inten = bed6.region.bpQuery.bigWig(loaded.bw.plus, loaded.bw.minus, df)
        inten.df = cbind(inten.df, mod.inten)
    }
    colnames(inten.df) = vec.names
    r.names = paste(df[,1], ':', df[,2], '-', df[,3],'_', df[,4], sep='')
    row.names(inten.df) = r.names
    return(inten.df)
}


df.3t3 = get.raw.counts.interval.pro(gene.file, direc, file.prefix = 'G')


colnames(df.3t3) = sapply(strsplit(colnames(df.3t3), '3T3_t'), '[', 2)
colnames(df.3t3) = sapply(strsplit(colnames(df.3t3), '_pro'), '[', 1)

save(df.3t3, file = 'adipogenesis.df.3t3.Rdata')

counts.pta = df.3t3[,1:21]
                                      #follow this: https://hbctraining.github.io/DGE_workshop/lessons/08_DGE_LRT.html
sample.conditions = factor(sapply(strsplit(as.character(colnames(counts.pta)), '_'), '[', 1))

sample.conditions = factor(sample.conditions, levels=c("0","20min","40min","60min","2h","3h","4h"))

deseq.counts.table = DESeqDataSetFromMatrix(counts.pta, as.data.frame(sample.conditions), ~ sample.conditions)

dds = DESeq(deseq.counts.table)

                                        #exploratory
rld = rlog(dds, blind=TRUE)

#
x = plotPCA(rld, intgroup="sample.conditions", returnData=TRUE)
plotPCAlattice(x, file = 'PCA_adipogenesis_pro_pTA_lattice_guertin.pdf')

                                        #lrt

dds.lrt = DESeq(dds, test="LRT", reduced = ~ 1)

#use for outside R
normalized.counts = counts(dds, normalized=TRUE)

res.lrt = results(dds.lrt)

padj.cutoff = 0.00001 #1e-25

siglrt.re = res.lrt[res.lrt$padj < padj.cutoff & !is.na(res.lrt$padj),]

dim(siglrt.re)                                    



rld_mat <- assay(rld)
cluster_rlog.pTA = rld_mat[rownames(siglrt.re),]
meta.pTA = as.data.frame(sample.conditions)
rownames(meta.pTA) = colnames(cluster_rlog.pTA)

save(cluster_rlog.pTA, meta.pTA, sample.conditions, file = 'pTA_pro_cluster_rlog_pval_1e5.Rdata')

library(DESeq2)
library(DEGreport)
library(tibble)
library(lattice)

#clusters.all.1e25 <- degPatterns(cluster_rlog.pTA, metadata = meta.pTA, minc = 100, time = "sample.conditions", col=NULL, eachStep = TRUE)
#save(clusters.all.1e25, 'clusters.all.1e25.Rdata')

                                        #rivanna

#cd /home/mjg7y/pro_pTA
#sbatch LRT_clustering_pTA.slurm
#LRT_largemem_clustering_pTA.R

load('/Volumes/GUERTIN_2/adipogenesis/pro_dREG/191230_clusters.pTA.pro.minc10.pval1e5.Rdata')

plot.df.pTA = clusters.all.1e5$normalized

plot.df.pTA$sample.conditions = as.character(plot.df.pTA$sample.conditions)
plot.df.pTA$sample.conditions[plot.df.pTA$sample.conditions == '20min'] = 20
plot.df.pTA$sample.conditions[plot.df.pTA$sample.conditions == '40min'] = 40
plot.df.pTA$sample.conditions[plot.df.pTA$sample.conditions == '60min'] = 60
plot.df.pTA$sample.conditions[plot.df.pTA$sample.conditions == '2h'] = 120
plot.df.pTA$sample.conditions[plot.df.pTA$sample.conditions == '3h'] = 180
plot.df.pTA$sample.conditions[plot.df.pTA$sample.conditions == '4h'] = 240
plot.df.pTA$sample.conditions = as.numeric(plot.df.pTA$sample.conditions)
plot.df.pTA = plot.df.pTA[order(plot.df.pTA$genes),]
plot.df.pTA = plot.df.pTA[order(plot.df.pTA$sample.conditions),]

plot.df.pTA$cluster = paste('cluster', as.character(plot.df.pTA$cluster), sep = '')




plot.df.pTA$chr = sapply(strsplit(plot.df.pTA$genes, ':'), '[', 1)
range.chr = sapply(strsplit(plot.df.pTA$genes, ':'), '[', 2)
plot.df.pTA$start = sapply(strsplit(range.chr, '-'), '[', 1)
plot.df.pTA$end = sapply(strsplit(range.chr, '-'), '[', 2)


for (i in unique(plot.df.pTA$cluster)) {
    print(i)
    write.table(plot.df.pTA[plot.df.pTA$cluster == i,
                        c('chr','start','end', 'value', 'cluster')][!duplicated(plot.df.pTA[plot.df.pTA$cluster == i,]$genes),],
                file = paste0('cluster_bed_pro_pTA',
                              gsub(" ", "", i, fixed = TRUE),'.bed'),
                quote = FALSE, row.names = FALSE, col.names = FALSE, sep = '\t')
}

library(data.table)


plot.df.pTA.cluster = dcast(plot.df.pTA, genes + cluster ~ sample.conditions, value.var="value")

avg.clusters = as.data.frame(matrix(nrow = 0, ncol = 7))
colnames(y) = colnames(plot.df.pTA.cluster[3:9])
for (i in unique(plot.df.pTA.cluster$cluster)) {
    z = data.frame(matrix(colMeans(plot.df.pTA.cluster[plot.df.pTA.cluster$cluster == i,3:9]), ncol = 7, nrow = 1))
    rownames(z) = c(i)
    colnames(z) = as.character(colnames(plot.df.pTA.cluster)[3:9])
    avg.clusters = rbind(avg.clusters, z)
}


dd = dist(avg.clusters)
hc = hclust(dd, method = "complete")


pdf(paste('dendrogram_HC_pTA','.pdf', sep=''), width=8, height=5)
plot(hc, xlab = "Clusters", main = ' ', hang = -1)
abline(h = 2.1, lty = 2)
dev.off()


plot.df.pTA.cluster[grep('Egr', plot.df.pTA.cluster[,1]),]





pdf(paste('Clusters_minc10_1e5_pvalue_pTA','.pdf', sep=''), width=11, height=32)

trellis.par.set(box.umbrella = list(lty = 1, col="black", lwd=1),
                box.rectangle = list( lwd=1.0, col="black", alpha = 1.0),
                plot.symbol = list(col="black", lwd=1.0, pch ='.'))
print(
    xyplot(value ~  sample.conditions | cluster, group = genes, data = plot.df.pTA, type = c('l'),#type = c('l','p'),
       scales=list(x=list(cex=1.0,relation = "free", rot = 45), y =list(cex=1.0, relation="free")),
       aspect=1.0,
#       layout = c(6,13),
       between=list(y=0.5, x=0.5),
#       index.cond=list(
#           c(10, 12, 24,
#             6, 39, 1,
#             15, 14,
#             27,26, 13, 30 , 23,
#             21,28, 19,
#             16,3,20,
#             22, 2, 9, 11,
#             33,
#             38, 34, 4,
#             25, 35, 17,
#             29, 5, 31, 7, 36, 8,
#             18,
#             37, 32)),
#       skip = c(F, F, F, T, T, T,
#                F, F, F, T, T, T,
#                F, F, T, T, T, T,
#                F, F, F, F, F, T,
#                F, F, F, T, T, T,
#                F, F, F, T, T, T,
#                F, F, F, F, T, T,
#                F, T, T, T, T, T,
#                F, F, F, T, T, T,
#                F, F, F, T, T, T ,
#                F, F, F, F, F, F,
#                F, T, T, T, T, T,
#               F, F, T, T, T, T),
       ylab = list(label = 'Normalized PRO signal', cex =1.0),
       xlab = list(label = 'Time (minutes)', cex =1.0),
       par.settings = list(superpose.symbol = list(pch = c(16),
                                                   col=c('grey20'), cex =0.5),
                           strip.background=list(col="grey80"),
                           superpose.line = list(col = c('#99999980'), lwd=c(1),
                                                 lty = c(1))),
       panel = function(x, y, ...) {
           panel.xyplot(x, y, ...)
           panel.bwplot(x, y, pch = '|', horizontal = FALSE, box.width = 15, do.out = FALSE)
           panel.loess(x, y, ..., col = "blue", lwd =2.0,  span = 1/2, degree = 1, family = c("gaussian"))
           
})

      )
dev.off()








                                        #currently doing de novo for each sub cluster and combined clusters.

                                        #dREG motif analysis

get.raw.counts.dREG <- function(df.prefix, path.to.bigWig, file.prefix = 'G') {
    vec.names = c()
    df.plus = read.table(paste0(df.prefix, '.plus.bed'))
    #dREG gives summits outside range for some reason
    df.plus[,2] = df.plus[,2] - 11
    df.plus[,3] = df.plus[,3] + 100
    df.plus[,2][df.plus[,2] < 0] = 1
    df.minus = read.table(paste0(df.prefix, '.minus.bed'))
    df.minus[,3] = df.minus[,3] + 11
    df.minus[,2] = df.minus[,2] -100
    df.minus[,2][df.minus[,2] < 0] = 1
    print('plus')
    print(head(df.plus))
    print('minus')
    print(head(df.minus))
    inten.df=data.frame(matrix(ncol = 0, nrow = nrow(df.plus)))
    for (mod.bigWig in Sys.glob(file.path(path.to.bigWig, 
                       paste(file.prefix, "*pro_plus.bigWig", sep ='')))) {
        factor.name = strsplit(strsplit(mod.bigWig, 
                        "/")[[1]][length(strsplit(mod.bigWig, "/")[[1]])], '_plus')[[1]][1]
        print(factor.name)
        vec.names = c(vec.names, factor.name)
        loaded.bw.plus = load.bigWig(mod.bigWig)
        print(mod.bigWig)
        print(paste(path.to.bigWig,'/',factor.name, '_minus.bigWig', sep=''))
        loaded.bw.minus = load.bigWig(paste(path.to.bigWig,'/',factor.name, 
                        '_minus.bigWig', sep=''))
        mod.inten.plus = bed.region.bpQuery.bigWig(loaded.bw.plus, df.plus)
        mod.inten.minus = bed.region.bpQuery.bigWig(loaded.bw.minus, df.minus)
        inten.df = cbind(inten.df, mod.inten.plus + mod.inten.minus)
    }
    colnames(inten.df) = vec.names
    r.names = paste(df.plus[,1], ':', df.plus[,2] + 11, '-', df.plus[,2] + 12, sep='')
    row.names(inten.df) = r.names
    return(inten.df)
}

setwd('/Volumes/GUERTIN_2/adipogenesis/pro_dREG')

direc = '/Volumes/GUERTIN_2/adipogenesis/pro_dREG'

counts.dreg = get.raw.counts.dREG('3T3_adipogenesis.dREG.peak', 
                                  direc, file.prefix = 'G')



colnames(counts.dreg) = sapply(strsplit(colnames(counts.dreg), '3T3_t'), '[', 2)
colnames(counts.dreg) = sapply(strsplit(colnames(counts.dreg), '_pro'), '[', 1)

save(counts.dreg, file = 'adipogenesis.counts.dreg.Rdata')
counts.dreg = counts.dreg[,1:21]
                                      #follow this: https://hbctraining.github.io/DGE_workshop/lessons/08_DGE_LRT.html
sample.conditions = factor(sapply(strsplit(as.character(colnames(counts.dreg)), '_'), '[', 1))

sample.conditions = factor(sample.conditions, levels=c("0","20min","40min","60min","2h","3h","4h"))

deseq.counts.table = DESeqDataSetFromMatrix(counts.dreg, as.data.frame(sample.conditions), ~ sample.conditions)

dds = DESeq(deseq.counts.table)

                                        #exploratory
rld = rlog(dds, blind=TRUE)

#
x = plotPCA(rld, intgroup="sample.conditions", returnData=TRUE)
plotPCAlattice(x, file = 'PCA_adipogenesis_pro_lattice_guertin.pdf')

                                        #lrt

dds.lrt = DESeq(dds, test="LRT", reduced = ~ 1)

#use for outside R
normalized.counts = counts(dds, normalized=TRUE)

res.lrt = results(dds.lrt)

padj.cutoff = 0.00001 #1e-5

siglrt.re = res.lrt[res.lrt$padj < padj.cutoff & !is.na(res.lrt$padj),]

dim(siglrt.re)                                    

rld_mat <- assay(rld)
cluster_rlog = rld_mat[rownames(siglrt.re),]
meta = as.data.frame(sample.conditions)
rownames(meta) = colnames(cluster_rlog)

save(cluster_rlog, meta, sample.conditions, file = 'pro_cluster_rlog_pval_0.00001.Rdata')

                                        #rivanna

load('191130_clusters.all.pro.minc100.pval0.00001.Rdata')





plot.df.pro = clusters.all.test.0.00001$normalized

plot.df.pro$sample.conditions = as.character(plot.df.pro$sample.conditions)
plot.df.pro$sample.conditions[plot.df.pro$sample.conditions == '20min'] = 20
plot.df.pro$sample.conditions[plot.df.pro$sample.conditions == '40min'] = 40
plot.df.pro$sample.conditions[plot.df.pro$sample.conditions == '60min'] = 60
plot.df.pro$sample.conditions[plot.df.pro$sample.conditions == '2h'] = 120
plot.df.pro$sample.conditions[plot.df.pro$sample.conditions == '3h'] = 180
plot.df.pro$sample.conditions[plot.df.pro$sample.conditions == '4h'] = 240
plot.df.pro$sample.conditions = as.numeric(plot.df.pro$sample.conditions)
plot.df.pro = plot.df.pro[order(plot.df.pro$genes),]
plot.df.pro = plot.df.pro[order(plot.df.pro$sample.conditions),]

plot.df.pro$cluster = paste('cluster', as.character(plot.df.pro$cluster), sep = '')




plot.df.pro$chr = sapply(strsplit(plot.df.pro$genes, ':'), '[', 1)
range.chr = sapply(strsplit(plot.df.pro$genes, ':'), '[', 2)
plot.df.pro$start = sapply(strsplit(range.chr, '-'), '[', 1)
plot.df.pro$end = sapply(strsplit(range.chr, '-'), '[', 2)


for (i in unique(plot.df.pro$cluster)) {
    print(i)
    write.table(plot.df.pro[plot.df.pro$cluster == i,
                        c('chr','start','end', 'value', 'cluster')][!duplicated(plot.df.pro[plot.df.pro$cluster == i,]$genes),],
                file = paste0('cluster_bed_pro_',
                              gsub(" ", "", i, fixed = TRUE),'.bed'),
                quote = FALSE, row.names = FALSE, col.names = FALSE, sep = '\t')
}

library(data.table)

plot.df.pro.cluster = dcast(plot.df.pro, genes + cluster ~ sample.conditions, value.var="value")

avg.clusters = as.data.frame(matrix(nrow = 0, ncol = 7))
colnames(y) = colnames(plot.df.pro.cluster[3:9])
for (i in unique(plot.df.pro.cluster$cluster)) {
    z = data.frame(matrix(colMeans(plot.df.pro.cluster[plot.df.pro.cluster$cluster == i,3:9]), ncol = 7, nrow = 1))
    rownames(z) = c(i)
    colnames(z) = as.character(colnames(plot.df.pro.cluster)[3:9])
    avg.clusters = rbind(avg.clusters, z)
}


dd = dist(avg.clusters)
hc = hclust(dd, method = "complete")


pdf(paste('dendrogram_HC_pro','.pdf', sep=''), width=8, height=5)
plot(hc, xlab = "Clusters", main = ' ', hang = -1)
abline(h = 1.8, lty = 2)
dev.off()


pdf(paste('Clusters_minc100_stringent_pvalue_post_HC_thirteen','.pdf', sep=''), width=11, height=32)

trellis.par.set(box.umbrella = list(lty = 1, col="black", lwd=1),
                box.rectangle = list( lwd=1.0, col="black", alpha = 1.0),
                plot.symbol = list(col="black", lwd=1.0, pch ='.'))
print(
    xyplot(value ~  sample.conditions | cluster, group = genes, data = plot.df.pro, type = c('l'),#type = c('l','p'),
       scales=list(x=list(cex=1.0,relation = "free", rot = 45), y =list(cex=1.0, relation="free")),
       aspect=1.0,
       layout = c(6,13),
       between=list(y=0.5, x=0.5),
       index.cond=list(
           c(10, 12, 24,
             6, 39, 1,
             15, 14,
             27,26, 13, 30 , 23,
             21,28, 19,
             16,3,20,
             22, 2, 9, 11,
             33,
             38, 34, 4,
             25, 35, 17,
             29, 5, 31, 7, 36, 8,
             18,
             37, 32)),
       skip = c(F, F, F, T, T, T,
                F, F, F, T, T, T,
                F, F, T, T, T, T,
                F, F, F, F, F, T,
                F, F, F, T, T, T,
                F, F, F, T, T, T,
                F, F, F, F, T, T,
                F, T, T, T, T, T,
                F, F, F, T, T, T,
                F, F, F, T, T, T ,
                F, F, F, F, F, F,
                F, T, T, T, T, T,
                F, F, T, T, T, T),

       ylab = list(label = 'Normalized PRO signal', cex =1.0),
       xlab = list(label = 'Time (minutes)', cex =1.0),
       par.settings = list(superpose.symbol = list(pch = c(16),
                                                   col=c('grey20'), cex =0.5),
                           strip.background=list(col="grey80"),
                           superpose.line = list(col = c('#99999980'), lwd=c(1),
                                                 lty = c(1))),
       panel = function(x, y, ...) {
           panel.xyplot(x, y, ...)
           panel.bwplot(x, y, pch = '|', horizontal = FALSE, box.width = 15, do.out = FALSE)
           panel.loess(x, y, ..., col = "blue", lwd =2.0,  span = 1/2, degree = 1, family = c("gaussian"))
           
})

      )
dev.off()

#cloose number of clusters for de novo

                                        #useful for network inference section 5.6 two wave network
                                        #https://kateto.net/netscix2016.html
#layered graph bipartite
#made a cyclical https://rdrr.io/cran/igraph/man/layout_with_sugiyama.html
#term sugiyama 
                                        #https://stackoverflow.com/questions/51194653/r-is-not-taking-the-parameter-hgap-in-layout-with-sugiyama
load.bigWig('/Volumes/GUERTIN_2/adipogenesis/pro_dREG/GSM3900969_3T3_t6d_rep2_pro_plus.bigWig')
coverage = o$basesCovered*o$mean
