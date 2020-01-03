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

pdf(paste('Clusters_GRE_de_novo','.pdf', sep=''), width=3.5, height=3.5)

trellis.par.set(box.umbrella = list(lty = 1, col="black", lwd=1),
                box.rectangle = list( lwd=1.0, col="black", alpha = 1.0),
                plot.symbol = list(col="black", lwd=1.0, pch ='.'))
print(
    xyplot(value ~  sample.conditions , group = genes, data =
                       plot.df[plot.df$cluster == 'cluster5' | 
    plot.df$cluster == 'cluster10' | 
    plot.df$cluster == 'cluster2' | 
    plot.df$cluster == 'cluster21' | 
    plot.df$cluster == 'cluster16' | 
    plot.df$cluster == 'cluster17' | 
    plot.df$cluster == 'cluster13' | 
   plot.df$cluster == 'cluster15',]
         , type = c('l'),#type = c('l','p'),
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


pdf(paste('Clusters_TWIST_de_novo','.pdf', sep=''), width=3.5, height=3.5)

trellis.par.set(box.umbrella = list(lty = 1, col="black", lwd=1),
                box.rectangle = list( lwd=1.0, col="black", alpha = 1.0),
                plot.symbol = list(col="black", lwd=1.0, pch ='.'))
print(
    xyplot(value ~  sample.conditions , group = genes, data =
                       plot.df[plot.df$cluster == 'cluster5' | 
    plot.df$cluster == 'cluster3' | 
    plot.df$cluster == 'cluster4' | 
    plot.df$cluster == 'cluster7' | 
    plot.df$cluster == 'cluster11',]
         , type = c('l'),#type = c('l','p'),
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
normalized.counts.pro = counts(dds, normalized=TRUE)

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

#I am losing the strand somewhere I need to to maintain this

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
plot.df.pTA$end = sapply(strsplit(sapply(strsplit(plot.df.pTA$genes, '-'), '[', 2), '_'), '[', 1) #x=
plot.df.pTA$gene = sapply(strsplit(sapply(strsplit(plot.df.pTA$genes, '-'), '[', 2), '_'), '[', 2) #x=


for (i in unique(plot.df.pTA$cluster)) {
    print(i)
    write.table(plot.df.pTA[plot.df.pTA$cluster == i,
                        c('chr','start','end', 'gene', 'cluster')][!duplicated(plot.df.pTA[plot.df.pTA$cluster == i,]$genes),],
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


                                        #plot the raw PRO signal as opposed to standard deviations
format.plot = plot.df.pTA[grep('Twist2',plot.df.pTA$genes),]

twist2.pro = normalized.counts.pro[grep('Twist2', rownames(normalized.counts.pro)),]
format.plot[1,'value'] = mean(twist2.pro[1:3])
format.plot[2,'value'] = mean(twist2.pro[4:6])
format.plot[5,'value'] = mean(twist2.pro[7:9])
format.plot[6,'value'] = mean(twist2.pro[10:12])
format.plot[3,'value'] = mean(twist2.pro[13:15])
format.plot[7,'value'] = mean(twist2.pro[16:18])
format.plot[4,'value'] = mean(twist2.pro[19:21])


pdf(paste('Clusters_Twist2','.pdf', sep=''), width=3.5, height=3.5)

trellis.par.set(box.umbrella = list(lty = 1, col="black", lwd=1),
                box.rectangle = list( lwd=1.0, col="black", alpha = 1.0),
                plot.symbol = list(col="black", lwd=1.0, pch ='.'))
print(
    xyplot(value ~  sample.conditions , group = genes, data = plot.df.pTA[grep('Twist2',plot.df.pTA$genes),], type = c('l'),#type = c('l','p'),
       scales=list(x=list(cex=1.0,relation = "free", rot = 45), y =list(cex=1.0, relation="free")),
       aspect=1.0,
       between=list(y=0.5, x=0.5),
       ylab = list(label = 'Normalized PRO signal', cex =1.0),
       xlab = list(label = 'Time (minutes)', cex =1.0),
       par.settings = list(superpose.symbol = list(pch = c(16),
                                                   col=c('grey20'), cex =2.5),
                           strip.background=list(col="grey80"),
                           superpose.line = list(col = c('#99999980'), lwd=c(2),
                                                 lty = c(1))),
       panel = function(x, y, ...) {
           panel.xyplot(x, y, ...)
#           panel.bwplot(x, y, pch = '|', horizontal = FALSE, box.width = 15, do.out = FALSE)
           #panel.loess(x, y, ..., col = "blue", lwd =2.0,  span = 1/2, degree = 1, family = c("gaussian"))
           
})
      )
dev.off()



pdf(paste('Clusters_Twist2_raw_PRO','.pdf', sep=''), width=3.5, height=3.5)

trellis.par.set(box.umbrella = list(lty = 1, col="black", lwd=1),
                box.rectangle = list( lwd=1.0, col="black", alpha = 1.0),
                plot.symbol = list(col="black", lwd=1.0, pch ='.'))
print(
    xyplot(value ~  sample.conditions , group = genes, data = format.plot, type = c('l'),#type = c('l','p'),
       scales=list(x=list(cex=1.0,relation = "free", rot = 45), y =list(cex=1.0, relation="free")),
       aspect=1.0,
       between=list(y=0.5, x=0.5),
       ylab = list(label = 'Normalized PRO signal', cex =1.0),
       xlab = list(label = 'Time (minutes)', cex =1.0),
       par.settings = list(superpose.symbol = list(pch = c(16),
                                                   col=c('grey20'), cex =2.5),
                           strip.background=list(col="grey80"),
                           superpose.line = list(col = c('#99999980'), lwd=c(2),
                                                 lty = c(1))),
       panel = function(x, y, ...) {
           panel.xyplot(x, y, ...)
#           panel.bwplot(x, y, pch = '|', horizontal = FALSE, box.width = 15, do.out = FALSE)
           #panel.loess(x, y, ..., col = "blue", lwd =2.0,  span = 1/2, degree = 1, family = c("gaussian"))
           
})
      )
dev.off()


                                        #need to change normailzed count name to normal.ized.counts.pro

format.plot.Smyd3 = plot.df.pTA[grep("Smyd3",plot.df.pTA$genes),]

Smyd3.pro = normalized.counts.pro[grep("Smyd3", rownames(normalized.counts.pro)),]
format.plot.Smyd3[1,'value'] = mean(Smyd3.pro[1:3])
format.plot.Smyd3[2,'value'] = mean(Smyd3.pro[4:6])
format.plot.Smyd3[5,'value'] = mean(Smyd3.pro[7:9])
format.plot.Smyd3[6,'value'] = mean(Smyd3.pro[10:12])
format.plot.Smyd3[3,'value'] = mean(Smyd3.pro[13:15])
format.plot.Smyd3[7,'value'] = mean(Smyd3.pro[16:18])
format.plot.Smyd3[4,'value'] = mean(Smyd3.pro[19:21])



pdf(paste('Clusters_Smyd3_raw_PRO','.pdf', sep=''), width=3.5, height=3.5)

trellis.par.set(box.umbrella = list(lty = 1, col="black", lwd=1),
                box.rectangle = list( lwd=1.0, col="black", alpha = 1.0),
                plot.symbol = list(col="black", lwd=1.0, pch ='.'))
print(
    xyplot(value ~  sample.conditions , group = genes, data = format.plot.Smyd3, type = c('l'),#type = c('l','p'),
       scales=list(x=list(cex=1.0,relation = "free", rot = 45), y =list(cex=1.0, relation="free")),
       aspect=1.0,
       between=list(y=0.5, x=0.5),
       ylab = list(label = 'Normalized PRO signal', cex =1.0),
       xlab = list(label = 'Time (minutes)', cex =1.0),
       par.settings = list(superpose.symbol = list(pch = c(16),
                                                   col=c('grey20'), cex =2.5),
                           strip.background=list(col="grey80"),
                           superpose.line = list(col = c('#4c4c4c80'), lwd=c(8),
                                                 lty = c(1))),
       panel = function(x, y, ...) {
           panel.xyplot(x, y, ...)
#           panel.bwplot(x, y, pch = '|', horizontal = FALSE, box.width = 15, do.out = FALSE)
           #panel.loess(x, y, ..., col = "blue", lwd =2.0,  span = 1/2, degree = 1, family = c("gaussian"))
           
})
      )
dev.off()







                                        #plot GRE
load('../atac/191119_clusters.all.minc100.pval0.00000001.Rdata')

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


pdf(paste('GRE_near_twist','.pdf', sep=''), width=3.5, height=3.5)

trellis.par.set(box.umbrella = list(lty = 1, col="black", lwd=1),
                box.rectangle = list( lwd=1.0, col="black", alpha = 1.0),
                plot.symbol = list(col="black", lwd=1.0, pch ='.'))
print(
xyplot(value ~  sample.conditions , group = genes, data = plot.df[grep('chr1:91927918',plot.df$genes),], type = c('l'),#type = c('l','p'),
       scales=list(x=list(cex=1.0,relation = "free", rot = 45), y =list(cex=1.0, relation="free")),
       aspect=1.0,
       between=list(y=0.5, x=0.5),
       main= 'chr1:91927918-91928118',
       cex.main = 0.2,
       ylab = list(label = 'Normalized ATAC signal', cex =1.0),
       xlab = list(label = 'Time (minutes)', cex =1.0),
       par.settings = list(superpose.symbol = list(pch = c(16),
                                                   col=c('grey20'), cex =25),
                           strip.background=list(col="grey80"),
                           superpose.line = list(col = c('#99999980'), lwd=c(2),
                                                 lty = c(1))),
       panel = function(x, y, ...) {
           panel.xyplot(x, y, ...)
#           panel.bwplot(x, y, pch = '|', horizontal = FALSE, box.width = 15, do.out = FALSE)
#           panel.loess(x, y, ..., col = "blue", lwd =2.0,  span = 1/2, degree = 1, family = c("gaussian"))
           
})

      )
dev.off()

format.plot.atac =plot.df[grep('chr1:179023571-179023771',plot.df$genes),]

Smyd3_1.atac = normalized.counts[grep('chr1:179023571-179023771', rownames(normalized.counts)),]
format.plot.atac[2,'value'] = mean(Smyd3_1.atac[1:3])
format.plot.atac[5,'value'] = mean(Smyd3_1.atac[4:6])
format.plot.atac[6,'value'] = mean(Smyd3_1.atac[7:9])
format.plot.atac[3,'value'] = mean(Smyd3_1.atac[10:12])
format.plot.atac[7,'value'] = mean(Smyd3_1.atac[13:15])
format.plot.atac[4,'value'] = mean(Smyd3_1.atac[16:18])
format.plot.atac[1,'value'] = mean(Smyd3_1.atac[19:21])


pdf(paste('TWIST_near_Smyd3_1','.pdf', sep=''), width=3.5, height=3.5)

trellis.par.set(box.umbrella = list(lty = 1, col="black", lwd=1),
                box.rectangle = list( lwd=1.0, col="black", alpha = 1.0),
                plot.symbol = list(col="black", lwd=1.0, pch ='.'))
print(
xyplot(value ~  sample.conditions , group = genes, data = format.plot.atac, type = c('l'),#type = c('l','p'),
       scales=list(x=list(cex=1.0,relation = "free", rot = 45), y =list(cex=1.0, relation="free")),
       aspect=1.0,
       between=list(y=0.5, x=0.5),
       main= 'chr1:179023571-179023771',
       cex.main = 0.2,
       ylab = list(label = 'Normalized ATAC signal', cex =1.0),
       xlab = list(label = 'Time (minutes)', cex =1.0),
       par.settings = list(superpose.symbol = list(pch = c(16),
                                                   col=c('grey20'), cex =25),
                           strip.background=list(col="grey80"),
                           superpose.line = list(col = c('#4c4c4c80'), lwd=c(8),
                                                 lty = c(1))),
       panel = function(x, y, ...) {
           panel.xyplot(x, y, ...)
#           panel.bwplot(x, y, pch = '|', horizontal = FALSE, box.width = 15, do.out = FALSE)
#           panel.loess(x, y, ..., col = "blue", lwd =2.0,  span = 1/2, degree = 1, family = c("gaussian"))
           
})

      )
dev.off()



format.plot.atac =plot.df[grep('chr1:179097150-179097350',plot.df$genes),]

Smyd3_2.atac = normalized.counts[grep('chr1:179097150-179097350', rownames(normalized.counts)),]
format.plot.atac[2,'value'] = mean(Smyd3_2.atac[1:3])
format.plot.atac[5,'value'] = mean(Smyd3_2.atac[4:6])
format.plot.atac[6,'value'] = mean(Smyd3_2.atac[7:9])
format.plot.atac[3,'value'] = mean(Smyd3_2.atac[10:12])
format.plot.atac[7,'value'] = mean(Smyd3_2.atac[13:15])
format.plot.atac[4,'value'] = mean(Smyd3_2.atac[16:18])
format.plot.atac[1,'value'] = mean(Smyd3_2.atac[19:21])


pdf(paste('TWIST_near_Smyd3_2','.pdf', sep=''), width=3.5, height=3.5)

trellis.par.set(box.umbrella = list(lty = 1, col="black", lwd=1),
                box.rectangle = list( lwd=1.0, col="black", alpha = 1.0),
                plot.symbol = list(col="black", lwd=1.0, pch ='.'))
print(
xyplot(value ~  sample.conditions , group = genes, data = format.plot.atac, type = c('l'),#type = c('l','p'),
       scales=list(x=list(cex=1.0,relation = "free", rot = 45), y =list(cex=1.0, relation="free")),
       aspect=1.0,
       between=list(y=0.5, x=0.5),
       main= 'chr1:179097150-179097350',
       cex.main = 0.2,
       ylab = list(label = 'Normalized ATAC signal', cex =1.0),
       xlab = list(label = 'Time (minutes)', cex =1.0),
       par.settings = list(superpose.symbol = list(pch = c(16),
                                                   col=c('grey20'), cex =25),
                           strip.background=list(col="grey80"),
                           superpose.line = list(col = c('#4c4c4c80'), lwd=c(8),
                                                 lty = c(1))),
       panel = function(x, y, ...) {
           panel.xyplot(x, y, ...)
#           panel.bwplot(x, y, pch = '|', horizontal = FALSE, box.width = 15, do.out = FALSE)
#           panel.loess(x, y, ..., col = "blue", lwd =2.0,  span = 1/2, degree = 1, family = c("gaussian"))
           
})

      )
dev.off()








load("../atac/191119_adipogenesis.Rdata")
format.plot.atac =plot.df[grep('chr1:91857151',plot.df$genes),]

twist2.atac = normalized.counts[grep('chr1:91857151', rownames(normalized.counts)),]
format.plot.atac[2,'value'] = mean(twist2.atac[1:3])
format.plot.atac[5,'value'] = mean(twist2.atac[4:6])
format.plot.atac[6,'value'] = mean(twist2.atac[7:9])
format.plot.atac[3,'value'] = mean(twist2.atac[10:12])
format.plot.atac[7,'value'] = mean(twist2.atac[13:15])
format.plot.atac[4,'value'] = mean(twist2.atac[16:18])
format.plot.atac[1,'value'] = mean(twist2.atac[19:21])


pdf(paste('GRE_site_near_twist2','.pdf', sep=''), width=3.5, height=3.5)

trellis.par.set(box.umbrella = list(lty = 1, col="black", lwd=1),
                box.rectangle = list( lwd=1.0, col="black", alpha = 1.0),
                plot.symbol = list(col="black", lwd=1.0, pch ='.'))
           print(
 xyplot(value ~  sample.conditions , group = genes, data = format.plot.atac, type = c('l'),#type = c('l','p'),              chr1:91875333-91875533
#               xyplot(value ~  sample.conditions , group = genes, data = plot.df[grep('chr1:91857151',plot.df$genes),], type = c('l'),#type = c('l','p'),
#xyplot(value ~  sample.conditions , group = genes, data = plot.df[grep('91927918-91928118',plot.df$genes),], type = c('l'),#type = c('l','p'),
       scales=list(x=list(cex=1.0,relation = "free", rot = 45), y =list(cex=1.0, relation="free")),
       aspect=1.0, 
       between=list(y=0.5, x=0.5),
       main= list(
           label='chr1:91857151-91857351',
           cex=0.75),
       ylab = list(label = 'Normalized ATAC signal', cex =1.0),
       xlab = list(label = 'Time (minutes)', cex =1.0),
       par.settings = list(superpose.symbol = list(pch = c(16),
                                                   col=c('grey50'), cex =25),
                           cex.main = 0.2,
                           strip.background=list(col="grey80"),
                           superpose.line = list(col = c('#4c4c4c80'), lwd=c(2),
                                                 lty = c(1))),
       panel = function(x, y, ...) {
           panel.xyplot(x, y, ...)
#           panel.bwplot(x, y, pch = '|', horizontal = FALSE, box.width = 15, do.out = FALSE)
#           panel.loess(x, y, ..., col = "blue", lwd =2.0,  span = 1/2, degree = 1, family = c("gaussian"))
           
})

      )
dev.off()



pdf(paste('Atf3_GRE_site_near_twist','.pdf', sep=''), width=3.5, height=3.5)

trellis.par.set(box.umbrella = list(lty = 1, col="black", lwd=1),
                box.rectangle = list( lwd=1.0, col="black", alpha = 1.0),
                plot.symbol = list(col="black", lwd=1.0, pch ='.'))
           print(
 xyplot(value ~  sample.conditions , group = genes, data = plot.df[grep('chr1:91857151',plot.df$genes),], type = c('l'),#type = c('l','p'),              chr1:91875333-91875533
#               xyplot(value ~  sample.conditions , group = genes, data = plot.df[grep('chr1:91857151',plot.df$genes),], type = c('l'),#type = c('l','p'),
#xyplot(value ~  sample.conditions , group = genes, data = plot.df[grep('91927918-91928118',plot.df$genes),], type = c('l'),#type = c('l','p'),
       scales=list(x=list(cex=1.0,relation = "free", rot = 45), y =list(cex=1.0, relation="free")),
       aspect=1.0, 
       between=list(y=0.5, x=0.5),
       main= list(
           label='chr1:91857151-91857351',
           cex=0.75),
       ylab = list(label = 'Normalized ATAC signal', cex =1.0),
       xlab = list(label = 'Time (minutes)', cex =1.0),
       par.settings = list(superpose.symbol = list(pch = c(16),
                                                   col=c('grey50'), cex =25),
                           cex.main = 0.2,
                           strip.background=list(col="grey80"),
                           superpose.line = list(col = c('#4c4c4c80'), lwd=c(2),
                                                 lty = c(1))),
       panel = function(x, y, ...) {
           panel.xyplot(x, y, ...)
#           panel.bwplot(x, y, pch = '|', horizontal = FALSE, box.width = 15, do.out = FALSE)
#           panel.loess(x, y, ..., col = "blue", lwd =2.0,  span = 1/2, degree = 1, family = c("gaussian"))
           
})

      )
dev.off()



                                        #get coposite ATAC at gres
source('https://raw.githubusercontent.com/guertinlab/seqOutBias/master/docs/R/seqOutBias_functions.R')

bed.window <- function(bed, half.window) {
    bed[,2] = (bed[,2] + bed[,3])/2 - half.window
    bed[,3] = bed[,2] + 2 * half.window
    return(bed)
}



get.ATAC.interval <- function(bed, tf, path.to.bigWig, stp = 20, half.win, file.suffix = '.bigWig') {
    all.fimo = data.frame(matrix(ncol = 4, nrow = 0))
    colnames(all.fimo) = c('density', 'tf', 'cond', 'range')
    bed.win = bed.window(bed, half.win)
    for (mod.bigWig in Sys.glob(file.path(path.to.bigWig,
                                          paste('*', file.suffix, sep ='')))) {
        cond = strsplit(strsplit(mod.bigWig, "/")[[1]][length(strsplit(mod.bigWig, '/')[[1]])], file.suffix)[[1]][1]
        print(cond)
        loaded.bw = load.bigWig(mod.bigWig)
        inten = bed.step.probeQuery.bigWig(loaded.bw, bed.win,
                                           gap.value = 0, step = stp, as.matrix = TRUE)
        query.df = data.frame(cbind(colMeans(inten), tf, cond,
                                 seq(-half.win+(0.5*stp), half.win-(0.5*stp), stp)), stringsAsFactors=F)
        colnames(query.df) = c('density', 'tf', 'cond', 'range')
        all.fimo = rbind(all.fimo, query.df)
    }
    all.fimo[,1] = as.numeric(all.fimo[,1])
    all.fimo[,4] = as.numeric(all.fimo[,4])
    return(all.fimo)
}

composites.func.panels.lattice <- function(dat, fact = 'ATAC', summit = 'Motif Summit', class= '', num.m = -400, num.p =400, y.low =0, y.high = 0.2,
                                           col.lines = c("#3B82AE", "#547294", "#6D617A", "#865160", "#9F4046", "#B8302C", "#D11F12", "#DE1705"),
                                           fill.poly = c(rgb(0,0,1,1/4), 
                                                                                                        rgb(1,0,0,1/4), rgb(0.1,0.5,0.05,1/4),rgb(0,0,0,1/4), rgb(1/2,0,1/2,1/4))) {
    pdf(paste('composite_', fact, '_signals_', summit, '_peaks', class, '.pdf', sep=''), width=8.83, 
        height=4) 
    print(xyplot(density ~ range|tf, group = cond, data = dat,
                 type = 'l',
               scales=list(x=list(cex=0.8,relation = "free"), y =list(cex=0.8, relation="free")),
                 xlim=c(num.m,num.p),
#                 ylim = c(y.low, y.high),
                 col = col.lines,
                 #main=list(label=class, cex=0.6),
                 auto.key = list(points=F, lines=T, cex=0.8),
               par.settings = list(superpose.symbol = list(pch = c(16), col=col.lines, cex =0.7), 
                   superpose.line = list(col = col.lines, lwd=c(2,2,2,2,2,2), 
                       lty = c(1,1,1,1,1,1,1,1,1))),
                 cex.axis=1.0,
                 par.strip.text=list(cex=0.9, font=1, col='black'),
                 aspect=1.0,
                 between=list(y=0.5, x=0.5),
#                 index.cond = list(c(3, 2, 1)),
                 lwd=2,
                 ylab = list(label = paste(fact," Signal", sep=''), cex =1),
                 xlab = list(label = paste("Distance from ", summit, sep=''), cex =1),
                 #upper = dat$upper,
                 #fill = fill.poly,
                 #lower = dat$lower,
                 strip = function(..., which.panel, bg) {
                     bg.col = c("grey90")#,"#ce228e" , "#2290cf","grey60")
                 strip.default(..., which.panel = which.panel, bg = rep(bg.col, length = which.panel)[which.panel])
                 },
                 #panel = function(x, y, ...){
                     #panel.superpose(x, y, panel.groups = 'my.panel.bands', ...)
                     #panel.xyplot(x, y, ...)
       #}
                 ))
    dev.off()
}

nr3c1.bed = read.table('/Volumes/GUERTIN_2/adipogenesis/atac/NR3C1_motifs_in_dynamic.bed')
twist.bed = read.table('/Volumes/GUERTIN_2/adipogenesis/atac/TWIST1_motifs_in_dynamic.bed')

path.bw = '/Volumes/GUERTIN_2/adipogenesis/atac/4hour/'



nrc31.comp = get.ATAC.interval(nr3c1.bed, tf = 'GR', path.bw, stp = 20, half.win = 500, file.suffix = '.merged.bigWig')
nrc31.comp$cond = factor(nrc31.comp$cond, levels = c("3T3_t0", "3T3_20min", "3T3_40min", "3T3_60min", "3T3_2hr", "3T3_3hr", "3T3_4hr"))


twist.comp = get.ATAC.interval(twist.bed, tf = 'TWIST', path.bw, stp = 20, half.win = 500, file.suffix = '.merged.bigWig')
twist.comp$cond = factor(nrc31.comp$cond, levels = c("3T3_t0", "3T3_20min", "3T3_40min", "3T3_60min", "3T3_2hr", "3T3_3hr", "3T3_4hr"))

#need to incorporate a bed file loop into the composite function
composites.tfs = rbind(nrc31.comp, twist.comp)

composites.func.panels.lattice(composites.tfs, num.m = -400, num.p =400)

                                        #GRE de novo clusters


pdf(paste('Clusters_with_GRE_de_novo','.pdf', sep=''), width=9, height=10)

trellis.par.set(box.umbrella = list(lty = 1, col="black", lwd=1),
                box.rectangle = list( lwd=1.0, col="black", alpha = 1.0),
                plot.symbol = list(col="black", lwd=1.0, pch ='.'))
print(
    xyplot(value ~  sample.conditions | cluster, group = genes, data = plot.df, type = c('l'),#type = c('l','p'),
       scales=list(x=list(cex=1.0,relation = "free", rot = 45), y =list(cex=1.0, relation="free")),
       aspect=1.0,
       layout = c(5,3),
       between=list(y=2.0, x=0.5),
       index.cond=list(
           c(4,
             12,8,9,5,7,
             17, 2, 11)),
       skip = c(F, T, T, T, T,
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



pdf(paste('Clusters_with_TWIST_de_novo','.pdf', sep=''), width=9, height=10)

trellis.par.set(box.umbrella = list(lty = 1, col="black", lwd=1),
                box.rectangle = list( lwd=1.0, col="black", alpha = 1.0),
                plot.symbol = list(col="black", lwd=1.0, pch ='.'))
print(
    xyplot(value ~  sample.conditions | cluster, group = genes, data = plot.df, type = c('l'),#type = c('l','p'),
       scales=list(x=list(cex=1.0,relation = "free", rot = 45), y =list(cex=1.0, relation="free")),
       aspect=1.0,
       layout = c(5,1),
       between=list(y=2.0, x=0.5),
      index.cond=list(
           c(
             19, 15,16,3
             )),
       skip = c(
                F, F, F, F, T
              ) ,
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




pdf(paste('Clusters_early_TWIST','.pdf', sep=''), width=3.5, height=3.5)

trellis.par.set(box.umbrella = list(lty = 1, col="black", lwd=1),
                box.rectangle = list( lwd=1.0, col="black", alpha = 1.0),
                plot.symbol = list(col="black", lwd=1.0, pch ='.'))
print(
    xyplot(value ~  sample.conditions, group = genes, data =
                plot.df.pTA[plot.df.pTA$cluster == 'cluster3' |
                           plot.df.pTA$cluster ==  'cluster13' |
                           plot.df.pTA$cluster ==  'cluster4' |
                          plot.df.pTA$cluster ==   'cluster25' |
                        plot.df.pTA$cluster ==     'cluster31' |
                          plot.df.pTA$cluster ==   'cluster23' |
                          plot.df.pTA$cluster ==   'cluster7' |
                          plot.df.pTA$cluster ==   'cluster15' |
                          plot.df.pTA$cluster ==   'cluster26' |
                           plot.df.pTA$cluster ==  'cluster42' ,]
         , type = c('l'),#type = c('l','p'),
       scales=list(x=list(cex=1.0,relation = "free", rot = 45), y =list(cex=1.0, relation="free")),
       aspect=1.0,
#       layout = c(6,13),
       between=list(y=0.5, x=0.5),
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




pdf(paste('Clusters_decrease_late_SMYD3','.pdf', sep=''), width=3.5, height=3.5)

trellis.par.set(box.umbrella = list(lty = 1, col="black", lwd=1),
                box.rectangle = list( lwd=1.0, col="black", alpha = 1.0),
                plot.symbol = list(col="black", lwd=1.0, pch ='.'))
print(
    xyplot(value ~  sample.conditions, group = genes, data =
                plot.df.pTA[plot.df.pTA$cluster == 'cluster9' ,]
         , type = c('l'),#type = c('l','p'),
       scales=list(x=list(cex=1.0,relation = "free", rot = 45), y =list(cex=1.0, relation="free")),
       aspect=1.0,
#       layout = c(6,13),
       between=list(y=0.5, x=0.5),
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














####working ith iGraph

nr3c1.twist2.node = read.table('/Volumes/GUERTIN_2/adipogenesis/atac/NR3C1_RE_to_Twist2_TU.bed')
re.nodes = paste(paste(nr3c1.twist2.node$V6, nr3c1.twist2.node$V7, sep =':'), nr3c1.twist2.node$V8, sep = '-')

first.cis.wave.edges = data.frame(cbind(re.nodes, as.character(nr3c1.twist2.node$V4)), stringsAsFactors =FALSE)

colnames(first.cis.wave.edges) = c('from', 'to')

#I do not carry over the Twist2 TU in this data frame, so I need to work on
#how to do this systematically. We can do this manually for now because
#the set of early modulated TFs that have de novo motifs identified in late response REs is probably very limited

twist2.re.node = read.table('/Volumes/GUERTIN_2/adipogenesis/atac/REs_w_TWIST_motifs_in_dynamic_ATAC_clusters.bed', stringsAsFactors =FALSE)
tu.nodes = paste(paste(twist2.re.node$V1, twist2.re.node$V2, sep =':'), twist2.re.node$V3, sep = '-')

first.trans.wave.edges = data.frame(cbind('Twist2', tu.nodes), stringsAsFactors =FALSE)
colnames(first.trans.wave.edges) = c('from', 'to')


                                        #next set of cis edges


twist2.tus.node = read.table('/Volumes/GUERTIN_2/adipogenesis/atac/decrease_late_TUs_near_REs_w_TWIST.bed')

#I reduced it here to make it visible on the graph
twist2.tus.node.cluster9 = twist2.tus.node[twist2.tus.node[,4] == 'cluster9',]

tu.nodes.2 = paste(paste(twist2.tus.node.cluster9$V6, twist2.tus.node.cluster9$V7, sep =':'), twist2.tus.node.cluster9$V8, sep = '-')

second.cis.wave.edges = data.frame(cbind(tu.nodes.2, as.character(twist2.tus.node$V4)), stringsAsFactors =FALSE)

colnames(second.cis.wave.edges) = c('from', 'to')



#full network
full.network = rbind(first.cis.wave.edges, first.trans.wave.edges, second.cis.wave.edges)

#this list is something I googled for.
#we want a complete one from TFClass
list.tfs = read.table('/Volumes/GUERTIN_2/adipogenesis/pro_dREG/list_tfs_putative.txt', sep = '\t')


tu.w.cis = second.cis.wave.edges[ toupper(second.cis.wave.edges$to) %in% toupper(list.tfs[,1]) ,]

#further reducing the $ nodes


#I want to prune RE leaf nodes
re.w.cis = first.trans.wave.edges[first.trans.wave.edges$to %in% tu.w.cis$from,]




#includes leaf nodes in final layer only
reduced.network = rbind(first.cis.wave.edges, re.w.cis, tu.w.cis )


library(igraph)
                                        #create the graph variable
list.tfs 
g=graph.data.frame(reduced.network, directed=TRUE)

vec.gene.names = read.table('/Volumes/GUERTIN_2/adipogenesis/primary_transcript_annotation.bed', stringsAsFactors =FALSE)[,4]
V(g)$type <- ifelse(V(g)$name %in% vec.gene.names, TRUE, FALSE)





V(g)$color <- ifelse(V(g)$type, "#73F992", "#FF5AF3")
V(g)$shape <- ifelse(V(g)$type, "square", "circle")
E(g)$color <- "gray30"
V(g)$size <- ifelse(V(g)$type, 16, 7)
V(g)$label <- ifelse(V(g)$type, V(g)$name, ' ')
V(g)$label.cex <- 0.75

x = layout_with_sugiyama(g, attributes="all", hgap = 2, vgap = 1,)

pdf(paste('Twist_network','.pdf', sep=''), width=10, height=10, useDingbats=FALSE)

plot(x$extd_graph,  vertex.label.family = "Courier", vertex.label.font= 2, vertex.label.color = "black", edge.arrow.size=1.0)

dev.off()



threecol=read.csv("fake_network_input.txt",
header=F,stringsAsFactors = F,sep='\t')
colnames(threecol)=c('from','to','e_value')
threecol$weight=abs(log(threecol$e_value))
library(igraph)
                                        #create the graph variable

g=graph.data.frame(threecol,directed=TRUE)

vec.gene.names = c('TU1', 'TU2','TU3','TU4')
V(g)$type <- ifelse(V(g)$name %in% vec.gene.names, TRUE, FALSE)





V(g)$color <- ifelse(V(g)$type, "lightblue", "salmon")
V(g)$shape <- ifelse(V(g)$type, "square", "circle")
E(g)$color <- "gray30"
V(g)$size <- ifelse(V(g)$type, 20, 10)
V(g)$label <- ifelse(V(g)$type, V(g)$name, ' ')
x = layout_with_sugiyama(g, attributes="all")



pdf(paste('Twist_network','.pdf', sep=''), width=5, height=5)

plot(x$extd_graph, vertex.label.family = "Courier", vertex.label.font= 2, vertex.label.color = "black")

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
