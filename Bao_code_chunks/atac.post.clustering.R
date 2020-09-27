library(lattice)
library(data.table)

source('plot.traces.R')

load('clusters.all.minc100.1e8.Rdata')

#generate 'plot.df' object
plot.df = clusters.all.test.1e8$normalized

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

plot.df$chr = sapply(strsplit(plot.df$genes, '[.]'), '[', 1)
plot.df$start = sapply(strsplit(plot.df$genes, '[.]'), '[', 2)
plot.df$end = sapply(strsplit(plot.df$genes, '[.]'), '[', 3)

save(plot.df,file='plot.df.Rdata')

write.table(cbind(plot.df$chr, plot.df$start, plot.df$end), file = 'dynamic_peaks.bed', quote = FALSE, sep = '\t', col.names=FALSE, row.names=FALSE)

#plot all clusters 
for (i in unique(plot.df$cluster)) {
    print(i)
    write.table(plot.df[plot.df$cluster == i,
                        c('chr','start','end', 'value', 'cluster')][!duplicated(plot.df[plot.df$cluster == i,]$genes),],
                file = paste0('cluster_bed_',
                              gsub(" ", "", i, fixed = TRUE),'.bed'),
                quote = FALSE, row.names = FALSE, col.names = FALSE, sep = '\t')
}

pdf('atac_clusters.pdf', width=11, height=15)

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
           panel.spline(x, y, col = 'blue', lwd =2.0, ...)            
})

      )
dev.off()

#up.flat - 9,23,17,11
#grad.up - 5,8,10,1
#up.down - 6,2,12
#grad.down - 7,3,4,13,14
#down.up - 24

#multi-colored loess plot (not very helpful)
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

colors.clusters = rep('white', 17)
colors.clusters[c(9,23,17,11)] = 'red'
colors.clusters[c(5,8,10,1)] = 'blue'
colors.clusters[c(6,2,12)] = 'green'
colors.clusters[c(7,3,4,13,14)] = 'black'
colors.clusters[c(24)] = 'purple'
colors.clusters[c(14)] = 'orange'

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

#dendrogram
x = as.data.table(plot.df)
plot.df.cluster = dcast(x, genes + cluster ~ sample.conditions, value.var="value")

avg.clusters = as.data.frame(matrix(nrow = 0, ncol = 7))
for (i in unique(plot.df.cluster$cluster)) {
    z = data.frame(matrix(colMeans(plot.df.cluster[plot.df.cluster$cluster == i,3:9]), ncol = 7, nrow = 1))
    rownames(z) = c(i)
    colnames(z) = as.character(colnames(plot.df.cluster)[3:9])
    avg.clusters = rbind(avg.clusters, z)
}


dd = dist(avg.clusters)
hc = hclust(dd, method = "complete")

pdf('dendrogram.pdf', width=8, height=5)
plot(hc, xlab = "Clusters", main = ' ', hang = -1)
abline(h = 2, lty = 2)
dev.off()

pdf(paste('dendrogram_HC_eight','.pdf', sep=''), width=8, height=5)
plot(hc, xlab = "Clusters", main = ' ', hang = -1)
abline(h = 1.38, lty = 2)
dev.off()

#plot clusters again this time organized by supercluster

df = data.frame(index=1:17,cluster=unique(plot.df$cluster)[order(unique(plot.df$cluster))])
df$cluster.num = as.integer(sapply(strsplit(df$cluster, 'cluster'), '[', 2))
df = df[order(df$cluster.num),]
df = df[reorder(df$cluster.num,c(23,17,9,11,5,8,1,10,12,2,6,24,14,13,7,4,3)),]

pdf('atac_clusters_org_by_sc.pdf', width=11, height=15)

trellis.par.set(box.umbrella = list(lty = 1, col="black", lwd=1),
                box.rectangle = list( lwd=1.0, col="black", alpha = 1.0),
                plot.symbol = list(col="black", lwd=1.0, pch ='.'))
print(
    xyplot(value ~  sample.conditions | cluster, group = genes, data = plot.df, type = c('l'),#type = c('l','p'),
       scales=list(x=list(cex=1.0,relation = "free", rot = 45), y =list(cex=1.0, relation="free")),
       aspect=1.0,
       layout = c(5,5),
       between=list(y=0.5, x=0.5),
       index.cond=list(rev(df$index)),
       skip = c(F,F,F,F,F,
                F,T,T,T,T,
                F,F,F,T,T,
                F,F,F,F,T,
                F,F,F,F,T),
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
           panel.spline(x, y, col = 'blue', lwd =2.0, ...)            
})

      )
dev.off()

#plot supercluster traces
gradual.down = plot.df[plot.df$cluster == 'cluster7' |
                    plot.df$cluster == 'cluster3' |
                    plot.df$cluster == 'cluster4' |
                    plot.df$cluster == 'cluster13' |
                    plot.df$cluster == 'cluster14',]

down.up = plot.df[plot.df$cluster == 'cluster24',]

gradual.up = plot.df[plot.df$cluster == 'cluster5' |
                    plot.df$cluster == 'cluster8' |
                    plot.df$cluster == 'cluster10' |
                    plot.df$cluster == 'cluster1',]

up.flat = plot.df[plot.df$cluster == 'cluster9' |
                    plot.df$cluster == 'cluster23' |
                    plot.df$cluster == 'cluster17' |
                    plot.df$cluster == 'cluster11',]

up.down = plot.df[plot.df$cluster == 'cluster6' |
                    plot.df$cluster == 'cluster12' |
                    plot.df$cluster == 'cluster2',]

gradual.down$supercluster = 'gradual.down'
down.up$supercluster = 'down.up'
gradual.up$supercluster = 'gradual.up'
up.flat$supercluster = 'up.flat'
up.down$supercluster = 'up.down'

nrow(gradual.down)/7
nrow(down.up)/7
nrow(gradual.up)/7                       
nrow(up.flat)/7
nrow(up.down)/7

plot.df.atac = rbind(gradual.down,
      down.up,
      gradual.up,
      up.flat,
      up.down)
plot.df.atac = plot.df.atac[,-(7:26)]
plot.df.atac$genes = paste0(plot.df.atac$chr,':',plot.df.atac$start,'-',plot.df.atac$end)
save(plot.df.atac,file='plot.df.atac.Rdata')

pdf('atac_superclusters.pdf', width=6.83, height=5)

trellis.par.set(box.umbrella = list(lty = 1, col="black", lwd=1),
                box.rectangle = list( lwd=1.0, col="black", alpha = 1.0),
                plot.symbol = list(col="black", lwd=1.0, pch ='.'))
print(
xyplot(value ~  sample.conditions | supercluster, group = genes, data = plot.df.atac, type = c('l'),#type = c('l','p'),
       scales=list(x=list(cex=1.0,relation = "free", rot = 45), y =list(cex=1.0, relation="free")),
       aspect=1.0,
       between=list(y=0.5, x=0.5),
       layout = c(5,1),
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
           panel.spline(x, y, col = 'blue', lwd =2.0, ...) 
           
})

      )
dev.off()

for (i in unique(plot.df.atac$supercluster)) {
    print(i)
    write.table(plot.df.atac[plot.df.atac$supercluster == i,
                              c('chr','start','end', 'value', 'supercluster')][!duplicated(plot.df.atac[plot.df.atac$supercluster == i,]$genes),],
                file = paste0(i,'.bed'),
                quote = FALSE, row.names = FALSE, col.names = FALSE, sep = '\t')
}

####

df = data.frame(sc = c("down.up","gradual.down","gradual.up","up.down","up.flat"),
                col = c('orange1','#4169E1','#FF0000','green2','#A020F0'))

for(sc in unique(plot.df.atac$supercluster)) {
    col = df[df$sc == sc,]$col

    y = plot.df.atac[plot.df.atac$supercluster == sc,]
    
    pdf(file=paste0(sc,'.traces.pdf'))
    
    trellis.par.set(box.umbrella = list(lty = 1, col="black", lwd=1),
                    box.rectangle = list( lwd=1.0, col="black", alpha = 1.0),
                    plot.symbol = list(col="black", lwd=1.0, pch ='.'))
    
    print(
        xyplot(value ~  sample.conditions, group = genes, data = y, type = c('l'),
               scales=list(x=list(cex=1.0,relation = "free", rot = 45), y = list(cex=1.0, relation="free")),
               aspect=1.0,
               between=list(y=0.5, x=0.5),
               main = list(label = paste0(sc, ' traces'), cex = 1.5),
               ylab = list(label = 'Normalized ATAC signal', cex =1.0),
               xlab = list(label = 'Time (minutes)', cex =1.0),
               par.settings = list(superpose.symbol = list(pch = c(16),col=c('grey20'), cex =0.5),
                                   strip.background=list(col="grey80"),
                                   superpose.line = list(col = c('#99999980'), lwd=c(1),lty = c(1))),
               panel = function(x, y, ...) {
                   panel.xyplot(x, y, ...)
                   panel.bwplot(x, y, pch = '|', horizontal = FALSE, box.width = 15, do.out = FALSE)
                   panel.spline(x, y, col = 'blue', lwd = 3.5, ...)
                   #replace col = 'blue' with col = col for different sc colors
               })
    )
        
    dev.off()
    
}
