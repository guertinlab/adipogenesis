                                        #how I got: dynamic_peaks.bed

chrom = vector(mode="character", length = length(unique(plot.df.atac$genes)))
start = vector(mode="integer", length = length(unique(plot.df.atac$genes)))
end = vector(mode="integer", length = length(unique(plot.df.atac$genes)))

count = 0
for (j in unique(plot.df.atac$genes)) {
    count = count +1
    chrom[count] = strsplit(j, ":")[[1]][1]
    start[count] = as.numeric(strsplit(strsplit(j, ":")[[1]][2], "-")[[1]][1])
    end[count] = as.numeric(strsplit(strsplit(j, ":")[[1]][2], "-")[[1]][2])
}













#get fimo info
res = data.frame(cbind(chrom, start, end))
colnames(res) = c('chr', 'start', 'end')
res$start = as.numeric(as.character(res$start))
res$end = as.numeric(as.character(res$end))

write.table(res, file='dynamic_peaks.bed', sep='\t', quote=F, row.names=F, col.names=F)
rownames(res) = paste0(res[,1], ':', res[,2], '-', res[,3])


res = head(res, 40)

count = -1
vec.names = c()
for(bed.file in Sys.glob(file.path('/Users/guertinlab/Desktop/fimo_stuff/*family_fimo.bed'))) {
    print(bed.file)
    count = count + 1
    factor.name = strsplit(bed.file, "/")[[1]]
    factor.name = strsplit(factor.name[length(factor.name)],
                           '_fimo.bed')[[1]][1]
    print(factor.name)
    x = read.table(bed.file, stringsAsFactors=FALSE)
    x = x[x[,6] != -1,]
    y = aggregate(as.numeric(V8)~V1+V2+V3, data=x, FUN=sum)
    colnames(y) = c('chr', 'start', 'end', factor.name)
    rownames(y) = paste0(y[,1], ':', y[,2], '-', y[,3])
    print(head(y))
    res = merge(res, y, by="row.names", all = TRUE)
    rnames = res[,1]
    print(head(res))
    vec.names = c(vec.names, factor.name)
    res = data.frame(res[,-c((ncol(res) -3 ):(ncol(res) -1))])
    res = data.frame(res[,(ncol(res) - count):ncol(res)])
    colnames(res) = vec.names
    rownames(res) = rnames
    print(head(res))
}

fimo.scores.atac = res
save(fimo.scores.atac, file = 'fimo.scores.atac.Rdata')


sapply(fimo.scores.atac, function(x) sum(!is.na(x)))

barplot(sapply(fimo.scores.atac, function(x) sum(!is.na(x))))


#add this tho super cluster plot:
for (i in 1:nrow(plot.df.atac)) {
    re = plot.df.atac$genes[i] 
    plot.df.atac[i, 11:15] = fimo.scores.atac[re,]
}

plot.df.atac = x
save(plot.df.atac, file = 'plot.df.atac.fimo.Rdata')

plot.df.atac$supercluster <- as.factor(plot.df.atac$supercluster)
plot.df.atac$genes <- as.factor(plot.df.atac$genes)

cat.colours = plot.df.atac[plot.df.atac$merge == 'one_groupt0', c(1,10)]
cat.colours$genes <- as.factor(cat.colours$genes)
cat.colours$supercluster <- as.factor(cat.colours$supercluster)


cat.colours$colour[cat.colours$supercluster == 'up.flat'] <- '#FF000008'
cat.colours$colour[cat.colours$supercluster == 'up.down'] <- '#FF000008'
cat.colours$colour[cat.colours$supercluster == 'gradual.up'] <- '#FF000008'
cat.colours$colour[cat.colours$supercluster == 'gradual.down'] <- '#0000FF08'
cat.colours$colour[cat.colours$supercluster == 'down.up'] <- '#0000FF08'

cat.colours$colour <- as.factor(cat.colours$colour)

cat.colours <- cat.colours[match(levels(plot.df.atac$genes), cat.colours$genes), ]



pdf(file = paste0('supercluster.traces.transcol.pdf'),width=14,height=4)

trellis.par.set(box.umbrella = list(lty = 1, col="black", lwd=1),
                box.rectangle = list( lwd=1.0, col="black", alpha = 1.0),
                plot.symbol = list(col="black", lwd=1.0, pch ='.'))

print(
    xyplot(value ~  sample.conditions | supercluster, group = genes, data = plot.df.atac, type = c('l'),#type = c('l','p'),
           scales=list(x=list(cex=1.0,relation = "free", rot = 45), y = list(cex=1.0, relation="free")),
           aspect=1.0,
           layout = c(5,1),
           between=list(y=0.5, x=0.5),
           index.cond=list(c(5:1)),
           #main = list(label = 'Peaks Separated by Supercluster', cex = 1.5),
           ylab = list(label = 'Normalized ATAC signal', cex =1.0),
           xlab = list(label = 'Time (minutes)', cex =1.0),
           par.settings = list(superpose.symbol = list(pch = c(16),col=c('grey20'), cex =0.5),
                               strip.background=list(col="grey80"),
                               superpose.line = list(col = as.character(cat.colours$colour), lwd=c(1),lty = c(1))),
           panel = function(x, y, ...) {
               panel.xyplot(x, y, ...)
               #panel.violin(x, y, horizontal = FALSE, col = "transparent",
               #         box.width = 20, do.out = FALSE)
               panel.bwplot(x, y, pch = '|', horizontal = FALSE, box.width = 10, do.out = FALSE)
               panel.spline(x, y, col = 'grey70', lwd =3.0, ...) 
               #panel.loess(x, y, ..., col = "blue", lwd =2.0,  span = 1/2, degree = 1, family = c("gaussian"))
           })
)

dev.off()



plot.df.atac$AP1 <- ifelse(is.na(plot.df.atac$AP1_family), FALSE, TRUE)
plot.df.atac$GR <- ifelse(is.na(plot.df.atac$GR_family), FALSE, TRUE)
plot.df.atac$SP <- ifelse(is.na(plot.df.atac$SP_family), FALSE, TRUE)
plot.df.atac$TWIST <- ifelse(is.na(plot.df.atac$TWIST_family), FALSE, TRUE)
plot.df.atac$bHLH <- ifelse(is.na(plot.df.atac$bHLH_family), FALSE, TRUE)




ap1 = plot.df.atac[plot.df.atac$AP1 == TRUE,]
ap1$motif = 'AP1_family'
gr = plot.df.atac[plot.df.atac$GR == TRUE,]
gr$motif = 'GR_family'
sp = plot.df.atac[plot.df.atac$SP == TRUE,]
sp$motif = 'SP_family'
twist = plot.df.atac[plot.df.atac$TWIST == TRUE,]
twist$motif = 'TWIST_family'
bhlh = plot.df.atac[plot.df.atac$bHLH == TRUE,]
bhlh$motif = 'bHLH_family'


fig.2j = rbind(ap1, gr, sp, twist, bhlh)

pdf(file = paste0('2j.traces.pdf'),width=14,height=4)
trellis.par.set(box.umbrella = list(lty = 1, col="black", lwd=1),
                box.rectangle = list( lwd=1.0, col="black", alpha = 1.0),
                plot.symbol = list(col="black", lwd=1.0, pch ='.'))
print(
    xyplot(value ~  sample.conditions | motif, group = genes, data = fig.2j, type = c('l'),#type = c('l','p'),
           scales=list(x=list(cex=1.0,relation = "free", rot = 45), y = list(cex=1.0, relation="free")),
           aspect=1.0,
           between=list(y=0.5, x=0.5),
           ylab = list(label = 'Normalized ATAC signal', cex =1.0),
           index.cond=list(c(1,3,4,5,3)),
           xlab = list(label = 'Time (minutes)', cex =1.0),
           par.settings = list(superpose.symbol = list(pch = c(16),col=c('grey20'), cex =0.5),
                               #I want to change the background strip to the corresponding motif color
                               strip.background=list(col=c("grey80", "green", "purple", "yellow","teal")),
                               superpose.line = list(col = as.character(cat.colours$colour), lwd=c(1),lty = c(1))),
)
)
dev.off()


#this needs to pare down because there is 7x redundancy 
pdf(file = paste0('2j.bw.pdf'),width=14,height=4)
trellis.par.set(box.umbrella = list(lty = 1, col="black", lwd=1),
                box.rectangle = list( lwd=1.0, col="black", alpha = 1.0),
                plot.symbol = list(col="black", lwd=1.0, pch ='.'))
print(
    bwplot(fig.2j$supercluster, fig.2j$AP1_family | fig.2j$motif, data = fig.2j, pch = '|', horizontal = FALSE, box.width = 10, do.out = FALSE)
)
dev.off()


