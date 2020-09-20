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


plot.df.atac$status[plot.df.atac$supercluster == 'up.flat'] <- 'Activated'
plot.df.atac$status[plot.df.atac$supercluster == 'up.down'] <- 'Activated'
plot.df.atac$status[plot.df.atac$supercluster == 'gradual.up'] <- 'Activated'
plot.df.atac$status[plot.df.atac$supercluster == 'gradual.down'] <- 'Repressed'
plot.df.atac$status[plot.df.atac$supercluster == 'down.up'] <- 'Repressed'



ap1 = plot.df.atac[plot.df.atac$AP1 == TRUE & plot.df.atac$GR == FALSE & plot.df.atac$SP == FALSE & plot.df.atac$TWIST == FALSE & plot.df.atac$bHLH == FALSE,]
ap1$motif = 'AP1 family'
gr = plot.df.atac[plot.df.atac$AP1 == FALSE & plot.df.atac$GR == TRUE & plot.df.atac$SP == FALSE & plot.df.atac$TWIST == FALSE & plot.df.atac$bHLH == FALSE,]
gr$motif = 'GR family'
sp = plot.df.atac[plot.df.atac$AP1 == FALSE & plot.df.atac$GR == FALSE & plot.df.atac$SP == TRUE & plot.df.atac$TWIST == FALSE & plot.df.atac$bHLH == FALSE,]
sp$motif = 'SP family'
twist = plot.df.atac[plot.df.atac$AP1 == FALSE & plot.df.atac$GR == FALSE & plot.df.atac$SP == FALSE & plot.df.atac$TWIST == TRUE & plot.df.atac$bHLH == FALSE,]
twist$motif = 'TWIST family'
bhlh = plot.df.atac[plot.df.atac$AP1 == FALSE & plot.df.atac$GR == FALSE & plot.df.atac$SP == FALSE & plot.df.atac$TWIST == FALSE & plot.df.atac$bHLH == TRUE,]
bhlh$motif = 'bHLH family'


fig.2j = rbind(ap1, gr, sp, twist, bhlh)

#added these value in AI
nrow(ap1[ap1$status == 'Activated',])/7
nrow(ap1[ap1$status == 'Repressed',])/7

nrow(gr[gr$status == 'Activated',])/7
nrow(gr[gr$status == 'Repressed',])/7

nrow(sp[sp$status == 'Activated',])/7
nrow(sp[sp$status == 'Repressed',])/7

nrow(twist[twist$status == 'Activated',])/7
nrow(twist[twist$status == 'Repressed',])/7

nrow(bhlh[bhlh$status == 'Activated',])/7
nrow(bhlh[bhlh$status == 'Repressed',])/7


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
           index.cond=list(c(1,3,4,5,2)),
           xlab = list(label = 'Time (minutes)', cex =1.0),
           par.settings = list(superpose.symbol = list(pch = c(16),col=c('grey20'), cex =0.5),
                               #I want to change the background strip to the corresponding motif color
                               strip.background=list(col=c("grey80", "green", "purple", "yellow","teal")),
                               superpose.line = list(col = as.character(cat.colours$colour), lwd=c(1),lty = c(1))),
)
)
dev.off()


ap1.gr = plot.df.atac[(plot.df.atac$AP1 == TRUE | plot.df.atac$GR == TRUE) & plot.df.atac$SP == FALSE & plot.df.atac$TWIST == FALSE & plot.df.atac$bHLH == FALSE,]
ap1.gr$motif = 'AP1 or GR only'

sp.twist.bhlh = plot.df.atac[plot.df.atac$AP1 == FALSE & plot.df.atac$GR == FALSE & (plot.df.atac$SP == TRUE | plot.df.atac$TWIST == TRUE | plot.df.atac$bHLH == TRUE),]
sp.twist.bhlh$motif = 'SP or TWIST or bHLH only'

#twist.bhlh = plot.df.atac[plot.df.atac$AP1 == FALSE & plot.df.atac$GR == FALSE & plot.df.atac$SP == FALSE & (plot.df.atac$TWIST == TRUE | plot.df.atac$bHLH == TRUE),]
#twist.bhlh$motif = 'TWIST and/or bHLH '

nrow(ap1.gr[ap1.gr$status == 'Activated',])/7
nrow(ap1.gr[ap1.gr$status == 'Repressed',])/7

nrow(sp.twist.bhlh[sp.twist.bhlh$status == 'Activated',])/7
nrow(sp.twist.bhlh[sp.twist.bhlh$status == 'Repressed',])/7

#nrow(twist.bhlh[twist.bhlh$status == 'Activated',])/7
#nrow(twist.bhlh[twist.bhlh$status == 'Repressed',])/7


fig.2j.2 = rbind(ap1.gr, sp.twist.bhlh)


pdf(file = paste0('2j.traces.act.rep.pdf'),width=14,height=4.0)
trellis.par.set(box.umbrella = list(lty = 1, col="black", lwd=1),
                box.rectangle = list( lwd=1.0, col="black", alpha = 1.0),
                plot.symbol = list(col="black", lwd=1.0, pch ='.'))
print(
    xyplot(value ~  sample.conditions | motif, group = genes, data = fig.2j.2, type = c('l'),#type = c('l','p'),
           scales=list(x=list(cex=1.0,relation = "free", rot = 45), y = list(cex=1.0, relation="free")),
           aspect=1.0,
           between=list(y=0.5, x=0.5),
           ylab = list(label = 'Normalized ATAC signal', cex =1.0),
           xlab = list(label = 'Time (minutes)', cex =1.0),
           par.settings = list(superpose.symbol = list(pch = c(16),col=c('grey20'), cex =0.5),
                               #I want to change the background strip to the corresponding motif color
                               strip.background=list(col=c("grey80", "green", "purple", "yellow","teal")),
                               superpose.line = list(col = as.character(cat.colours$colour), lwd=c(1),lty = c(1))),
)
)
dev.off()

                                        #are GREs and AP1 motifs in the activated?
act.triple = unique(sp.twist.bhlh[sp.twist.bhlh$status == 'Activated',]$genes)

chrom = vector(mode="character", length = length(act.triple))
start = vector(mode="integer", length = length(act.triple))
end = vector(mode="integer", length = length(act.triple))

count = 0
for (j in act.triple) {
    count = count +1
    chrom[count] = strsplit(j, ":")[[1]][1]
    start[count] = as.numeric(strsplit(strsplit(j, ":")[[1]][2], "-")[[1]][1])
    end[count] = as.numeric(strsplit(strsplit(j, ":")[[1]][2], "-")[[1]][2])
}

res = data.frame(cbind(chrom, start, end))
colnames(res) = c('chr', 'start', 'end')
res$start = as.numeric(as.character(res$start))
res$end = as.numeric(as.character(res$end))

write.table(res, file='activated_triple.bed', sep='\t', quote=F, row.names=F, col.names=F)

#fastaFromBed -fi /Volumes/GUERTIN_2/adipogenesis/atac/mm10.fa -bed activated_triple.bed -fo activated_triple.fasta

#meme -oc activated_triple.meme_output -nmotifs 5 -objfun classic -evt 0.01 -searchsize 0 -minw 6 -maxw 15 -revcomp -dna -markov_order 3 -maxsize 100000000 activated_triple.fasta



#makes the box whiskey plot of motif scores.
fig.2j.uniq = fig.2j[!duplicated(paste0(fig.2j$genes, fig.2j$motif)),]

x = rbind(tail(head(fig.2j.uniq, 1000)), tail(head(fig.2j.uniq, 100)), tail(head(fig.2j.uniq, 7000)), tail(head(fig.2j.uniq, 6500)))
                                        #

fig.2j.uniq$new[fig.2j.uniq$AP1 == TRUE] <- fig.2j.uniq$AP1_family[fig.2j.uniq$AP1 == TRUE]
fig.2j.uniq$new[fig.2j.uniq$GR == TRUE] <- fig.2j.uniq$GR_family[fig.2j.uniq$GR == TRUE]
fig.2j.uniq$new[fig.2j.uniq$SP == TRUE] <- fig.2j.uniq$SP_family[fig.2j.uniq$SP == TRUE]
fig.2j.uniq$new[fig.2j.uniq$TWIST == TRUE] <- fig.2j.uniq$TWIST_family[fig.2j.uniq$TWIST == TRUE]
fig.2j.uniq$new[fig.2j.uniq$bHLH == TRUE] <- fig.2j.uniq$bHLH_family[fig.2j.uniq$bHLH == TRUE]

fig.2j.uniq$new <- as.numeric(as.character(fig.2j.uniq$new))



#this needs to pare down because there is 7x redundancy 
pdf(file = paste0('2j.bw.pdf'),width=14,height=4)
trellis.par.set(box.umbrella = list(lty = 1, col="black", lwd=1),
                box.rectangle = list( lwd=1.0, col="black", alpha = 1.0),
                plot.symbol = list(col="black", lwd=1.0, pch ='.'))
print(
    
bwplot(log(as.numeric(as.character(new))) ~ status | motif, data = fig.2j.uniq, horizontal = FALSE, pch = '|', do.out = FALSE,
       scales=list(x=list(cex=1.0, relation = "free", rot = 45), y = list(cex=1.0, relation="free")),
       aspect=2.0,
       between=list(y=0.5, x=0.5),
       index.cond=list(c(1,3,4,5,2)),
       ylab = expression('log'[10]* paste(Sigma,'(motif scores)')),
       xlab = expression('ATAC Peak category'),
#manually settign to avoid outlier, since do.out = FALSE       
       ylim = list(c(2.5, 3.9), c(2.4, 4.1), c(2.1, 4.1), c(2.3, 5.0), c(2.3, 3.1)),
       par.settings = list(superpose.symbol = list(pch = c(16),col=c('grey20'), cex =0.5),
                               #I want to change the background strip to the corresponding motif color
                               strip.background=list(col=c("grey80")),
                               superpose.line = list(col = as.character(cat.colours$colour), lwd=c(1),lty = c(1))))

)
dev.off()


