library(lattice)
library(MASS)

setwd('/scratch/bhn9by/ATAC/SP_KLF_split')

load('/scratch/bhn9by/ATAC/plot.df.atac.Rdata')

sp.mast = read.table('output_sp1.txt', stringsAsFactors=FALSE, sep = '\t', header = TRUE)
klf.mast = read.table('output_klf.txt', stringsAsFactors=FALSE, sep = '\t', header = TRUE)

sp.mast = sp.mast[sp.mast[,6] != -1,]
sp.mast = aggregate(as.numeric(score)~sequence_name, data=sp.mast, FUN=max)
rownames(sp.mast) = sp.mast[,1]

klf.mast = klf.mast[klf.mast[,6] != -1,]
klf.mast = aggregate(as.numeric(score)~sequence_name, data=klf.mast, FUN=max)
rownames(klf.mast) = klf.mast[,1]

sp.data = merge(sp.mast, klf.mast, by = 1)
colnames(sp.data) = c('seq', 'sp', 'klf')

xyplot(sp.data$sp ~ sp.data$klf, pch = 16, cex = 0.5, xlim = c(-5,25), ylim = c(-5,25), aspect = 1, panel = function(x, y, ...) {
         panel.xyplot(x, y, ...)
         panel.abline(0, 1)
})

out.klf <- strsplit(as.character(klf.mast[,1]), ':') 
all.klf.motifs <- data.frame(do.call(rbind, out.klf, quote=FALSE), klf.mast[,c(2)])
colnames(all.klf.motifs) = c(colnames(all.klf.motifs)[1:9], 'score')
all.klf.temp = all.klf.motifs

all.klf.temp$peaks = paste0(all.klf.temp$X1,':',all.klf.temp$X2,'-',all.klf.temp$X3)

for (peak in unique(all.klf.temp$peaks)) {
    temp = all.klf.temp[all.klf.temp$peaks == peak,]
    if(nrow(temp) != 1) {
        all.klf.temp = all.klf.temp[-which(all.klf.temp$peaks == peak),]
    }
}

out.sp <- strsplit(as.character(sp.mast[,1]), ':') 
all.sp.motifs <- data.frame(do.call(rbind, out.sp, quote=FALSE), sp.mast[,c(2)])
colnames(all.sp.motifs) = c(colnames(all.sp.motifs)[1:9], 'score')
all.sp.temp = all.sp.motifs

all.sp.temp$peaks = paste0(all.sp.temp$X1,':',all.sp.temp$X2,'-',all.sp.temp$X3)

for (peak in unique(all.sp.temp$peaks)) {
    temp = all.sp.temp[all.sp.temp$peaks == peak,]
    if(nrow(temp) != 1) {
        all.sp.temp = all.sp.temp[-which(all.sp.temp$peaks == peak),]
    }
}

sp.motifs.lda = aggregate(as.numeric(score)~X1+X2+X3, data=all.sp.temp, FUN=max)
klf.motifs.lda = aggregate(as.numeric(score)~X1+X2+X3, data=all.klf.temp, FUN=max)
colnames(sp.motifs.lda) = c('chr', 'start', 'end', 'sp')
colnames(klf.motifs.lda) = c('chr', 'start', 'end', 'klf')

rownames(klf.motifs.lda) = paste0(klf.motifs.lda[,1], ':', klf.motifs.lda[,2], '-', klf.motifs.lda[,3])
rownames(sp.motifs.lda) = paste0(sp.motifs.lda[,1], ':', sp.motifs.lda[,2], '-', sp.motifs.lda[,3])

chr = sapply(strsplit(unique(plot.df.atac$genes), ':'), '[', 1)
x = sapply(strsplit(unique(plot.df.atac$genes), ':'), '[', 2)
start = as.numeric(sapply(strsplit(x, '-'), '[', 1))
end = as.numeric(sapply(strsplit(x, '-'), '[', 2))

res.lda = data.frame(chr, start, end)
rownames(res.lda) = unique(plot.df.atac$genes)

res.lda = merge(res.lda, klf.motifs.lda, by="row.names", all = TRUE)
rnames = res.lda[,1]

res.lda = data.frame(res.lda[,c(2:4, 8)])
rownames(res.lda) = rnames

res.lda = merge(res.lda, sp.motifs.lda, by="row.names", all = TRUE)
rnames = res.lda[,1]

res.lda = data.frame(res.lda[,c(2:5,9)])
rownames(res.lda) = rnames

save(res.lda, file = 'res.lda.Rdata')

res.lda = res.lda[!is.na(res.lda$sp) & !is.na(res.lda$klf),]

func <- function(peak) {
    status = 'Activated'
    sc = unique(plot.df.atac[plot.df.atac$genes == peak,]$supercluster)
    if(sc == 'gradual.down' | sc == 'down.up') {
        status = 'Repressed'
    }
    return(status)
}

res.lda$status = sapply(rownames(res.lda),func)

xyplot(sp ~ klf, groups = status, data= res.lda, pch = 16, col = c('red', 'blue'),
       panel = function(x, y, ...) {
         panel.xyplot(x, y, ...)
       })

fit <- lda(status ~ sp + klf, data= res.lda)

plot(fit)

gmean <- fit$prior%*%fit$means
const <- drop(gmean%*%fit$scaling)

slope <- -fit$scaling[1]/fit$scaling[2]
intercept <- const/fit$scaling[2]

slope
#1.551901
#1.910971
intercept
#-5.454026
#-10.8391

res.lda$status = as.factor(res.lda$status)

relevel(res.lda$status, "Repressed")

pdf(file = paste0('LDA_analysis.pdf'),width=3,height=3)

print(
    
    xyplot(klf ~ sp, groups = status, data= res.lda, pch = 16, cex = 0.6,
           col = c('#FF000066', '#0000FF40'),
#       ylim = c(-5, 25),
#       xlim = c(-5, 25),
       xlab = "SP score",
       ylab = "KLF score",
       scales=list(alternating=c(1,0)),
       panel = function(x, y, ...) {
         panel.xyplot(x, y, ...)
         panel.abline(intercept , slope, lty = 2, col = 'grey45', lwd = 2.3)
       })

)
dev.off()

klf.motifs = all.klf.motifs[(all.klf.motifs$score > (all.sp.motifs$score*slope + intercept)),]
sp.motifs = all.sp.motifs[(all.klf.motifs$score <= (all.sp.motifs$score*slope + intercept)),]

chr = sapply(strsplit(unique(plot.df.atac$genes), ':'), '[', 1)
x = sapply(strsplit(unique(plot.df.atac$genes), ':'), '[', 2)
start = as.numeric(sapply(strsplit(x, '-'), '[', 1))
end = as.numeric(sapply(strsplit(x, '-'), '[', 2))

res = data.frame(chr, start, end)
rownames(res) = unique(plot.df.atac$genes)

sp.motifs.1 = aggregate(as.numeric(score)~X1+X2+X3, data=sp.motifs, FUN=sum)
klf.motifs.1 = aggregate(as.numeric(score)~X1+X2+X3, data=klf.motifs, FUN=sum)

rownames(sp.motifs.1) = paste0(sp.motifs.1[,1], ':', sp.motifs.1[,2], '-', sp.motifs.1[,3])
rownames(klf.motifs.1) = paste0(klf.motifs.1[,1], ':', klf.motifs.1[,2], '-', klf.motifs.1[,3])

colnames(sp.motifs.1) = c('chr', 'start', 'end', 'sp')
colnames(klf.motifs.1) = c('chr', 'start', 'end', 'klf')

res = merge(res, sp.motifs.1, by="row.names", all = TRUE)
rnames = res[,1]

res = data.frame(res[,c(2:4, 8)])

colnames(res) = c('chr', 'start','end', 'SP_alone')
rownames(res) = rnames

res = merge(res, klf.motifs.1, by="row.names", all = TRUE)
rnames = res[,1]

res = data.frame(res[,c(2:5, 9)])

colnames(res) = c('chr', 'start','end', 'SP_alone', 'KLF_alone')
rownames(res) = rnames

sp.klf.scores.atac = res
save(sp.klf.scores.atac, file = 'sp.klf.scores.atac.Rdata')

#now all instances
sp.mast.nondyn = read.table('output_sp1_2M.txt', stringsAsFactors=FALSE, sep = '\t', header = TRUE)
klf.mast.nondyn = read.table('output_klf_2M.txt', stringsAsFactors=FALSE, sep = '\t', header = TRUE)

sp.mast.nondyn = aggregate(as.numeric(score)~sequence_name, data=sp.mast.nondyn, FUN=max)
rownames(sp.mast.nondyn) = sp.mast.nondyn[,1]

klf.mast.nondyn = aggregate(as.numeric(score)~sequence_name, data=klf.mast.nondyn, FUN=max)
rownames(klf.mast.nondyn) = klf.mast.nondyn[,1]

sp.data = merge(sp.mast.nondyn, klf.mast.nondyn, by = 1)
colnames(sp.data) = c('seq', 'sp', 'klf')

slope = 1.551901
intercept = -5.454026

klf.motifs = sp.data[(sp.data$klf > (sp.data$sp*slope + intercept)),]
sp.motifs = sp.data[(sp.data$klf <= (sp.data$sp*slope + intercept)),]

chr = sapply(strsplit(sp.motifs$seq, ':'), '[', 1)
start = sapply(strsplit(sp.motifs$seq, ':'), '[', 2)
end = sapply(strsplit(sp.motifs$seq, ':'), '[', 3)

sp.bed = data.frame(chr,start,end,sp.motifs$sp)

chr = sapply(strsplit(klf.motifs$seq, ':'), '[', 1)
start = sapply(strsplit(klf.motifs$seq, ':'), '[', 2)
end = sapply(strsplit(klf.motifs$seq, ':'), '[', 3)

klf.bed = data.frame(chr,start,end,klf.motifs$klf)

write.table(sp.bed,file='/scratch/bhn9by/ATAC/fimo_composites/main_figure_beds/SP_unsorted.bed',quote=F,sep='\t',row.names=F,col.names=F)

write.table(klf.bed,file='/scratch/bhn9by/ATAC/fimo_composites/main_figure_beds/KLF_unsorted.bed',quote=F,sep='\t',row.names=F,col.names=F)
