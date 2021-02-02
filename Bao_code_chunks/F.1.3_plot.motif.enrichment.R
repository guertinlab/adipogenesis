library(lattice)
library(bigWig)

dir = '/scratch/bhn9by/ATAC/fimo_composites/'
setwd(dir)

bed.window <- function(bed, half.window) {
    bed[,2] = (bed[,2] + bed[,3])/2 - half.window
    bed[,3] = bed[,2] + 2 * half.window
                                        #This shouldn't be necessary, but I used it as a workaround while troubleshooting
    bed = subset(bed, bed[,2] > 0)
    return(bed)
}

plot.fimo.lattice <- function(dat, fact = 'Motif', summit = 'Hypersensitivity Summit', class= '',
                              num.m = -200, num.p =90, y.low =0, y.high = 0.2,
                              col.lines = c(rgb(0,0,1,1/2), rgb(1,0,0,1/2),
                                            rgb(0.1,0.5,0.05,1/2), rgb(0,0,0,1/2),
                                            rgb(1/2,0,1/2,1/2), rgb(0,1/2,1/2,1/2), rgb(1/2,1/2,0,1/2)),
                              fill.poly = c(rgb(0,0,1,1/4), rgb(1,0,0,1/4), rgb(0.1,0.5,0.05,1/4),
                                            rgb(0,0,0,1/4), rgb(1/2,0,1/2,1/4))) {
    
    pdf('motif_enrichment_around_summits.pdf')#, width=6.83, height=3.5)
    print(xyplot(density ~ range|tf, groups = category, data = dat, type = 'l',
                 scales=list(x=list(cex=0.8,relation = "free"), y =list(cex=0.8,axs = 'i',relation = "free")),
                 xlim=c(num.m,num.p),                                       
                 col = col.lines,
                 auto.key = list(points=F, lines=T, cex=0.8, columns = 2),
                 par.settings = list(superpose.symbol = list(pch = c(16), col=col.lines, cex =0.7),
                                     superpose.line = list(col = col.lines, lwd=c(2,2,2,2,2,2),
                                                           lty = c(1,1,1,1,1,1,1,1,1))),
                 cex.axis=1.0,
                 par.strip.text=list(cex=0.9, font=1, col='black',font=2),
                 aspect=1.0,
                 between=list(y=0.5, x=0.5),
                 index.cond = list(c(4:6,1:3)),
                 lwd=2,
                 ylab = list(label = "Weighted Motif Density", cex =1,font=2),
                 xlab = list(label = 'Distance from ATAC-seq Peak Summit', cex =1,font=2),
                 strip = function(..., which.panel, bg) {
                     bg.col = 'grey'#c("blue","grey65","red")
                     strip.default(..., which.panel = which.panel, bg = rep(bg.col, length = which.panel)[which.panel])
                 }
                 ))
    dev.off()
}

load('/scratch/bhn9by/ATAC/plot.df.atac.Rdata')
load('/scratch/bhn9by/ATAC/fimo.scores.atac.Rdata')

plot.df = data.frame()
for(i in 1:ncol(fimo.scores.atac)) {
    temp = plot.df.atac[plot.df.atac$genes %in% rownames(fimo.scores.atac[!is.na(fimo.scores.atac[,i]),]),]
    temp$family = colnames(fimo.scores.atac)[i]
    plot.df = rbind(plot.df,temp)
}

plot.df$status = 'Activated'
plot.df[plot.df$supercluster == 'gradual.down' | plot.df$supercluster == 'down.up',]$status = 'Repressed'

all.fimo = data.frame(matrix(ncol = 4, nrow = 0))
colnames(all.fimo) = c('density', 'tf', 'category', 'range')

half.win = 600
file.suffix = '_mm10_instances.bigWig'
dir = '/scratch/bhn9by/ATAC/fimo_composites/main_figure_beds/'

decreased = plot.df[plot.df$status == 'Repressed',7:9]
decreased[,2] = as.numeric(decreased[,2])
decreased[,3] = as.numeric(decreased[,3])
decreased = bed.window(decreased,half.win)

increased = plot.df[plot.df$status == 'Activated',7:9]
increased[,2] = as.numeric(increased[,2])
increased[,3] = as.numeric(increased[,3])
increased = bed.window(increased,half.win)

not.different = read.table('/scratch/bhn9by/ATAC/nondynamic_peaks.bed')
not.different = not.different[not.different$V1 != 'chrM',]
not.different = bed.window(not.different,half.win)

all.fimo = data.frame()

for(i in 1:ncol(fimo.scores.atac)) {
    factor = colnames(fimo.scores.atac)[i]

    mod.bigWig = paste0(dir,factor,file.suffix)
    factor.name = factor
    print(factor.name)

    loaded.bw = load.bigWig(mod.bigWig)
    
    dec.inten = bed.step.probeQuery.bigWig(loaded.bw, decreased,
                                           gap.value = 0, step = 10, as.matrix = TRUE)
    dec.query.df = data.frame(cbind(colMeans(dec.inten), factor.name,
                                    'Closed', seq(-half.win, (half.win-10), 10)), stringsAsFactors=F)
    colnames(dec.query.df) = c('density', 'tf', 'category', 'range')
    
    inc.inten = bed.step.probeQuery.bigWig(loaded.bw, increased,
                                           gap.value = 0, step = 10, as.matrix = TRUE)
    inc.query.df = data.frame(cbind(colMeans(inc.inten), factor.name,
                                    'Opened', seq(-half.win,(half.win-10), 10)), stringsAsFactors=F)
    colnames(inc.query.df) = c('density', 'tf', 'category', 'range')
    
    ctrl.inten = bed.step.probeQuery.bigWig(loaded.bw, not.different,
                                            gap.value = 0, step = 10, as.matrix = TRUE)
    ctrl.query.df = data.frame(cbind(colMeans(ctrl.inten), factor.name,
                                     'Nondynamic', seq(-half.win, (half.win-10), 10)), stringsAsFactors=F)
    colnames(ctrl.query.df) = c('density', 'tf', 'category', 'range')
    
    tf.all = rbind(dec.query.df, inc.query.df, ctrl.query.df)

    all.fimo = rbind(all.fimo,tf.all)
}

all.fimo[,1] = as.numeric(all.fimo[,1])
all.fimo[,4] = as.numeric(all.fimo[,4])

plot.fimo.lattice(all.fimo, num.m = -500, num.p = 500,
                  col.lines = c('blue','grey65','red'))
