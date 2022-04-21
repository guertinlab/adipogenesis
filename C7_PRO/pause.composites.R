library(bigWig)
library(zoo)
library(lattice)

setwd('/Volumes/Dutta_external_SSD/C7_PRO/')

#isolate genes of interest from primary transcript annotation .bed and find TSS

list <- read.table('dex_activated_genes.txt',sep='\t')
act.genes <- sapply(strsplit(list[,1], '_'), '[', 2)

pTA <- read.table('primary_transcript_annotation.bed')
bed.activated <- pTA[pTA$V4 %in% act.genes,]
bed.activated <- bed.activated[,-5]
bed.activated[,6] <- bed.activated[,2]
bed.activated[bed.activated[,5] == '-',6] <- bed.activated[bed.activated[,5] == '-',3]
bed.activated[,7] <- bed.activated[,6] + 1
bed.activated[,8] <- 'Activated'

list <- read.table('dex_unchanged_genes.txt',sep='\t')
unchanged.genes <- sapply(strsplit(list[,1], '_'), '[', 2)

bed.unchanged <- pTA[pTA$V4 %in% unchanged.genes,]
bed.unchanged <- bed.unchanged[,-5]
bed.unchanged[,6] <- bed.unchanged[,2]
bed.unchanged[bed.unchanged[,5] == '-',6] <- bed.unchanged[bed.unchanged[,5] == '-',3]
bed.unchanged[,7] <- bed.unchanged[,6] + 1
bed.unchanged[,8] <- 'Unchanged'

bed.input <- rbind(bed.activated,bed.unchanged)
colnames(bed.input) = c('chr', 'start', 'end', 'gene', 'strand', 'TSS','TSS+1','type')

#define 1000 bp window centered around TSS
tss.df = bed.input[,c(1,6:7,4,8,5)]
tss.df[,3] = tss.df[,2] + 500
tss.df[,2] = tss.df[,2] - 500

#sum polymerase density for each position in window from merged .bigWigs 
combined.plus.bw = load.bigWig('/Volumes/Dutta_external_SSD/C7_PRO/C7_pro_plus_merged.bigWig')
combined.minus.bw = load.bigWig('/Volumes/Dutta_external_SSD/C7_PRO/C7_pro_minus_merged.bigWig')
dat.x = bed6.step.probeQuery.bigWig(combined.plus.bw, combined.minus.bw, tss.df, step = 1, as.matrix = TRUE, op = "sum", follow.strand = FALSE)
dat.x[is.na(dat.x)] <- 0
dat.x.win = t(apply(dat.x, 1, function(x){rollapply(x, width = 50, FUN = mean, by = 1, by.column = FALSE,align = "left")}))
index.max = apply(dat.x.win, 1, which.max)
the.max = apply(dat.x.win, 1, max)

#find highest value in window to use as pause summit and take a 50 bp window around it - this is the 'pause region'
pauseregion.df = tss.df
pauseregion.df[2][pauseregion.df$strand == '-',] = pauseregion.df[2][pauseregion.df$strand == '-',] + index.max[pauseregion.df$strand == '-']
pauseregion.df[2][pauseregion.df$strand == '+',] = pauseregion.df[2][pauseregion.df$strand == '+',] + index.max[pauseregion.df$strand == '+']
pauseregion.df[3] = pauseregion.df[2] + 1

#set values for composite calculation
upstream = 350
downstream = 1400
roll.avg = 10
step = 5

#define region set around pause summit
bedTSSwindow = fiveprime.bed(pauseregion.df, upstreamWindow = upstream,
                             downstreamWindow = downstream)

#find polymerase density in regions for each condition
composites.df <- data.frame()
for(type in unique(bedTSSwindow$type)) {
    df <- bedTSSwindow[bedTSSwindow$type == type,]
    for(bw in Sys.glob(file.path('*_plus_combined_normalized.bigWig'))) {
        cond <- strsplit(bw,'_plus')[[1]][1]
        print(cond)
        wiggle.plus = load.bigWig(bw)
        wiggle.minus = load.bigWig(paste0(cond,'_minus_combined_normalized.bigWig'))
        tss.matrix = bed6.step.bpQuery.bigWig(wiggle.plus, wiggle.minus, bedTSSwindow,
                                              step = step, as.matrix=TRUE, follow.strand=TRUE)

        coordin.start = (-upstream - 1)  + (step * roll.avg)/2
        coordin.end = downstream - (step * roll.avg)/2

        composite.lattice = data.frame(seq(coordin.start, coordin.end, by = step),
                                       rollmean(colMeans(tss.matrix), roll.avg),
                                       type = type,
                                       cond = cond,
                                       stringsAsFactors = FALSE)

        colnames(composite.lattice) = c('x', 'est', 'type', 'cond')
        composite.lattice$x = as.numeric(composite.lattice$x)

        composites.df = rbind(composites.df,composite.lattice)
    }
}
save(composites.df,file='composites.df.Rdata')

#plot 0 and 1 hr composites
plot.df <- composites.df[composites.df$type == 'Activated' & composites.df$cond %in% c('0hr','1hr_dex'),]

col.lines = c(rgb(0,0,1,1/2), rgb(1,0,0,1/2),rgb(0.1,0.5,0.05,1/2), rgb(0,0,0,1/2),rgb(1/2,0,1/2,1/2), rgb(0,1/2,1/2,1/2),rgb(1/2,1/2,0,1/2))
fill.poly = c(rgb(0,0,1,1/4),rgb(1,0,0,1/4), rgb(0.1,0.5,0.05,1/4),rgb(0,0,0,1/4), rgb(1/2,0,1/2,1/4))

pdf(file = 'composite_PolII_density_around_dex_activated_TSS_1hr.pdf',width=3.4, height=3.4) 
print(
    xyplot(est ~ x, group = cond, data = plot.df,
               type = 'l',
               scales=list(x=list(cex=1,relation = "free"), y =list(cex=1, relation="free",axs = 'i')),
               col = col.lines,
               auto.key = list(points=F, lines=T, cex=1.0,font=1),
               par.settings = list(strip.background=list(col="grey85"),
                                   superpose.symbol = list(pch = c(16),
                                                           col=col.lines, cex =0.5), 
                                   superpose.line = list(col = col.lines, lwd=c(3), 
                                                         lty = c(1))),
               par.strip.text=list(cex=0.9, font=1, col='black'),
               aspect=1.0,
               between=list(y=0.5, x=0.5),
               lwd=3,
           ylab = list(label = "RNA Polymerase Density", cex =1.2,font=2),
           xlab = list(label = "Distance from Pause Peak", cex =1.2,font=2),
           xlim = c(-150,950),
           ylim = c(0, 3.3)
           )
)
dev.off()

#plot 0 and 4 hr composites
plot.df <- composites.df[composites.df$type == 'Activated' & composites.df$cond %in% c('0hr','4hr_dex'),]

col.lines = c(rgb(0,0,1,1/2), rgb(1,0,0,1/2),rgb(0.1,0.5,0.05,1/2), rgb(0,0,0,1/2),rgb(1/2,0,1/2,1/2), rgb(0,1/2,1/2,1/2),rgb(1/2,1/2,0,1/2))
fill.poly = c(rgb(0,0,1,1/4),rgb(1,0,0,1/4), rgb(0.1,0.5,0.05,1/4),rgb(0,0,0,1/4), rgb(1/2,0,1/2,1/4))

pdf(file = 'composite_PolII_density_around_dex_activated_TSS_4hr.pdf',width=3.4, height=3.4) 
print(
    xyplot(est ~ x, group = cond, data = plot.df,
               type = 'l',
               scales=list(x=list(cex=1,relation = "free"), y =list(cex=1, relation="free",axs = 'i')),
               col = col.lines,
               auto.key = list(points=F, lines=T, cex=1.0,font=1),
               par.settings = list(strip.background=list(col="grey85"),
                                   superpose.symbol = list(pch = c(16),
                                                           col=col.lines, cex =0.5), 
                                   superpose.line = list(col = col.lines, lwd=c(3), 
                                                         lty = c(1))),
               par.strip.text=list(cex=0.9, font=1, col='black'),
               aspect=1.0,
               between=list(y=0.5, x=0.5),
               lwd=3,
           ylab = list(label = "RNA Polymerase Density", cex =1.2,font=2),
           xlab = list(label = "Distance from Pause Peak", cex =1.2,font=2),
           xlim = c(-150,950),
           ylim = c(0, 3.3)
           )
)
dev.off()
