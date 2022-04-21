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
pauseregion.df[3] = pauseregion.df[2] + 50  

#define 'body region' as end of pause region to the end of the gene (strand specific)
body.df = bed.input[,c(1:4,8,5)]
body.df[body.df$strand == '+',2] = pauseregion.df[pauseregion.df$strand == '+',3]+1
body.df[body.df$strand == '-',3] = pauseregion.df[pauseregion.df$strand == '-',2]-1

#measure polymerase density within pause and body regions (sum for pause, average for body)
#use .bigwigs from each condition and save density separately
pause.df <- data.frame()
for(bw in Sys.glob(file.path('*_plus_combined_normalized.bigWig'))) {
    cond <- strsplit(bw,'_plus')[[1]][1]
    wiggle.plus = load.bigWig(bw)
    wiggle.minus = load.bigWig(paste0(cond,'_minus_combined_normalized.bigWig'))
    
    temp.df = bed.input[,c(1:4,8,5)]
    temp.df[,7] = bed6.region.bpQuery.bigWig(wiggle.plus, wiggle.minus, pauseregion.df, op = "sum")
    temp.df[,8] = bed6.region.bpQuery.bigWig(wiggle.plus, wiggle.minus, body.df, op = "avg")
    temp.df[,9] = cond
    colnames(temp.df)[7:9] = c('pause_sum','body_avg','cond')
    pause.df <- rbind(pause.df,temp.df)
}

#find pause index for each gene of interest 
pause.df$pause.index = pause.df[,'pause_sum'] / pause.df[,'body_avg']
save(pause.df,file='pause.df.Rdata')

#define new plotting df with pause index fold changes b/w measurements
plot.df <- data.frame(genes = unique(pause.df$gene),
                      type = '',
                      fc.1hr = 0,
                      fc.4hr = 0)

for(i in 1:nrow(plot.df)) {
    gene <- plot.df$gene[i]
    plot.df$type[i] <- unique(pause.df[pause.df$gene == gene,'type'])
    plot.df$fc.1hr[i] <- pause.df[pause.df$gene == gene & pause.df$cond == '1hr_dex','pause.index'] / pause.df[pause.df$gene == gene & pause.df$cond == '0hr','pause.index']
    plot.df$fc.4hr[i] <- pause.df[pause.df$gene == gene & pause.df$cond == '4hr_dex','pause.index'] / pause.df[pause.df$gene == gene & pause.df$cond == '0hr','pause.index']
}
plot.df <- plot.df[!plot.df$fc.1hr %in% c(Inf,NaN),]
plot.df <- plot.df[!(plot.df$type == 'Unchanged'),]

#plot
pdf('dex.pause.index.fc.1hr.pdf', width=2.5, height=3)

trellis.par.set(box.umbrella = list(lty = 1, col="#93939380", lwd=2),
                box.rectangle = list(col = '#93939380', lwd=1.6),
                plot.symbol = list(col='#93939380', lwd=1.6, pch ='.'))

print(bwplot(log(fc.1hr) ~ type, data = plot.df,
             between=list(y=1.0, x = 1.0),
             scales=list(x=list(draw=FALSE),relation="free",rot = 45, alternating=c(1,1,1,1),cex=1,font=1),
             xlab = '',
            #main = "Pause Index (PI) Ratio",
             ylab = expression("log"[2]~"(60min PI / 0min PI)"),
             horizontal = FALSE,
             col= 'black',
             aspect = 2,
             ylim = c(-1.52, 1),
             par.settings=list(par.xlab.text=list(cex=1.2,font=1),
                               par.ylab.text=list(cex=1.2,font=1),
                               par.main.text=list(cex=1.2, font=1),
                               plot.symbol = list(col='black', lwd=1.6, pch =19, cex = 0.25)),
             strip = function(..., which.panel, bg) {
               #bg.col = c("#ce228e" ,"grey60", "#2290cf","grey90")
               strip.default(..., which.panel = which.panel,
                             bg = rep(bg.col, length = which.panel)[which.panel])
             },
             panel = function(..., box.ratio, col) {
                 panel.abline(h = 0, col = 'grey45', lty = 2)
                 panel.violin(..., col = 'light blue',
                              varwidth = FALSE, box.ratio = box.ratio, outer = FALSE)
                 panel.stripplot(..., col='#54545380', do.out=FALSE, jitter.data=TRUE, amount = 0.2, pch = 16)
                 panel.bwplot(..., pch = '|', do.out = FALSE)    
       }))
dev.off()

pdf('dex.pause.index.fc.4hr.pdf', width=2.5, height=3)

trellis.par.set(box.umbrella = list(lty = 1, col="#93939380", lwd=2),
                box.rectangle = list(col = '#93939380', lwd=1.6),
                plot.symbol = list(col='#93939380', lwd=1.6, pch ='.'))

print(bwplot(log(fc.4hr) ~ type, data = plot.df,
             between=list(y=1.0, x = 1.0),
             scales=list(x=list(draw=FALSE),relation="free",rot = 45, alternating=c(1,1,1,1),cex=1,font=1),
             xlab = '',
            #main = "Pause Index (PI) Ratio",
             ylab = expression("log"[2]~"(240min PI / 0min PI)"),
             horizontal = FALSE,
             col= 'black',
             aspect = 2,
             ylim = c(-1.7, 1.1),
             par.settings=list(par.xlab.text=list(cex=1.2,font=1),
                               par.ylab.text=list(cex=1.2,font=1),
                               par.main.text=list(cex=1.2, font=1),
                               plot.symbol = list(col='black', lwd=1.6, pch =19, cex = 0.25)),
             strip = function(..., which.panel, bg) {
               #bg.col = c("#ce228e" ,"grey60", "#2290cf","grey90")
               strip.default(..., which.panel = which.panel,
                             bg = rep(bg.col, length = which.panel)[which.panel])
             },
             panel = function(..., box.ratio, col) {
                 panel.abline(h = 0, col = 'grey45', lty = 2)
                 panel.violin(..., col = 'light blue',
                              varwidth = FALSE, box.ratio = box.ratio, outer = FALSE)
                 panel.stripplot(..., col='#54545380', do.out=FALSE, jitter.data=TRUE, amount = 0.2, pch = 16)
                 panel.bwplot(..., pch = '|', do.out = FALSE)
           
       }))
dev.off()
