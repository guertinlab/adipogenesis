
                                        #arun_PE1_H9.dREG.peak.full.bed
source('https://raw.githubusercontent.com/guertinlab/seqOutBias/master/docs/R/seqOutBias_functions.R')
source('https://raw.githubusercontent.com/mjg54/znf143_pro_seq_analysis/master/docs/ZNF143_functions.R')


                                        #bigwigs: 3T3_plus_dREG_merged.bigWig and 3T3_minus_dREG_merged.bigWig

create.composites.heatmaps <- function(path.dir, composite.input, region=20, step=1, grp = 'PRO', first=FALSE, file.prefix = 'sh') { 
    vec.names = c('chr','start','end')
    hmap.data = list()
    composite.df=data.frame(matrix(ncol = 6, nrow = 0))
    for (mod.bigWig.plus in Sys.glob(file.path(path.dir, paste(file.prefix, "*plus_dREG_merged.bigWig", sep='')))) {
        mod.bigWig.minus=(paste(strsplit(mod.bigWig.plus, 'plus')[[1]][1], 'minus_dREG_merged.bigWig', sep=''))
        factor.name = strsplit(strsplit(mod.bigWig.plus, "/")[[1]][length(strsplit(mod.bigWig.plus, "/")[[1]])], '_plus')[[1]][1]
        print(factor.name)
        print(mod.bigWig.minus)
        print(mod.bigWig.plus)
        vec.names = c(vec.names, factor.name)
        wiggle.plus = load.bigWig(mod.bigWig.plus)
        wiggle.minus = load.bigWig(mod.bigWig.minus)
        bpqueryPM = window.step.matrix(composite.input, wiggle.plus, wiggle.minus, region, step)
        subsample = subsampled.quantiles.metaprofile.adjust((bpqueryPM[[1]]))
        mult.row = ncol(bpqueryPM[[1]])
        hmap.data[[paste(factor.name,'_plus', sep='')]] = list(colMeans(bpqueryPM[[1]]), subsample$top, subsample$bottom, colMeans(bpqueryPM[[1]]), paste(factor.name,'_plus', sep=''), bpqueryPM[[1]])
        df.up <- data.frame(matrix(ncol = 6, nrow = mult.row))
        df.up[, 1] <- colMeans(bpqueryPM[[1]])
        df.up[, 2] <- seq((-1 * region) + 0.5 * step, region - 0.5 * step, by = step)
        df.up[, 3] <- matrix(data = paste(factor.name,'', sep=''), nrow=mult.row, ncol=1)
        df.up[, 4] <- subsample$top
        df.up[, 5] <- subsample$bottom
        df.up[, 6] <- matrix(data = 'Plus', nrow=mult.row, ncol=1)
        composite.df = rbind(composite.df, df.up)
        bpqueryMP = window.step.matrix(composite.input, wiggle.minus, wiggle.plus, region, step)
        subsample = subsampled.quantiles.metaprofile.adjust((bpqueryMP[[1]]))
        mult.row = ncol(bpqueryMP[[1]])
        hmap.data[[paste(factor.name,'_minus', sep='')]] = list(colMeans(bpqueryMP[[1]]), subsample$top, subsample$bottom, colMeans(bpqueryMP[[1]]), paste(factor.name,'_minus', sep=''), bpqueryMP[[1]])
        df.up <- data.frame(matrix(ncol = 6, nrow = mult.row))
        df.up[, 1] <- colMeans(bpqueryMP[[1]])
        df.up[, 2] <- seq((-1 * region) + 0.5 * step, region - 0.5 * step, by = step)
        df.up[, 3] <- matrix(data = paste(factor.name,'', sep=''), nrow=mult.row, ncol=1)
        df.up[, 4] <- subsample$top
        df.up[, 5] <- subsample$bottom
        df.up[, 6] <- matrix(data = 'Minus', nrow=mult.row, ncol=1)
        composite.df = rbind(composite.df, df.up)
        unload.bigWig(wiggle.plus)
        unload.bigWig(wiggle.minus)

    }
    colnames(composite.df) <- c('est', 'x', 'cond', 'upper', 'lower', 'grp')
    for (cond in (1:length(hmap.data))) {
    rownames(hmap.data[[cond]][[6]]) = paste(composite.input[,1], ':',
                composite.input[,2], '-', composite.input[,3], sep='')
    colnames(hmap.data[[cond]][[6]]) = seq((-1 * region) + 0.5 * step, region - 0.5 * step, by = step)
  }
    return(list(composite.df, hmap.data))
}

draw.heatmap.both.pro.overlap <- function(var.name, ord = var.name, filename = "hmap.data.pdf",
                              width.vp = 3, width.space = 0.8, height.vp = NULL,
                              height.legend = 1, upstream = -1200, downstream = 1200,
                              sub.set = 1000, order.layout = c(1:(length(var.name)/2)), avg.cols = 5, avg.rows = 40, order.buffer = 6, 
                              order.indx = 1, arrangement = c("left", "center", "custom", "asis", "pro"), max.range = -Inf, min.range = Inf,
                              colors.ramp = c("#ffffff80", "#FCDBDB80", "#FCC0C080",
                                              "#F98F8F80", "#F95A5A80", "#ff000080", "#cd000080", "#8b000080","#3F020280"),
                              colors.ramp.2 = c("#ffffff80", "#DBDBFC80", "#C0C0FC80",
                                              "#8F8FF980", "#5A5AF980",  "#0000ff80", "#0000cd80", "#00008b80","#02023F80"),
                              convert.zeros = TRUE) {
    arrangement <- match.arg(arrangement)
    plot.num = length(var.name)/2
    print('plot.num')
    print(plot.num)
    if (is.null(height.vp)) {
        height.vp = nrow(var.name[[1]][[6]])*0.0002
    }
    pdf(paste('', filename, sep=''),w=((width.vp+width.space)*plot.num)+(width.space),
        h=(height.vp+height.legend*2))
    grid.newpage()
    step = abs(as.numeric(as.character(colnames(var.name[[1]][[6]])[1])) -
        as.numeric(as.character(colnames(var.name[[1]][[6]])[2])))
    probes = ((downstream-upstream)/step)/2
        playout = grid.layout(3,(2*plot.num)+1, widths=unit(c(width.space,
                         rep(c(width.vp, width.space), plot.num)),
                         c("inches", rep(c("inches", "inches"), plot.num))),
                         heights=unit(c(height.legend, height.vp,height.legend),
                         c("inches","inches", "inches")), respect =
        matrix(data =1, nrow=3, ncol=(plot.num*2)+1))
    pushViewport(viewport(layout=playout))
    #grid.rect()
    plot.list = list()
    x.ord = apply(var.name[[order.indx]][[6]],1,which.max)
    #average over columns:
    if (arrangement == 'pro') {probes = (probes / avg.cols)}
    
    for (cond in (1:length(var.name))) {
        if (arrangement == 'left') {
            var.name[[cond]][[6]] = var.name[[cond]][[6]][order(x.ord),]
        }
        if (arrangement == 'center') {
            var.name[[cond]][[6]] = var.name[[cond]][[6]][
                                order(-ord[[order.indx]][[6]][,ncol(
                                    var.name[[order.indx]][[6]])/2]), ]
        }
                                        #this is offsetting by 150, if the window is less than 150 stuff gets messed up
        
        if (arrangement == 'custom') {
            var.name[[cond]][[6]] = var.name[[cond]][[6]][
                                order(-ord[[order.indx]][[6]][,ceiling((ncol(
                                    var.name[[order.indx]][[6]])/2) - 100/((upstream - downstream)/ncol(
                                    var.name[[order.indx]][[6]])))]), ]
        }
        if (arrangement == 'asis') {
            print('order of rows is predefined')
        }
        if (arrangement == 'pro') {
#             plus.order = order(-ord[[order.indx]][[6]][,ceiling((ncol(
#                                    var.name[[order.indx]][[6]])/2) - 100/((upstream - downstream)/ncol(
                                        #                                                                                                       var.name[[order.indx]][[6]])))])[rep(c(TRUE, FALSE))]
                                    minus.order = order(-rowSums(ord[[order.indx +1 ]][[6]][,(ceiling((ncol(
                                 var.name[[order.indx + 1]][[6]])/2) + 150/((upstream - downstream)/ncol(
                                   var.name[[order.indx +1]][[6]]))) - order.buffer): (ceiling((ncol(
                                                    var.name[[order.indx +1]][[6]])/2) + 150/((upstream - downstream)/ncol(
                                                                                                                       var.name[[order.indx +1 ]][[6]]))) + order.buffer)]))[rep(c(TRUE, FALSE))]
              plus.order = order(-rowSums(ord[[order.indx]][[6]][,(ceiling((ncol(
                                 var.name[[order.indx]][[6]])/2) - 150/((upstream - downstream)/ncol(
                                   var.name[[order.indx]][[6]]))) - order.buffer): (ceiling((ncol(
                                                    var.name[[order.indx]][[6]])/2) - 150/((upstream - downstream)/ncol(
                                                                                                                       var.name[[order.indx]][[6]]))) + order.buffer)]))
#             minus.order = order(-ord[[order.indx + 1]][[6]][,ceiling((ncol(
#                                    var.name[[order.indx + 1]][[6]])/2) + 100/((upstream - downstream)/ncol(
#                                                                                                           var.name[[order.indx +1]][[6]])))])

             plus.not.in.minus = plus.order[!(plus.order %in% minus.order)]
             final.order = as.vector(t(data.frame(minus.order, plus.not.in.minus)))
             var.name[[cond]][[6]] = var.name[[cond]][[6]][final.order,]
             var.name[[cond]][[6]] = zeroToOne.scale(var.name[[cond]][[6]])
             #averages every 40 rows to make the heatmap a reasonable size
             var.name[[cond]][[6]] = aggregate(var.name[[cond]][[6]],list(rep(1:(nrow(var.name[[cond]][[6]])%/%avg.rows+1),each=avg.rows,len=nrow(var.name[[cond]][[6]]))),mean)[-1];
             var.name[[cond]][[6]] = t(aggregate(t(var.name[[cond]][[6]]),list(rep(1:(nrow(t(var.name[[cond]][[6]]))%/%avg.cols+1),each=avg.cols,len=nrow(t(var.name[[cond]][[6]])))),mean)[-1])
                                        #             var.name[[cond]][[6]] = t(mm)
             print(dim(var.name[[cond]][[6]]))
        }
        
        if (nrow(var.name[[cond]][[6]]) > sub.set) {
            first.k = var.name[[cond]][[6]][1:sub.set,
                ((ncol(var.name[[cond]][[6]])/2)-probes):
                ((ncol(var.name[[cond]][[6]])/2)+probes)]
        } else {
        first.k = var.name[[cond]][[6]][,
            ((ncol(var.name[[cond]][[6]])/2)-probes):
            ((ncol(var.name[[cond]][[6]])/2)+probes)]
    }
        plot.list[[cond]] = first.k
    }
#   min.range = Inf
#    max.range = -Inf
    if (max.range == -Inf & min.range == Inf ) {
        for (cond in (1:length(var.name))) {
            min.cond = min(plot.list[[cond]][plot.list[[cond]] != 0], na.rm=TRUE)
            max.cond = max(plot.list[[cond]])
            if (min.cond < min.range) {
                min.range = min.cond
            }
            if (max.cond > max.range) {
                max.range = max.cond
            }
        }
    }
    print(min.range)
    print(max.range)
  #plot it!
    count = 0
    #print(length(var.name))
    for (cond in (1:(length(var.name)/2))) {
        count = count +1 
        print(cond)
        print("order.layout[count]*2")
        print(order.layout[count]*2)
        pushViewport(viewport(layout.pos.col=(order.layout[count]*2), layout.pos.row=2))
        
        first.k = plot.list[[(2*cond)-1]][nrow(plot.list[[(2*cond)-1]]):1,]
        first.k.2 = plot.list[[(2*cond)]][nrow(plot.list[[(2*cond)]]):1,]
        if (convert.zeros) {
            first.k[first.k==0] <- min.range
        }
                                        #      first.k = scale(first.k, center = FALSE, scale = TRUE)
                                        #                print(levelplot(t(log(first.k, base = 2)),
        print(levelplot(t(first.k),
                        aspect = height.vp/width.vp,
                        col.regions = colorRampPalette(colors.ramp, bias=1, alpha = TRUE)(150),
                        at = seq(0, 1, length=150), 
                        xlab="",
                        axes = FALSE,
                        par.settings = list(axis.line = list(lty = 0)),
                        ylab="",
                        main='',
                        sub="",
                        colorkey = FALSE,
                        region = TRUE,
                        scales = list(draw = FALSE)),
              panel.width = list(width.vp, "inches"),
              panel.height=list(height.vp, "inches"),
              newpage = FALSE)
                                        #        print(levelplot(t(log(first.k.2, base = 2)),
        print(levelplot(t(first.k.2),

                        aspect = height.vp/width.vp,
                        col.regions = colorRampPalette(colors.ramp.2, bias=1, alpha = TRUE)(150),
                        at = seq(0, 1, length=150), 
                        xlab="",
                        axes = FALSE,
#                        par.settings = list(axis.line = list(lty = 0)),
                       ylab="",
                        main='',
                        sub="",
                        colorkey = FALSE,
                        region = TRUE,
                        scales = list(draw = FALSE)),
              panel.width = list(width.vp, "inches"),
              panel.height=list(height.vp, "inches"),
              newpage = FALSE)


        upViewport()
        
                                        #scale
        pushViewport(viewport(layout.pos.col=(order.layout[count]*2), layout.pos.row=3))
        #grid.rect()
        #grid.text(cond)
        grid.lines(x = unit(c(0,1), "npc"),
                   y = unit(c(0.9,0.9), "npc"), gp=gpar(lwd=3))
        grid.lines(x = unit(c(0,0), "npc"),
                   y = unit(c(0.8, 0.9), "npc"), gp=gpar(lwd=3))
        grid.lines(x = unit(c(1,1), "npc"),
                   y = unit(c(0.8, 0.9), "npc"), gp=gpar(lwd=3))
        grid.lines(x = unit(c(0.5,0.5), "npc"),
                   y = unit(c(0.8, 0.9), "npc"), gp=gpar(lwd=3))
        grid.text(as.character(upstream), x = unit(0, "npc"),
                  y = unit(0.6, "npc"), gp = gpar(fontsize = 12, fontface = "bold"))
        grid.text(as.character(downstream), x = unit(1, "npc"),
                  y = unit(0.6, "npc"), gp = gpar(fontsize = 12, fontface = "bold"))
        grid.text("0", x = unit(0.5, "npc"),
                  y = unit(0.6, "npc"), gp = gpar(fontsize = 12, fontface = "bold"))
        
        upViewport()
    
        pushViewport(viewport(layout.pos.col=(order.layout[count]*2), layout.pos.row=1))

        grid.text(strsplit(var.name[[(2*cond)-1]][[5]][1], '_plus')[[1]], x = unit(0.5, "npc"),
                  y = unit(0.5, "npc"), gp = gpar(fontsize = 14, fontface = "bold"))

        upViewport()
    
    }
  
    dev.off()

    pdf(paste("scale_", filename, sep=''), width=6.83, height=4.83)
    grid.newpage()
    playout = grid.layout(1,1, widths=unit(c(6.83),c("inches")),
        heights=unit(c(4.35),c("inches")), respect = matrix(data =1, nrow=1, ncol=1))
    pushViewport(viewport(layout=playout))
    pushViewport(viewport(layout.pos.col=1, layout.pos.row=1))
    pushViewport(viewport(x=unit(0.5, "npc"), width=unit(0.5, "npc"), 
                          y=unit(0.5, "npc"), height=unit(0.25, "npc")))
    hmcols = colorRampPalette(colors.ramp)
    pal <- hmcols
    grid.imageFun(1, length(colors.ramp), pal(length(colors.ramp)))
    
#ths needs to be changed if "Error in validDetails.axis(x) :  'labels' and 'at' locations must have same length"
    if (length(seq(round(log(min.range, base = 2)), round(log(max.range, base =2)),by=round((round(log(max.range, base =2))-round(log(min.range, base =2)))/(length(colors.ramp)),1))) == length(colors.ramp)+1) {
        grid.xaxis(seq(0,1,by=1/(length(colors.ramp))), label=round(seq(round(log(min.range, base = 2)), round(log(max.range, base =2)),by=round((round(log(max.range, base =2))-round(log(min.range, base =2)))/(length(colors.ramp)),1)), 2), gp= gpar(cex = 1))
    } else {
        if (length(seq(round(log(min.range, base = 2)), round(log(max.range, base =2)),by=floor((round(log(max.range, base =2))-round(log(min.range, base =2)))/(length(colors.ramp))/0.1)*0.1)) == length(colors.ramp)+1) {
            grid.xaxis(seq(0,1,by=1/(length(colors.ramp))), label=round(seq(round(log(min.range, base = 2)), round(log(max.range, base =2)),by=floor((round(log(max.range, base =2))-round(log(min.range, base =2)))/(length(colors.ramp))/0.1)*0.1), 2), gp = gpar(cex = 1))
        }
    }
    
    grid.text(expression("log"[2]~"PRO Intensity"), y = unit(-3, "lines"))
    grid.rect()
    dev.off()
    
}

composites.func.panels.pro.2 <- function(dat, fact = 'RNA polymerase', summit = 'Summit', class= '', num=90, 
                                   col.lines = c(rgb(0,0,1,1/2), rgb(1,0,0,1/2),  rgb(0.1,0.5,0.05,1/2), rgb(0,0,0,1/2),  
                                       rgb(1/2,0,1/2,1/2), rgb(0,1/2,1/2,1/2), rgb(1/2,1/2,0,1/2)), fill.poly = c(rgb(0,0,1,1/4), 
              rgb(1,0,0,1/4), rgb(0.1,0.5,0.05,1/4),rgb(0,0,0,1/4), rgb(1/2,0,1/2,1/4))) {
    count = length(unique(dat$grp))
    ct.cons = 0
    lst.cons = list()
    unique(dat$grp)[order(unique(dat$grp))]
    #for (i in unique(dat$grp)[order(unique(dat$grp))]) {
    #    ct.cons= ct.cons + 1
    #  lst.cons[[ct.cons]] = c(min(dat[dat$grp == i,]$lower), max(dat[dat$grp == i,]$upper))
    #}
    pdf(paste('composite_', fact, '_signals_', summit, '_peaks', class, '.pdf', sep=''), width=6.83, 
        height=ceiling((count)) * 3.00) 
    print(xyplot(est ~ x|cond, group = grp, data = dat,
                 type = 'l',
               scales=list(x=list(cex=0.8,relation = "free"), y =list(cex=0.8, relation="free")),
                 xlim=c(-(num),(num)),
                 #ylim = lst.cons,
                 col = col.lines,
                 main=list(label=class, cex=0.6),
                 auto.key = list(points=F, lines=T, cex=0.8),
               par.settings = list(superpose.symbol = list(pch = c(16), col=col.lines, cex =0.5), 
                   superpose.line = list(col = col.lines, lwd=c(2,2,2,2,2,2), 
                       lty = c(1,1,1,1,1,1,1,1,1))),
                 cex.axis=1.0,
                 par.strip.text=list(cex=0.9, font=1, col='black'),
                 aspect=1.0,
               between=list(y=0.5, x=0.5),
                                        #lwd=2,
                 ylab = list(label = paste(fact," Density", sep=''), cex =0.8),
                 xlab = list(label = paste("Distance from ", summit, sep=''), cex =0.8),
                 upper = dat$upper,
                 fill = fill.poly,
                 lower = dat$lower,
                 strip = function(..., which.panel, bg) {
                     bg.col = c("grey85")
                 strip.default(..., which.panel = which.panel, bg = rep(bg.col, length = which.panel)[which.panel])
                 },
                 panel = function(x, y, ...){
                     #panel.superpose(x, y, panel.groups = 'my.panel.bands', ...)
                     panel.xyplot(x, y, ...)
       }
                 ))
    dev.off()
}

library(DESeq2)
library(lattice)
library(bigWig)
library(grid)
library(zoo)
library(latticeExtra)
library(data.table)


#processes

load('~/Desktop/fimo_stuff/final.dREG.df.Rdata')

hmap.chr.summit = sapply(strsplit(final.dREG.df$summit,":",fixed = TRUE), `[`, 1)
hmap.start.summit = sapply(strsplit((sapply(strsplit(final.dREG.df$summit,":",fixed = TRUE), `[`, 2)),"-",fixed = TRUE), `[`, 1)
hmap.end.summit = sapply(strsplit((sapply(strsplit(final.dREG.df$summit,":",fixed = TRUE), `[`, 2)),"-",fixed = TRUE), `[`, 2)

final.dREG.df$summitstart = as.numeric(hmap.start.summit)
final.dREG.df$summitend = as.numeric(hmap.end.summit)
final.dREG.df$range = final.dREG.df$end - final.dREG.df$start

#ad hoc filtering of very high signal 
#x = final.dREG.df[(final.dREG.df$min.0 < 200 &
#                              final.dREG.df$min.20 < 200 &
#                              final.dREG.df$min.40 < 200 &
#                              final.dREG.df$min.60 < 200 &
#                              final.dREG.df$min.120 < 200 &
#                              final.dREG.df$min.180 < 200 &
#                              final.dREG.df$min.240 < 200),] 

#x = final.dREG.df#[rev(order(final.dREG.df$range)),]
#this is just to keep the window on teh chromosome
x = final.dREG.df[final.dREG.df$summitstart > 5000 ,]

x = x[,c(1,12,13,2:3)]
x[,6] = '+'

                                        #where is the NA introduced?
#this helps make a window that is not off by one intended base
x[,3] = x[,3] -1

mjg = create.composites.heatmaps('/Users/guertinlab/Desktop/fimo_stuff/dREG_bigWigs/', x, region = 2000, step = 5, file.prefix = '3')

#this is the plus matrix:  y[[2]][[1]][[6]]
#this is the minus matrix: y[[2]][[2]][[6]]



draw.heatmap.both.pro.overlap(mjg[[2]], file = '3t3heatmap_MJG_full_final.pdf', upstream = -2000, downstream = 2000, sub.set = 30000, arrangement = 'pro', order.buffer = 25,
                              colors.ramp = c("#ffffff00","#FF0000FF", "#590000FF", "#510000FF", "#2F0000FF", "#2d0101FF","#2E0000FF"),
                              colors.ramp.2 = c("#ffffff00","#0000FF80", "#00105980", "#00115180", "#00072F80", "#00072E80", "#01072d80"),  width.vp = 2, height.vp = 4)

mjg2 = create.composites.heatmaps('/Users/guertinlab/Desktop/fimo_stuff/dREG_bigWigs/', x, region = 2000, step = 25, file.prefix = '3')

composites.func.panels.pro.2(mjg2[[1]], num = 2000, col.lines = c('blue', 'red'))
