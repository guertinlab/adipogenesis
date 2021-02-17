library(ggplot2)
library(lattice)
library(pheatmap)
library(fabricatr)
library(ggseqlogo)

load('/scratch/bhn9by/ATAC/plot.df.atac.Rdata')

dir = '/scratch/bhn9by/ATAC/fimo_composites/supp_figure_beds/'
setwd(dir)

#generate 'fimo.scores.atac.Rdata' object
res = unique(read.table('/scratch/bhn9by/ATAC/dynamic_peaks.bed'))
colnames(res) = c('chr', 'start', 'end')
res$start = as.numeric(as.character(res$start))
res$end = as.numeric(as.character(res$end))
rownames(res) = paste0(res[,1], ':', res[,2], '-', res[,3])

#supp figure families

for(bed.file in Sys.glob(file.path(paste0(dir,'*_fimo.bed')))) {
    
    factor.name = strsplit(bed.file, "/")[[1]]
    factor.name = strsplit(factor.name[length(factor.name)],
                           '_fimo.bed')[[1]][1]
    print(factor.name)
    x = read.table(bed.file, stringsAsFactors=FALSE)
    x = x[x[,6] != -1,]   
    y = aggregate(as.numeric(V8)~V1+V2+V3, data=x, FUN=sum)
    colnames(y) = c('chr', 'start', 'end', factor.name)
    rownames(y) = paste0(y[,1], ':', y[,2], '-', y[,3])
    
    func <- function(peak) {
        if(peak %in% rownames(y)) {
            return(y[rownames(y) == peak,ncol(y)])
        } else {
            return(NA)
        }
    }
    
    res[,ncol(res)+1] = sapply(rownames(res),func)
    colnames(res)[ncol(res)] = factor.name

}

res = res[,-c(1:3)]

fimo.scores.atac = res
save(fimo.scores.atac, file = '/scratch/bhn9by/ATAC/fimo_composites/supp_figure_beds/supp.fimo.scores.atac.Rdata')

#save as bed files
for(i in 1:ncol(fimo.scores.atac)) {
    temp = fimo.scores.atac[!is.na(fimo.scores.atac[,i]),]
    chr = sapply(strsplit(rownames(temp), ':'), '[', 1)
    x = sapply(strsplit(rownames(temp), ':'), '[', 2)
    start = sapply(strsplit(x, '-'), '[', 1)
    end = sapply(strsplit(x, '-'), '[', 2)
    bed = data.frame(chr, start, end, score = temp[,i])
    write.table(bed, file = paste0(colnames(fimo.scores.atac)[i],'_instances.bed'),row.names=F,col.names=F,quote=F,sep='\t')
}

#plot all family peaks
plot.df = data.frame()
for(i in 1:ncol(fimo.scores.atac)) {
    temp = plot.df.atac[plot.df.atac$genes %in% rownames(fimo.scores.atac[!is.na(fimo.scores.atac[,i]),]),]
    temp$family = colnames(fimo.scores.atac)[i]
    plot.df = rbind(plot.df,temp)
}

plot.df$genes = as.factor(plot.df$genes)
cat.colours = plot.df[plot.df$merge == 'one_groupt0', c(1,10)]
cat.colours$genes <- as.factor(cat.colours$genes)
cat.colours$supercluster <- as.factor(cat.colours$supercluster)

cat.colours$colour[cat.colours$supercluster == 'up.flat'] <- '#FF000008'
cat.colours$colour[cat.colours$supercluster == 'up.down'] <- '#FF000008'
cat.colours$colour[cat.colours$supercluster == 'gradual.up'] <- '#FF000008'
cat.colours$colour[cat.colours$supercluster == 'gradual.down'] <- '#0000FF08'
cat.colours$colour[cat.colours$supercluster == 'down.up'] <- '#0000FF08'

cat.colours$colour <- as.factor(cat.colours$colour)

cat.colours <- cat.colours[match(levels(plot.df$genes), cat.colours$genes), ]

pdf(file = 'all_families_dynamic_accessibility.pdf',width=14,height=14)

trellis.par.set(box.umbrella = list(lty = 1, col="black", lwd=1),
                box.rectangle = list( lwd=1.0, col="black", alpha = 1.0),
                plot.symbol = list(col="black", lwd=1.0, pch ='.'))

print(
    xyplot(value ~  sample.conditions | family, group = genes, data = plot.df, type = c('l'),#type = c('l','p'),
       scales=list(x=list(cex=1.0,relation = "free", rot = 45), y =list(cex=1.0, relation="free")),
       aspect=1.0,
       layout = c(3,2),
       between=list(y=0.5, x=0.5),
       index.cond=list(c(7:9,
                         4:6,
                         1:3)),
       main = list(label = 'All Family Motifs in Dynamic Peaks', cex = 1.5),
       ylab = list(label = 'Normalized ATAC signal', cex =1.0),
       xlab = list(label = 'Time (minutes)', cex =1.0),
       par.settings = list(superpose.symbol = list(pch = c(16),col=c('grey20'), cex =0.5),
                                   strip.background=list(col="grey80"),
                                   superpose.line = list(col = as.character(cat.colours$colour), lwd=c(1),lty = c(1))),
       panel = function(x, y, ...) {
           panel.xyplot(x, y, ...)
           #panel.bwplot(x, y, pch = '|', horizontal = FALSE, box.width = 15, do.out = FALSE)
           #panel.spline(x, y, col ='blue', lwd =2.0, ...) 
           
       })

)
dev.off()

#all families bar plot, not incl. no family
plot.df = data.frame()
for(i in 1:ncol(fimo.scores.atac)) {
    temp = plot.df.atac[plot.df.atac$genes %in% rownames(fimo.scores.atac[!is.na(fimo.scores.atac[,i]),]),]
    temp$family = colnames(fimo.scores.atac)[i]
    plot.df = rbind(plot.df,temp)
}

plot.df$status = 'Opened'
plot.df[plot.df$supercluster == 'gradual.down' | plot.df$supercluster == 'down.up',]$status = 'Closed'
activated = table(plot.df[plot.df$status == 'Opened',]$family) / 7
repressed = table(plot.df[plot.df$status == 'Closed',]$family) / 7
df = data.frame(names = c(names(activated),names(repressed)),
                num = c(as.vector(activated),as.vector(repressed)),
                cond = c(rep('Opened',length(activated)),rep('Closed',length(repressed))),
                index = rep(as.vector(table(plot.df$family)),2))

pdf(file='all.peaks.bar.pdf',width=7,height=5)
print(
        ggplot(df,aes(x = reorder(names,-index),y = num,fill=cond)) +
        geom_bar(stat='identity',position='stack',color='black') +
        #geom_text(aes(label=num),position=position_stack(vjust = 0.5),size=6) +
        labs(title = 'Number of Peaks For Each Motif Family', y = 'Number of Peaks', x = NULL, fill = NULL) +
        theme_minimal() +
        theme(panel.grid.minor = element_blank(),
              axis.ticks = element_blank(),
              axis.text.x = element_text(angle=45,size=20,hjust=.99,vjust=1,color='black',face='bold'),
              axis.text.y = element_text(size=20,face='bold',color='black'),
              axis.title.y = element_text(size=18,face='bold'),
              legend.text = element_text(size=18,face='bold'),
              plot.title = element_text(size=22,face='bold',hjust=0.5)) +
        scale_fill_manual(values = c("dodgerblue","indianred1"))
)
dev.off()

#all families bar plot, incl. no family
plot.df = data.frame()
for(i in 1:ncol(fimo.scores.atac)) {
    temp = plot.df.atac[plot.df.atac$genes %in% rownames(fimo.scores.atac[!is.na(fimo.scores.atac[,i]),]),]
    temp$family = colnames(fimo.scores.atac)[i]
    plot.df = rbind(plot.df,temp)
}

plot.df$genes = as.factor(plot.df$genes)
no.fam = rownames(fimo.scores.atac[which(rowSums(is.na(fimo.scores.atac)) == ncol(fimo.scores.atac)),])
temp = plot.df.atac[plot.df.atac$genes %in% no.fam,]
temp$family = 'No Family'
plot.df = rbind(plot.df,temp)

plot.df$status = 'Activated'
plot.df[plot.df$supercluster == 'gradual.down' | plot.df$supercluster == 'down.up',]$status = 'Repressed'
activated = table(plot.df[plot.df$status == 'Activated',]$family) / 7
repressed = table(plot.df[plot.df$status == 'Repressed',]$family) / 7
df = data.frame(names = c(names(activated),names(repressed)),
                num = c(as.vector(activated),as.vector(repressed)),
                cond = c(rep('Activated',length(activated)),rep('Repressed',length(repressed))),
                index = rep(as.vector(table(plot.df$family)),2))
df[df$names == 'No Family',]$index = min(df$index) - 1

pdf(file='all.peaks.with.no.fam.bar.pdf',width=12,height=7)
print(
        ggplot(df,aes(x = reorder(names,-index),y = num,fill=cond)) +
        geom_bar(stat='identity',position='stack',color='black') +
        #geom_text(aes(label=num),position=position_stack(vjust = 0.5),size=6) +
        labs(title = 'Number of Peaks For Each Motif Family', y = 'Number of Peaks', x = NULL, fill = NULL) +
        theme_minimal() +
        theme(panel.grid.minor = element_blank(),
              axis.ticks = element_blank(),
              axis.text.x = element_text(angle=45,size=20,hjust=.99,vjust=1,color='black',face='bold'),
              axis.text.y = element_text(size=20,face='bold',color='black'),
              axis.title.y = element_text(size=18,face='bold'),
              legend.text = element_text(size=18,face='bold'),
              plot.title = element_text(size=22,face='bold',hjust=0.5)) +
        scale_fill_manual(values = c("indianred1","dodgerblue"))
)
dev.off()

#now isolated peaks
plot.df = data.frame()

for(i in 1:ncol(fimo.scores.atac)) {
      scores.temp = fimo.scores.atac[!is.na(fimo.scores.atac[,i]),]
      scores.temp = scores.temp[,-i]
      scores.temp = scores.temp[which(rowSums(is.na(scores.temp)) == ncol(scores.temp)),]
      temp = plot.df.atac[plot.df.atac$genes %in% rownames(scores.temp),]
      temp$family = colnames(fimo.scores.atac)[i]
      plot.df = rbind(plot.df,temp)
}

plot.df$genes = as.factor(plot.df$genes)
cat.colours = plot.df[plot.df$merge == 'one_groupt0', c(1,10)]
cat.colours$genes <- as.factor(cat.colours$genes)
cat.colours$supercluster <- as.factor(cat.colours$supercluster)

cat.colours$colour[cat.colours$supercluster == 'up.flat'] <- '#FF00001A'
cat.colours$colour[cat.colours$supercluster == 'up.down'] <- '#FF00001A'
cat.colours$colour[cat.colours$supercluster == 'gradual.up'] <- '#FF00001A'
cat.colours$colour[cat.colours$supercluster == 'gradual.down'] <- '#0000FF1A'
cat.colours$colour[cat.colours$supercluster == 'down.up'] <- '#0000FF1A'

cat.colours$colour <- as.factor(cat.colours$colour)

cat.colours <- cat.colours[match(levels(plot.df$genes), cat.colours$genes), ]

pdf(file = 'all_isolated_families_dynamic_accessibility.pdf',width=14,height=14)

trellis.par.set(box.umbrella = list(lty = 1, col="black", lwd=1),
                box.rectangle = list( lwd=1.0, col="black", alpha = 1.0),
                plot.symbol = list(col="black", lwd=1.0, pch ='.'))

print(
    xyplot(value ~  sample.conditions | family, group = genes, data = plot.df, type = c('l'),#type = c('l','p'),
       scales=list(x=list(cex=1.5,relation = "free", rot = 45), y =list(cex=1.5, relation="free")),
       aspect=1.0,
       layout = c(3,2),
       between=list(y=0.5, x=0.5),
       index.cond=list(c(7:9,
                         4:6,
                         1:3)),
       main = list(label = 'All Isolated Family Motifs in Dynamic Peaks', cex = 2.0),
       ylab = list(label = 'Normalized ATAC signal', cex =2.0),
       xlab = list(label = 'Time (minutes)', cex =2.0),
       par.settings = list(superpose.symbol = list(pch = c(16),col=c('grey20'), cex = 2.0),
                           strip.background=list(col="grey80"),
                           superpose.line = list(col = as.character(cat.colours$colour), lwd=c(1),lty = c(1))),
       panel = function(x, y, ...) {
           panel.xyplot(x, y, ...)
           #panel.bwplot(x, y, pch = '|', horizontal = FALSE, box.width = 15, do.out = FALSE)
           #panel.spline(x, y, col ='blue', lwd =2.0, ...) 
           
       })

)
dev.off()

#isolated peaks bar plot
plot.df$status = 'Activated'
plot.df[plot.df$supercluster == 'gradual.down' | plot.df$supercluster == 'down.up',]$status = 'Repressed'
activated = table(plot.df[plot.df$status == 'Activated',]$family) / 7
repressed = table(plot.df[plot.df$status == 'Repressed',]$family) / 7
df = data.frame(names = c(names(activated),names(repressed)),
                num = c(as.vector(activated),as.vector(repressed)),
                cond = c(rep('Activated',length(activated)),rep('Repressed',length(repressed))),
                index = rep(as.vector(table(plot.df$family)),2))

pdf(file='isolated.peaks.bar.pdf',width=12,height=7)
print(
        ggplot(df,aes(x = reorder(names,-index),y = num,fill=cond)) +
        geom_bar(stat='identity',position='stack',color='black') +
        #geom_text(aes(label=num),position=position_stack(vjust = 0.5),size=6) +
        labs(title = 'Number of Isolated Peaks For Each Motif Family',y = 'Number of Peaks', x = NULL, fill = NULL) +
        theme_minimal() +
        theme(panel.grid.minor = element_blank(),
              axis.ticks = element_blank(),
              axis.text.x = element_text(angle=45,size=12,hjust=.99,vjust=1,color='black',face='bold'),
              axis.text.y = element_text(size=12,face='bold',color='black'),
              axis.title.y = element_text(size=14,face='bold'),
              legend.text = element_text(size=12,face='bold'),
              plot.title = element_text(size=16,face='bold',hjust=0.5)) +
        scale_fill_manual(values = c("indianred1","dodgerblue"))
)
dev.off()

#heatmaps
prop.fam.df = data.frame(matrix(nrow=ncol(fimo.scores.atac),ncol=ncol(fimo.scores.atac)),row.names=colnames(fimo.scores.atac))
total.fam.df = data.frame(matrix(nrow=ncol(fimo.scores.atac),ncol=ncol(fimo.scores.atac)),row.names=colnames(fimo.scores.atac))
colnames(prop.fam.df) = colnames(fimo.scores.atac)
colnames(total.fam.df) = colnames=colnames(fimo.scores.atac)

for(i in 1:ncol(fimo.scores.atac)) {

    x = fimo.scores.atac[!is.na(fimo.scores.atac[,i]),]
    
    for(j in 1:ncol(x)) {
        z = x[!is.na(x[,j]),]
        total.fam.df[j,i] = nrow(z)
        prop.fam.df[j,i] = 100*(nrow(z)/nrow(fimo.scores.atac[!is.na(fimo.scores.atac[,j]),]))  
    }
}

rownames(prop.fam.df) = colnames(fimo.scores.atac)
colnames(prop.fam.df) = colnames(fimo.scores.atac)

prop.fam.df = prop.fam.df[order(colnames(prop.fam.df)),]

pdf(file='prop.fam.heatmap.pdf')
pheatmap(prop.fam.df,cluster_rows=FALSE, show_rownames=TRUE,cluster_cols=FALSE,fontsize=12)
dev.off()

pdf(file='total.fam.heatmap.pdf')
pheatmap(total.fam.df,cluster_rows=FALSE, show_rownames=TRUE,cluster_cols=FALSE,fontsize=12)
dev.off()

#up/down split occurrences
plot.df = plot.df.atac
plot.df$status = 'Activated'
plot.df[plot.df$supercluster == 'gradual.down' | plot.df$supercluster == 'down.up',]$status = 'Repressed'

for(i in 1:ncol(fimo.scores.atac)) {
    up.num = c()
    down.num = c()
    
    for(j in 1:ncol(fimo.scores.atac)) {
        
        x = fimo.scores.atac[!is.na(fimo.scores.atac[,j]),]
        y = plot.df[plot.df$genes %in% rownames(x),]
        up = fimo.scores.atac[rownames(fimo.scores.atac) %in% y[y$status == 'Activated',]$genes,]
        down = fimo.scores.atac[rownames(fimo.scores.atac) %in% y[y$status == 'Repressed',]$genes,]

        up.num = append(up.num, (nrow(up[!is.na(up[,i]),])/nrow(up))*100)
        down.num = append(down.num, (nrow(down[!is.na(down[,i]),])/nrow(down))*100) 
    }

    df = data.frame(names = rep(colnames(fimo.scores.atac),2),
                    #names = rep(c('NRF','TWIST','STAT','TFAP2','ZNF263','AP1','bHLH','GR','CTCFL','SP/KLF','ELK/ETV','GFX/ZBTB33','Maz','NFY'),2),
                    num = c(up.num,down.num),
                    cond = c(rep('Activated',length(up.num)),rep('Repressed',length(down.num))))
                    #index = rep(c(10:14,1:9),2))
    
    name = df[i,]$names
    
    df = df[-which(df$names == name),]

    df$names = as.factor(df$names)
    
    pdf(file=paste0('percent.',name,'.in.other.families.pdf'),width=8,height=5)
    print(
        ggplot(df,aes(x = names,y = num,fill=cond)) +
        geom_bar(stat='identity',position='dodge',color='black') + 
        labs(title = paste0('% of ',name,' Occurrences in Other Family Motifs'),
             y = paste0('% of Occurrences with ',name),
             x = 'Motif Family',
             fill = 'Direction') +
        theme_minimal() +
        scale_x_discrete(labels= levels(df$names)) +
        theme(panel.grid.minor = element_blank(),
              plot.title = element_text(hjust = 0.5,face='bold'),   
              axis.text.x = element_text(angle=45,size=12,hjust=.99,vjust=1,color='black',face='bold'),
              axis.title.x = element_text(size=14,face='bold'),
              axis.text.y = element_text(size=12,face='bold',color='black'),
              axis.title.y = element_text(size=14,face='bold'),
              legend.title = element_text(size=14,face='bold',color='black'),
              legend.text = element_text(size=12,face='bold',color='black'),
              axis.ticks = element_blank()) +
        scale_fill_manual(values = c("indianred1","dodgerblue"))
    )
    dev.off() 
}

#fimo scores box and whisker plots - isolated peaks
plot.df = data.frame()

for(i in 1:ncol(fimo.scores.atac)) {
      scores.temp = fimo.scores.atac[!is.na(fimo.scores.atac[,i]),]
      scores.temp = scores.temp[,-i]
      scores.temp = scores.temp[which(rowSums(is.na(scores.temp)) == ncol(scores.temp)),]
      temp = plot.df.atac[plot.df.atac$genes %in% rownames(scores.temp),]
      temp$family = colnames(fimo.scores.atac)[i]
      plot.df = rbind(plot.df,temp)
}

plot.df$status = 'Activated'
plot.df[plot.df$supercluster == 'gradual.down' | plot.df$supercluster == 'down.up',]$status = 'Repressed'
plot.df = unique(plot.df[,c(1,11,12)])

func <- function(peak) {
    family = plot.df[plot.df$genes == peak,]$family    
    colnum = which(colnames(fimo.scores.atac) == family)
    score = fimo.scores.atac[rownames(fimo.scores.atac) == peak,colnum]
    return(score)
}

plot.df$score = sapply(plot.df$genes,func)

pdf(file = paste0('fimo.scores.isolated.peaks.bw.pdf'),width=12,height=8)

trellis.par.set(box.umbrella = list(lty = 1, col=c("red", "blue"), lwd=2),
                box.rectangle = list( lwd=2.0, col=c("red", "blue"), alpha = 1.0),
                plot.symbol = list(col=c("red", "blue"), lwd=2.0, pch ='.'))

print(
    bwplot(log(as.numeric(as.character(score)), base = 10) ~ status | family, data = plot.df, horizontal = FALSE, pch = '|', do.out = FALSE,
           scales=list(x=list(cex=1.0, relation = "free", rot = 45), y = list(cex=1.0, relation="free")),
           aspect=2.0,
           between=list(y=0.5, x=0.5),
           index.cond=list(c(7:9,4:6,1:3)),
           ylab = expression('log'[10]* paste(Sigma,'(motif scores)')),
           xlab = expression('ATAC Peak category'),
                                        #manually settign to avoid outlier, since do.out = FALSE       
           #ylim = list(c(0.9, 1.9),  c(0.9, 1.7), c(0.8, 1.8), c(1, 2.15),  c(0.96, 1.4)),
           par.settings = list(superpose.symbol = list(pch = c(16),col=c('grey20'), cex =0.5),
                                        #I want to change the background strip to the corresponding motif color
                               strip.background=list(col=c("grey80"))))
    
)
dev.off()

#fimo scores box and whisker plots - all peaks
plot.df = data.frame()
for(i in 1:ncol(fimo.scores.atac)) {
    temp = plot.df.atac[plot.df.atac$genes %in% rownames(fimo.scores.atac[!is.na(fimo.scores.atac[,i]),]),]
    temp$family = colnames(fimo.scores.atac)[i]
    plot.df = rbind(plot.df,temp)
}

plot.df$status = 'Activated'
plot.df[plot.df$supercluster == 'gradual.down' | plot.df$supercluster == 'down.up',]$status = 'Repressed'
plot.df = unique(plot.df[,c(1,11,12)])

func <- function(peak,family) {
    #family = plot.df[plot.df$genes == peak,]$family
    colnum = which(colnames(fimo.scores.atac) == family)
    score = fimo.scores.atac[rownames(fimo.scores.atac) == peak,colnum]
    return(score)
}

plot.df$score = mapply(func,plot.df$genes,family=plot.df$family)

pdf(file = paste0('fimo.scores.all.peaks.bw.pdf'),width=12,height=8)

trellis.par.set(box.umbrella = list(lty = 1, col=c("red", "blue"), lwd=2),
                box.rectangle = list( lwd=2.0, col=c("red", "blue"), alpha = 1.0),
                plot.symbol = list(col=c("red", "blue"), lwd=2.0, pch ='.'))

print(
    bwplot(log(as.numeric(as.character(score)), base = 10) ~ status | family, data = plot.df, horizontal = FALSE, pch = '|', do.out = FALSE,
           scales=list(x=list(cex=1.0, relation = "free", rot = 45), y = list(cex=1.0, relation="free")),
           aspect=2.0,
           between=list(y=0.5, x=0.5),
           index.cond=list(c(7:9,4:6,1:3)),
           ylab = expression('log'[10]* paste(Sigma,'(motif scores)')),
           xlab = expression('ATAC Peak category'),
                                        #manually settign to avoid outlier, since do.out = FALSE       
           #ylim = list(c(0.9, 1.9),  c(0.9, 1.7), c(0.8, 1.8), c(1, 2.15),  c(0.96, 1.4)),
           par.settings = list(superpose.symbol = list(pch = c(16),col=c('grey20'), cex =0.5),
                                        #I want to change the background strip to the corresponding motif color
                               strip.background=list(col=c("grey80"))))
    
)
dev.off()

#split on fimo scores - isolated peaks
plot.df = data.frame()

for(i in 1:ncol(fimo.scores.atac)) {
      scores.temp = fimo.scores.atac[!is.na(fimo.scores.atac[,i]),]
      scores.temp = scores.temp[,-i]
      scores.temp = scores.temp[which(rowSums(is.na(scores.temp)) == ncol(scores.temp)),]
      temp = plot.df.atac[plot.df.atac$genes %in% rownames(scores.temp),]
      temp$family = colnames(fimo.scores.atac)[i]
      plot.df = rbind(plot.df,temp)
}

plot.df$status = 'Activated'
plot.df[plot.df$supercluster == 'gradual.down' | plot.df$supercluster == 'down.up',]$status = 'Repressed'

func <- function(peak) {
    family = unique(plot.df[plot.df$genes == peak,]$family)  
    colnum = which(colnames(fimo.scores.atac) == family)
    score = fimo.scores.atac[rownames(fimo.scores.atac) == peak,colnum]
    return(score)
}

plot.df$score = sapply(plot.df$genes,func)

for(fam in unique(plot.df$family)) {
    df = plot.df[plot.df$family == fam,]

    fam = gsub('/','.',fam)
    print(fam)
    
    df$quantile = paste0('Quantile ',split_quantile(x = as.numeric(df$score),type = 5))
    
    df$genes = as.factor(df$genes)
    cat.colours = df
    cat.colours$genes <- as.factor(cat.colours$genes)
    cat.colours$status <- as.factor(cat.colours$status)
    
    cat.colours$colour[cat.colours$status == 'Activated'] <- '#FF000008'
    cat.colours$colour[cat.colours$status == 'Repressed'] <- '#0000FF08'
    
    cat.colours$colour <- as.factor(cat.colours$colour)
    
    cat.colours <- cat.colours[match(levels(df$genes), cat.colours$genes),]
    
    pdf(file = paste0(fam,'_isolated_split_on_FIMO_scores.pdf'),width=11,height=4)
    
    trellis.par.set(box.umbrella = list(lty = 1, col="black", lwd=1),
                    box.rectangle = list( lwd=1.0, col="black", alpha = 1.0),
                    plot.symbol = list(col="black", lwd=1.0, pch ='.'))
    
    print(
        xyplot(value ~  sample.conditions | quantile, group = genes, data = df, type = c('l'),#type = c('l','p'),
               scales=list(x=list(cex=1.0,relation = "free", rot = 45), y = list(cex=1.0, relation="free")),
               aspect=1.0,
               layout = c(5,1),
               between=list(y=0.5, x=0.5),
               main = list(label = paste0(fam, ' Isolated Motifs Stratified By FIMO Score'), cex = 1.5),
               ylab = list(label = 'Normalized ATAC signal', cex =1.0),
               xlab = list(label = 'Time (minutes)', cex =1.0),
               par.settings = list(superpose.symbol = list(pch = c(16),col=c('grey20'), cex =0.5),
                                   strip.background=list(col="grey80"),
                                   superpose.line = list(col = as.character(cat.colours$colour), lwd=c(1),lty = c(1))),
               panel = function(x, y, ...) {
                   panel.xyplot(x, y, ...)
                   panel.bwplot(x, y, pch = '|', horizontal = FALSE, box.width = 15, do.out = FALSE)
                   panel.spline(x, y, col = 'grey70', lwd =3.0, ...) 
                   
               })
    )
    
    dev.off()
    
}

#split on fimo scores - all peaks
plot.df = data.frame()

for(i in 1:ncol(fimo.scores.atac)) {
    temp = plot.df.atac[plot.df.atac$genes %in% rownames(fimo.scores.atac[!is.na(fimo.scores.atac[,i]),]),]
    temp$family = colnames(fimo.scores.atac)[i]
    plot.df = rbind(plot.df,temp)
}

plot.df$status = 'Activated'
plot.df[plot.df$supercluster == 'gradual.down' | plot.df$supercluster == 'down.up',]$status = 'Repressed'

func <- function(peak,family) {
    #family = plot.df[plot.df$genes == peak,]$family
    colnum = which(colnames(fimo.scores.atac) == family)
    score = fimo.scores.atac[rownames(fimo.scores.atac) == peak,colnum]
    return(score)
}

plot.df$score = mapply(func,plot.df$genes,family=plot.df$family)

for(fam in unique(plot.df$family)) {
    df = plot.df[plot.df$family == fam,]

    print(fam)
    
    df$quantile = paste0('Quantile ',split_quantile(x = as.numeric(df$score),type = 5))
    
    df$genes = as.factor(df$genes)
    cat.colours = df
    cat.colours$genes <- as.factor(cat.colours$genes)
    cat.colours$status <- as.factor(cat.colours$status)
    
    cat.colours$colour[cat.colours$status == 'Activated'] <- '#FF000008'
    cat.colours$colour[cat.colours$status == 'Repressed'] <- '#0000FF08'
    
    cat.colours$colour <- as.factor(cat.colours$colour)
    
    cat.colours <- cat.colours[match(levels(df$genes), cat.colours$genes),]
    
    pdf(file = paste0(fam,'_all_split_on_FIMO_scores.pdf'),width=11,height=4)
    
    trellis.par.set(box.umbrella = list(lty = 1, col="black", lwd=1),
                    box.rectangle = list( lwd=1.0, col="black", alpha = 1.0),
                    plot.symbol = list(col="black", lwd=1.0, pch ='.'))
    
    print(
        xyplot(value ~  sample.conditions | quantile, group = genes, data = df, type = c('l'),#type = c('l','p'),
               scales=list(x=list(cex=1.0,relation = "free", rot = 45), y = list(cex=1.0, relation="free")),
               aspect=1.0,
               layout = c(5,1),
               between=list(y=0.5, x=0.5),
               main = list(label = paste0(fam, ' All Motifs Stratified By FIMO Score'), cex = 1.5),
               ylab = list(label = 'Normalized ATAC signal', cex =1.0),
               xlab = list(label = 'Time (minutes)', cex =1.0),
               par.settings = list(superpose.symbol = list(pch = c(16),col=c('grey20'), cex =0.5),
                                   strip.background=list(col="grey80"),
                                   superpose.line = list(col = as.character(cat.colours$colour), lwd=c(1),lty = c(1))),
               panel = function(x, y, ...) {
                   panel.xyplot(x, y, ...)
                   panel.bwplot(x, y, pch = '|', horizontal = FALSE, box.width = 15, do.out = FALSE)
                   panel.spline(x, y, col = 'grey70', lwd =3.0, ...) 
                   
               })
    )
    
    dev.off()
    
}

#get composites for activated / repressed (all peaks)
pswm.fimo <- function(fimo.in, out = 'outfilename', nm = 'bHLH_Activated', rc = FALSE) {
    posnum = nchar(fimo.in$sequence[1])
    fimo.in$sequence = toupper(fimo.in$sequence)
    col.matrix = matrix()
    for (g in 1:posnum){
        itnum = lapply(strsplit(as.character(fimo.in$sequence), ''), "[", g)
        if (g == 1) {
            col.matrix = itnum
        } else {
            col.matrix = cbind(col.matrix, itnum)
        }
    }  
  a.nuc = sapply(1:posnum, function(x) sum(col.matrix[,x] == "A"))
  t.nuc = sapply(1:posnum, function(x) sum(col.matrix[,x] == "T"))
  c.nuc = sapply(1:posnum, function(x) sum(col.matrix[,x] == "C"))
  g.nuc = sapply(1:posnum, function(x) sum(col.matrix[,x] == "G"))
  
  pswm = cbind(a.nuc, c.nuc, g.nuc, t.nuc)
  print(pswm)
  outfile = file(paste0(out, '.txt'))
  on.exit(close(outfile))
  writeLines(c("MEME version 4", "ALPHABET= ACGT", "strands: + -", " ", 
               "Background letter frequencies (from uniform background):", 
               "A 0.30000 C 0.20000 G 0.20000 T 0.30000", paste("MOTIF", out), " ",
               paste("letter-probability matrix: alength= 4 w=", posnum)), outfile)
  pswm = pswm/rowSums(pswm)
  if (rc == "TRUE"){
    pswm<- pswm[nrow(pswm):1,ncol(pswm):1]
  } else {}   
    write.table(pswm, file = paste0(out, '.txt'), append = TRUE, quote=FALSE, row.names =FALSE, col.names = FALSE)
    pswm = t(pswm)
    rownames(pswm) = c('A', 'C', 'G', 'T')
    return(pswm)
#    
#    system(paste0('/Users/guertinlab/meme/libexec/meme-5.1.1/ceqlogo -i ', out, '.txt -m Composite -o ', nm, '.eps'))
#    system(paste0('/Users/guertinlab/meme/libexec/meme-5.1.1/ceqlogo -i ', out, '.txt -m Composite -o ', nm, '.rc.eps -r'))
}

plot.df = plot.df.atac
plot.df$status = 'Activated'
plot.df[plot.df$supercluster == 'gradual.down' | plot.df$supercluster == 'down.up',]$status = 'Repressed'

count = -1
vec.names = c()

df.seq = as.data.frame(matrix(ncol=5, nrow=0),stringsAsFactors = FALSE)
colnames(df.seq) = c('chr', 'start', 'end', 'sequence', 'factor')

for(bed.file in Sys.glob(file.path(paste0(dir,'main_figure_beds/*fimo.bed')))) {
    #print(bed.file)
    count = count + 1
    factor.name = strsplit(bed.file, "/")[[1]]
    factor.name = strsplit(factor.name[length(factor.name)],
                           '_fimo.bed')[[1]][1]
    print(factor.name)
    x = read.table(bed.file, stringsAsFactors=FALSE)
    x = x[x[,6] != -1,]
    x = x[,c(1,2,3,10)]
    x$factor = factor.name
    colnames(x) = c('chr', 'start', 'end', 'sequence', 'factor')
    x$re = paste0(x[,1], ':', x[,2], '-', x[,3])
    df.seq = rbind(df.seq, x) 
}

#now add SP and KLF
sp = read.table('/scratch/bhn9by/ATAC/SP_KLF_split/output_sp1.txt', stringsAsFactors=FALSE, sep = '\t', header = TRUE)[,c(3,10)]

out.sp <- strsplit(as.character(sp[,1]), ':') 
sp <- data.frame(do.call(rbind, out.sp, quote=FALSE), sp[,c(2)])
sp <- sp[,c(1:3,10)]
colnames(sp) <- c('chr','start','end','sequence')
sp$factor = 'SP'
sp$re = paste0(sp$chr,':',sp$start,'-',sp$end)

df.seq = rbind(df.seq,sp)

klf = read.table('/scratch/bhn9by/ATAC/SP_KLF_split/output_klf.txt', stringsAsFactors=FALSE, sep = '\t', header = TRUE)[,c(3,10)]

out.klf <- strsplit(as.character(klf[,1]), ':') 
klf <- data.frame(do.call(rbind, out.klf, quote=FALSE), klf[,c(2)])
klf <- klf[,c(1:3,10)]
colnames(klf) <- c('chr','start','end','sequence')
klf$factor = 'KLF'
klf$re = paste0(klf$chr,':',klf$start,'-',klf$end)

df.seq = rbind(df.seq,klf)

func <- function(peak) {
     return(unique(plot.df[plot.df$genes == peak,]$status))
}

df.seq$status = sapply(df.seq$re,func)

for(i in unique(df.seq$factor)) {
      act = df.seq[df.seq$factor == i & df.seq$status == 'Activated',]
      rep = df.seq[df.seq$factor == i & df.seq$status == 'Repressed',]
      pswm.fimo(act, out = paste0(i,'_activated'))
      pswm.fimo(rep, out = paste0(i,'_repressed'))
}     
                                        

