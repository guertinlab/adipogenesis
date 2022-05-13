library(ggplot2)

setwd('~/Adipogenesis/ChIPseq')

load('~/Adipogenesis/PRO_analysis_redo/trans.links.210923.Rdata')
load('~/Adipogenesis/ATAC_analysis_redo/final.atac.all.df.Rdata')

factor.df <- data.frame(file.name = c('CEBPbeta','GR','JunB','KLF4','KLF5','Tw2_GM','cJun'),
                        factor.name = c('CEBP','GR','AP1','KLF','KLF','TWIST','AP1'))

plot.df <- data.frame()

for(bed in Sys.glob('*.atac.peaks.overlap.bed')) {
    name <- strsplit(bed,'[.]')[[1]][1]
    print(name)
    
    x <- unique(read.table(bed))
    x$peak <- paste0(x$V1,':',x$V2,'-',x$V3)
    x$RE <- paste0(x$V7,':',x$V8,'-',x$V9)
    x[x$V8 == -1,]$RE <- NA
    x$w.motif <- FALSE
    x[x$V6 != 0,]$w.motif <- TRUE

    factor <- factor.df[factor.df$file.name == name,'factor.name']
    linked.REs <- unique(trans.links[trans.links$source == factor,'target'])

    x$linked <- FALSE
    x[x$RE %in% linked.REs,]$linked <- TRUE

    pdf(file = paste0('ChIP.peak.score.v.num.',name,'.motifs.pdf'),width=6,height=4)
    print(
        ggplot(x, aes(x=factor(V6), y=log(V5,10))) + 
        geom_boxplot(color="black",outlier.shape=NA) +
#geom_violin(alpha=0.5) +
        geom_jitter(shape=16, position=position_jitter(0.2),alpha=0.1) +
        theme_minimal() +
        labs(y = 'log10 Peak Score', 
             x = 'Num. Motifs',
             title= paste0('ChIP Peak Score v. Num. ',name,' Motifs')) +
        theme(
#panel.grid.minor = element_blank(),
#legend.position = "none",
            plot.title = element_text(size=14,face='bold',hjust = 0.5),
            axis.ticks = element_blank(),
            axis.text.y = element_text(size=11,color='black',face='bold'),
            axis.text.x = element_text(size=11,color='black',angle = 45,face='bold',vjust=0.95,hjust=0.9),
            axis.title = element_text(size=12,color='black',face='bold'),
            legend.text = element_text(size=10,color='black',face='bold'),
            legend.title = element_text(size=11,color='black',face='bold'))
#coord_cartesian(ylim = c(-10, 10))
    )
    dev.off()
    
    x$cat <- 'Overlap w/ Non-linked RE'
    x[x$linked == TRUE,]$cat <- 'Overlap w/ Linked RE'
    x[is.na(x$RE),]$cat <- 'No Overlap w/ RE'
  
    x$factor <- name

    plot.df <- rbind(plot.df,x)
    
}

plot.df$factor <- factor(plot.df$factor,levels=c('JunB','cJun','CEBPbeta','GR','KLF4','KLF5','Tw2_GM'))

plot.df$cat <- factor(plot.df$cat,levels=c("No Overlap w/ RE",
                                           "Overlap w/ Non-linked RE",
                                           "Overlap w/ Linked RE"))

pdf(file = 'ChIP.peak.score.v.peak.cat.pdf',width=9,height=4)
print(
    ggplot(plot.df[plot.df$w.motif == TRUE,], aes(x=factor, y=log(V5,10),fill=cat)) + 
    geom_boxplot(color="black") +
#geom_violin(alpha=0.5) +
    #geom_jitter(shape=16, position=position_jitter(0.2),alpha=0.2) +
    theme_minimal() +
    labs(y = 'ChIP Peak Score', 
         x = 'Factor',
         title= 'All ChIP Peaks w/ Motif',
         fill = NULL) +
    theme(
#panel.grid.minor = element_blank(),
#legend.position = "none",
        plot.title = element_text(size=14,face='bold',hjust = 0.5),
        axis.ticks = element_blank(),
        axis.text.y = element_text(size=11,color='black',face='bold'),
        axis.text.x = element_text(size=11,color='black',angle = 45,face='bold',vjust=0.95,hjust=0.9),
        axis.title = element_text(size=12,color='black',face='bold'),
        legend.text = element_text(size=10,color='black',face='bold'),
        legend.title = element_text(size=11,color='black',face='bold'))
#coord_cartesian(ylim = c(-10, 10))
)
dev.off()

#Num. motifs in ChIP peak v. whether overlapping RE is linked or not

pdf(file = 'num.motifs.in.ChIP.peak.v.peak.cat.pdf',width=6,height=4)
print(
    ggplot(plot.df[plot.df$w.motif == TRUE,], aes(x=factor, y=V6,fill=cat)) + 
    geom_boxplot(color="black") +
#geom_violin(alpha=0.5) +
    #geom_jitter(shape=16, position=position_jitter(0.2),alpha=0.2) +
    theme_minimal() +
    labs(y = 'Num. Motifs', 
         x = 'Factor',
         title= NULL,
         fill = NULL) +
    theme(
#panel.grid.minor = element_blank(),
#legend.position = "none",
        plot.title = element_text(size=14,face='bold',hjust = 0.5),
        axis.ticks = element_blank(),
        axis.text.y = element_text(size=11,color='black',face='bold'),
        axis.text.x = element_text(size=11,color='black',angle = 45,face='bold',vjust=0.95,hjust=0.9),
        axis.title = element_text(size=12,color='black',face='bold'),
        legend.text = element_text(size=10,color='black',face='bold'),
        legend.title = element_text(size=11,color='black',face='bold'))
#coord_cartesian(ylim = c(-10, 10))
)
dev.off()

#num. peaks w/ and w/out motif for each peak set
bar.plot.df <- data.frame(file.name = rep(c('CEBPbeta','GR','JunB','KLF4','KLF5','Tw2_GM','cJun'),each = 2),
                      factor.name = rep(c('CEBP','GR','AP1','KLF','KLF','TWIST','AP1'),each=2),
                      type = rep(c('w.motif','wout.motif'),7),
                      num=0)

for(factor in unique(bar.plot.df$file.name)) {
    bar.plot.df[bar.plot.df$file.name == factor & bar.plot.df$type == 'w.motif','num'] <-
        length(unique(plot.df[plot.df$factor == factor & plot.df$w.motif == TRUE,'peak']))
    bar.plot.df[bar.plot.df$file.name == factor & bar.plot.df$type == 'wout.motif','num'] <-
        length(unique(plot.df[plot.df$factor == factor & plot.df$w.motif == FALSE,'peak']))
}

bar.plot.df$file.name <- factor(bar.plot.df$file.name,levels=c('JunB','cJun','CEBPbeta','GR','KLF4','KLF5','Tw2_GM'))

pdf(file='ChIP.peaks.w.motif.bar.pdf',height=5,width=6) 
print(
    ggplot(bar.plot.df,aes(x = file.name,y = num,fill=type)) +
    geom_bar(stat='Identity',color='black',position='stack') +
    labs(title = 'Num. ChIP-seq Peaks w/ Motif', 
         y = 'Num. Peaks', 
         x = 'ChIP-seq Factor',
         fill = NULL) +
    theme_minimal() +
    theme(panel.grid.minor = element_blank(),
          axis.ticks = element_blank(),
          axis.text.x = element_text(angle=45,size=12,hjust=.99,vjust=1,color='black',face='bold'),
          axis.text.y = element_text(size=14,face='bold',color='black'),
          axis.title.x = element_text(size=14,face='bold'),
          axis.title.y = element_text(size=14,face='bold'),
          legend.title = element_text(size=12,face='bold'),
          legend.text = element_text(size=10,face='bold'),
          plot.title = element_text(size=16,face='bold',hjust=0.5)) +
    scale_fill_manual(values = c('dodgerblue','lightgrey'))
)
dev.off()

#ChIP peak score v. motif presence

pdf(file = 'ChIP.peak.score.v.motif.presence.pdf',width=6,height=4)
print(
    ggplot(plot.df, aes(x=factor, y=log(V5,10),fill=w.motif)) + 
    geom_boxplot(color="black") +
#geom_violin(alpha=0.5) +
#geom_jitter(aes(col = Sample), shape=16, position=position_jitter(0.2)) +
    theme_minimal() +
    labs(y = 'log10 Peak Score', 
         x = 'Factor',
         title= 'ChIP Peak Score v. Motif Presence',
         fill = 'w/ motif') +
    theme(
#panel.grid.minor = element_blank(),
#legend.position = "none",
        plot.title = element_text(size=14,face='bold',hjust = 0.5),
        axis.ticks = element_blank(),
        axis.text.y = element_text(size=11,color='black',face='bold'),
        axis.text.x = element_text(size=11,color='black',angle = 45,face='bold',vjust=0.95,hjust=0.9),
        axis.title = element_text(size=12,color='black',face='bold'),
        legend.text = element_text(size=10,color='black',face='bold'),
        legend.title = element_text(size=11,color='black',face='bold'))
#coord_cartesian(ylim = c(-10, 10))
)
dev.off()

#Perc. Peaks Overlapping REs
bar.plot.df <- data.frame(file.name = rep(c('CEBPbeta','GR','JunB','KLF4','KLF5','Tw2_GM','cJun'),each = 2),
                      factor.name = rep(c('CEBP','GR','AP1','KLF','KLF','TWIST','AP1'),each=2),
                      type = rep(c('linked','not.linked'),7),
                      perc=0)

for(factor in unique(bar.plot.df$file.name)) {
    temp.df <- plot.df[plot.df$factor == factor & plot.df$V6 != 0,]

    total.peaks <- length(unique(temp.df$peak))

    not.linked <- length(unique(temp.df[!is.na(temp.df$RE) & temp.df$linked == FALSE,'peak']))
    linked <- length(unique(temp.df[!is.na(temp.df$RE) & temp.df$linked == TRUE,'peak']))

    bar.plot.df[bar.plot.df$file.name == factor & bar.plot.df$type == 'linked','perc'] <- 100 * (linked/total.peaks)
    bar.plot.df[bar.plot.df$file.name == factor & bar.plot.df$type == 'not.linked','perc'] <- 100 * (not.linked/total.peaks)
    
    #prop.table(table(temp.df$w.motif,temp.df$cat),margin=2)
    #prop.table(table(temp.df$w.motif,temp.df$cat),margin=1)
}

bar.plot.df$file.name <- factor(bar.plot.df$file.name,levels=c('JunB','cJun','CEBPbeta','GR','KLF4','KLF5','Tw2_GM'))
bar.plot.df$type <- factor(bar.plot.df$type,levels=c('not.linked','linked'))

pdf(file='ChIP.peaks.w.motif.RE.overlap.bar.pdf',height=5,width=6) 
print(
    ggplot(bar.plot.df,aes(x = file.name,y = perc,fill=factor.name,alpha=type)) +
    geom_bar(stat='Identity',color='black',position='stack') +
    labs(title = 'All ChIP-seq Peaks w/ Motif That Overlap REs', 
         y = 'Percent of Motif-Containing Peaks', 
         x = 'ChIP-seq Factor',
         fill = NULL,
         alpha = NULL) +
    theme_minimal() +
    theme(panel.grid.minor = element_blank(),
          axis.ticks = element_blank(),
          axis.text.x = element_text(angle=45,size=12,hjust=.99,vjust=1,color='black',face='bold'),
          axis.text.y = element_text(size=14,face='bold',color='black'),
          axis.title.x = element_text(size=14,face='bold'),
          axis.title.y = element_text(size=14,face='bold'),
          legend.title = element_text(size=12,face='bold'),
          legend.text = element_text(size=10,face='bold'),
          plot.title = element_text(size=16,face='bold',hjust=0.5)) +
    scale_fill_manual(values = c('#ffcd30','#c469aa','#6aa3ce','#5caa62','#fbae23')) +
    scale_alpha_discrete(range=c(0.5, 1))
)
dev.off()

#Perc. REs Overlapping Peaks

bar.plot.df <- data.frame(file.name = rep(c('CEBPbeta','GR','JunB','KLF4','KLF5','Tw2_GM','cJun'),each = 2),
                      factor.name = rep(c('CEBP','GR','AP1','KLF','KLF','TWIST','AP1'),each=2),
                      type = rep(c('linked','not.linked'),7),
                      perc=0)

for(factor in unique(bar.plot.df$file.name)) {
    name <- factor.df[factor.df$file.name == factor,'factor.name']
    print(name)

    motif.REs <- rownames(final.atac.all.df[!is.na(final.atac.all.df[,name]),])
    total.REs <- length(unique(motif.REs))
    
    temp.df <- plot.df[plot.df$factor == factor & plot.df$RE %in% motif.REs,]
    
    not.linked <- length(unique(temp.df[temp.df$linked == FALSE,'RE']))
    linked <- length(unique(temp.df[temp.df$linked == TRUE,'RE']))

    bar.plot.df[bar.plot.df$file.name == factor & bar.plot.df$type == 'linked','perc'] <- 100 * (linked/total.REs)
    bar.plot.df[bar.plot.df$file.name == factor & bar.plot.df$type == 'not.linked','perc'] <- 100 * (not.linked/total.REs)
    
    #prop.table(table(temp.df$w.motif,temp.df$cat),margin=2)
    #prop.table(table(temp.df$w.motif,temp.df$cat),margin=1)
}

bar.plot.df$file.name <- factor(bar.plot.df$file.name,levels=c('JunB','cJun','CEBPbeta','GR','KLF4','KLF5','Tw2_GM'))
bar.plot.df$type <- factor(bar.plot.df$type,levels=c('not.linked','linked'))

pdf(file='RE.overlap.w.ChIP.peaks.bar.pdf',height=5,width=6) 
print(
    ggplot(bar.plot.df,aes(x = file.name,y = perc,fill=factor.name,alpha=type)) +
    geom_bar(stat='Identity',color='black',position='stack') +
    labs(title = 'All REs with Motif Overlapping ChIP-seq Peaks', 
         y = 'Percent of Motif-Containing REs', 
         x = 'ChIP-seq Factor',
         fill = NULL,
         alpha = NULL) +
    theme_minimal() +
    theme(panel.grid.minor = element_blank(),
          axis.ticks = element_blank(),
          axis.text.x = element_text(angle=45,size=12,hjust=.99,vjust=1,color='black',face='bold'),
          axis.text.y = element_text(size=14,face='bold',color='black'),
          axis.title.x = element_text(size=14,face='bold'),
          axis.title.y = element_text(size=14,face='bold'),
          legend.title = element_text(size=12,face='bold'),
          legend.text = element_text(size=10,face='bold'),
          plot.title = element_text(size=16,face='bold',hjust=0.5)) +
    scale_fill_manual(values = c('#ffcd30','#c469aa','#6aa3ce','#5caa62','#fbae23')) +
    scale_alpha_discrete(range=c(0.5, 1))
)
dev.off()


#Still need:
#compare baseMean of peaks overlapping ChIP-seq peaks v. not
