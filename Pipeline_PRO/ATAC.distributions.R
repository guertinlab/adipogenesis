library(ggplot2)
library(wesanderson)

setwd('~/Adipogenesis/ATAC_analysis_redo')

#dynamic

x = read.table('dynamic_ATAC_peaks_promoters_temp.bed')
x$location = paste0(x$V1,':',x$V2,'-',x$V3)

y = read.table('dynamic_ATAC_peaks_intragenic_temp.bed')
y$location = paste0(y$V1,':',y$V2,'-',y$V3)

z = read.table('dynamic_peaks.bed')
z$location = paste0(z$V1,':',z$V2,'-',z$V3)

promoters = x$location
intragenic = setdiff(y$location,x$location)
intergenic = setdiff(z$location,y$location)
intergenic = setdiff(intergenic,x$location)

#write promoter dynamic ATAC peaks to bed file
chr = sapply(strsplit(promoters, ':'), '[', 1)
x = sapply(strsplit(promoters, ':'), '[', 2)
start = sapply(strsplit(x, '-'), '[', 1)
end = sapply(strsplit(x, '-'), '[', 2)
bed = data.frame(chr,start,end)
write.table(bed,file='dynamic_ATAC_peaks_promoters.bed',sep='\t',quote=F,col.names = F,row.names = F)

#write intragenic dynamic ATAC peaks to bed file
chr = sapply(strsplit(intragenic, ':'), '[', 1)
x = sapply(strsplit(intragenic, ':'), '[', 2)
start = sapply(strsplit(x, '-'), '[', 1)
end = sapply(strsplit(x, '-'), '[', 2)
bed = data.frame(chr,start,end)
write.table(bed,file='dynamic_ATAC_peaks_intragenic.bed',sep='\t',quote=F,col.names = F,row.names = F)

#write intergenic dynamic ATAC peaks to bed file
chr = sapply(strsplit(intergenic, ':'), '[', 1)
x = sapply(strsplit(intergenic, ':'), '[', 2)
start = sapply(strsplit(x, '-'), '[', 1)
end = sapply(strsplit(x, '-'), '[', 2)
bed = data.frame(chr,start,end)
write.table(bed,file='dynamic_ATAC_peaks_intergenic.bed',sep='\t',quote=F,col.names = F,row.names = F)

system('intersectBed -wa -a dynamic_ATAC_peaks_promoters.bed -b ~/Adipogenesis/PRO_analysis_redo/dREG/dynamic_dREG_peaks.bed -u > dynamic_ATAC_atac_overlap_promoters.bed')
system('intersectBed -wa -a dynamic_ATAC_peaks_intragenic.bed -b ~/Adipogenesis/PRO_analysis_redo/dREG/dynamic_dREG_peaks.bed -u > dynamic_ATAC_atac_overlap_intragenic.bed')
system('intersectBed -wa -a dynamic_ATAC_peaks_intergenic.bed -b ~/Adipogenesis/PRO_analysis_redo/dREG/dynamic_dREG_peaks.bed -u > dynamic_ATAC_atac_overlap_intergenic.bed')

x = read.table('dynamic_ATAC_atac_overlap_promoters.bed')
x$location = paste0(x$V1,':',x$V2,'-',x$V3)

y = read.table('dynamic_ATAC_atac_overlap_intragenic.bed')
y$location = paste0(y$V1,':',y$V2,'-',y$V3)

z = read.table('dynamic_ATAC_atac_overlap_intergenic.bed')
z$location = paste0(z$V1,':',z$V2,'-',z$V3)

distributions = data.frame(
  type = rep(c('Promoters','Intragenic','Intergenic'),2),
  total = c(length(promoters)-nrow(x),length(intragenic)-nrow(y),length(intergenic)-nrow(z),nrow(x),nrow(y),nrow(z)),
  atac.overlap = rep(c('not.overlap','overlap'),each = 3)
  )
distributions$prop = distributions$total / sum(distributions$total)

x = read.table('dynamic.expected.distributions.txt',header=T)

distributions$z.score = NA
distributions[distributions$type == 'Promoters',]$z.score = round((sum(distributions[distributions$type == 'Promoters',]$total) - mean(x$promoters))/sd(x$promoters),2)
distributions[distributions$type == 'Intragenic',]$z.score = round((sum(distributions[distributions$type == 'Intragenic',]$total) - mean(x$intragenic))/sd(x$intragenic),2)
distributions[distributions$type == 'Intergenic',]$z.score = round((sum(distributions[distributions$type == 'Intergenic',]$total) - mean(x$intergenic))/sd(x$intergenic),2)

distributions$pos = NA
distributions[distributions$type == 'Promoters',]$pos = sum(distributions[distributions$type == 'Promoters',]$prop) + 0.02
distributions[distributions$type == 'Intragenic',]$pos = sum(distributions[distributions$type == 'Intragenic',]$prop) + 0.02
distributions[distributions$type == 'Intergenic',]$pos = sum(distributions[distributions$type == 'Intergenic',]$prop) + 0.02

pdf(file='prop.dyn.ATAC.peaks.distributions.pdf',width=6,height=6)
print(
  ggplot(data = distributions, aes(x = type,y = prop,fill=atac.overlap)) +
    geom_bar(stat='identity',position='stack',color='black') +
    geom_text(aes(y = pos,label=z.score),size=5) +
    labs(y = '% of All Dynamic ATAC Peaks', x = NULL, fill = 'Dynamic dREG Overlap') +
    theme_minimal() +
    theme(panel.grid.minor = element_blank(),
          axis.ticks = element_blank(),
          axis.text.x = element_text(size=12,face='bold',color='black',angle=45),
          axis.text.y = element_text(size=12,face='bold',color='black'),
          axis.title.y = element_text(size=14,face='bold'),
          legend.title = element_text(size=12,face='bold'),
          legend.text = element_text(size=12,face='bold')) +
    scale_fill_manual(values = wes_palette("GrandBudapest2")[c(3,4)])
)
dev.off()

for (file in Sys.glob(file.path('~/Adipogenesis/ATAC_analysis_redo/fimo_composites/*instances.bed'))) {
  name = strsplit(strsplit(file,'/')[[1]][7],'_instances.bed')[[1]][1]
  print(name)
  x = read.table(file)
  x$location = paste0(x[,1],':',x[,2],'-',x[,3])
  promoter.overlap = intersect(x$location,promoters)
  intragenic.overlap = intersect(x$location,intragenic)
  intergenic.overlap = intersect(x$location,intergenic)
  plot.df = data.frame(type = c('promoters','intragenic','intergenic'),
                       num = 100*c(length(promoter.overlap)/nrow(x),
                               length(intragenic.overlap)/nrow(x),
                               length(intergenic.overlap)/nrow(x)))
  
  pdf(file=paste0('prop.',name,'.ATAC.peaks.distributions.pdf'),width=6,height=6)
  print(
    ggplot(data = plot.df, aes(x = type,y = num)) +
      geom_bar(stat='identity',color='black',fill='black') +
      #geom_text(aes(y = pos,label=z.score),size=5) +
      labs(y = paste0('% of All ',name,' Dynamic ATAC Peaks'), x = NULL) +
      theme_minimal() +
      theme(panel.grid.minor = element_blank(),
            axis.ticks = element_blank(),
            axis.text.x = element_text(size=12,face='bold',color='black',angle=45),
            axis.text.y = element_text(size=12,face='bold',color='black'),
            axis.title.y = element_text(size=14,face='bold'),
            legend.title = element_text(size=12,face='bold'),
            legend.text = element_text(size=12,face='bold'))
  )
  dev.off()
}



#all

x = read.table('all_ATAC_peaks_promoters_temp.bed')
x$location = paste0(x$V1,':',x$V2,'-',x$V3)

y = read.table('all_ATAC_peaks_intragenic_temp.bed')
y$location = paste0(y$V1,':',y$V2,'-',y$V3)

z = read.table('all_peaks.bed')
z$location = paste0(z$V1,':',z$V2,'-',z$V3)

promoters = x$location
intragenic = setdiff(y$location,x$location)
intergenic = setdiff(z$location,y$location)
intergenic = setdiff(intergenic,x$location)

#write promoter all ATAC peaks to bed file
chr = sapply(strsplit(promoters, ':'), '[', 1)
x = sapply(strsplit(promoters, ':'), '[', 2)
start = sapply(strsplit(x, '-'), '[', 1)
end = sapply(strsplit(x, '-'), '[', 2)
bed = data.frame(chr,start,end)
write.table(bed,file='all_ATAC_peaks_promoters.bed',sep='\t',quote=F,col.names = F,row.names = F)

#write intragenic all ATAC peaks to bed file
chr = sapply(strsplit(intragenic, ':'), '[', 1)
x = sapply(strsplit(intragenic, ':'), '[', 2)
start = sapply(strsplit(x, '-'), '[', 1)
end = sapply(strsplit(x, '-'), '[', 2)
bed = data.frame(chr,start,end)
write.table(bed,file='all_ATAC_peaks_intragenic.bed',sep='\t',quote=F,col.names = F,row.names = F)

#write intergenic all ATAC peaks to bed file
chr = sapply(strsplit(intergenic, ':'), '[', 1)
x = sapply(strsplit(intergenic, ':'), '[', 2)
start = sapply(strsplit(x, '-'), '[', 1)
end = sapply(strsplit(x, '-'), '[', 2)
bed = data.frame(chr,start,end)
write.table(bed,file='all_ATAC_peaks_intergenic.bed',sep='\t',quote=F,col.names = F,row.names = F)

system('intersectBed -wa -a all_ATAC_peaks_promoters.bed -b ~/Adipogenesis/PRO_analysis_redo/dREG/all_dREG_peaks.bed -u > all_ATAC_atac_overlap_promoters.bed')
system('intersectBed -wa -a all_ATAC_peaks_intragenic.bed -b ~/Adipogenesis/PRO_analysis_redo/dREG/all_dREG_peaks.bed -u > all_ATAC_atac_overlap_intragenic.bed')
system('intersectBed -wa -a all_ATAC_peaks_intergenic.bed -b ~/Adipogenesis/PRO_analysis_redo/dREG/all_dREG_peaks.bed -u > all_ATAC_atac_overlap_intergenic.bed')

x = read.table('all_ATAC_atac_overlap_promoters.bed')
x$location = paste0(x$V1,':',x$V2,'-',x$V3)

y = read.table('all_ATAC_atac_overlap_intragenic.bed')
y$location = paste0(y$V1,':',y$V2,'-',y$V3)

z = read.table('all_ATAC_atac_overlap_intergenic.bed')
z$location = paste0(z$V1,':',z$V2,'-',z$V3)

distributions = data.frame(
  type = rep(c('Promoters','Intragenic','Intergenic'),2),
  total = c(length(promoters)-nrow(x),length(intragenic)-nrow(y),length(intergenic)-nrow(z),nrow(x),nrow(y),nrow(z)),
  atac.overlap = rep(c('not.overlap','overlap'),each = 3)
)
distributions$prop = distributions$total / sum(distributions$total)

x = read.table('all_expected.distributions.txt',header=T)

distributions$z.score = NA
distributions[distributions$type == 'Promoters',]$z.score = round((sum(distributions[distributions$type == 'Promoters',]$total) - mean(x$promoters))/sd(x$promoters),2)
distributions[distributions$type == 'Intragenic',]$z.score = round((sum(distributions[distributions$type == 'Intragenic',]$total) - mean(x$intragenic))/sd(x$intragenic),2)
distributions[distributions$type == 'Intergenic',]$z.score = round((sum(distributions[distributions$type == 'Intergenic',]$total) - mean(x$intergenic))/sd(x$intergenic),2)

distributions$pos = NA
distributions[distributions$type == 'Promoters',]$pos = sum(distributions[distributions$type == 'Promoters',]$prop) + 0.02
distributions[distributions$type == 'Intragenic',]$pos = sum(distributions[distributions$type == 'Intragenic',]$prop) + 0.02
distributions[distributions$type == 'Intergenic',]$pos = sum(distributions[distributions$type == 'Intergenic',]$prop) + 0.02

pdf(file='prop.all.ATAC.peaks.distributions.pdf',width=6,height=6)
print(
  ggplot(data = distributions, aes(x = type,y = prop,fill=atac.overlap)) +
    geom_bar(stat='identity',position='stack',color='black') +
    geom_text(aes(y = pos,label=z.score),size=5) +
    labs(y = '% of All ATAC Peaks', x = NULL, fill = 'All dREG Overlap') +
    theme_minimal() +
    theme(panel.grid.minor = element_blank(),
          axis.ticks = element_blank(),
          axis.text.x = element_text(size=12,face='bold',color='black',angle=45),
          axis.text.y = element_text(size=12,face='bold',color='black'),
          axis.title.y = element_text(size=14,face='bold'),
          legend.title = element_text(size=12,face='bold'),
          legend.text = element_text(size=12,face='bold')) +
    scale_fill_manual(values = wes_palette("GrandBudapest2")[c(3,4)])
)
dev.off()

