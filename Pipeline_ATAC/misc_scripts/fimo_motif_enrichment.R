Args=commandArgs(TRUE)
motif = Args[1]

library(lattice)

setwd('/scratch/abd3x/Adipogenesis/ATAC/fimo_motif_enrichment')

supercluster.key = data.frame(row.names = c('cluster9','cluster23','cluster17','cluster11',
                       'cluster5','cluster8','cluster10','cluster1',
                       'cluster6','cluster2','cluster12',
                       'cluster7','cluster3','cluster4','cluster13','cluster14',
                       'cluster24',
                       'nondynamic'),
                       supercluster = c(rep('up.flat',4),
                                        rep('grad.up',4),
                                        rep('up.down',3),
                                        rep('grad.down',5),
                                        rep('down.up',1),
                                        'na')
                       )

colors = c('#000000','#dddddd')

motif.name = strsplit(motif,'.txt')[[1]][1]
print(motif.name)
table = t(read.table(motif,sep='\t',header=T,row.names=1))

#result = 100*sweep(table, 2, colSums(table), "/")
result = table
result = result[,c('cluster9','cluster23','cluster17','cluster11',
                   'cluster5','cluster8','cluster10','cluster1',
                   'cluster6','cluster2','cluster12',
                   'cluster7','cluster3','cluster4','cluster13','cluster14',
                   'cluster24',
                   'nondynamic')]

sig = FALSE
sig.clusters = c()
for (cluster in colnames(result)) {
    small.table = result[,c(cluster,'nondynamic')]
    output = chisq.test(small.table)
    if (output$p.value < 0.001) {
        
        change = ''
        if ((small.table[1,1] / small.table[2,1]) > (small.table[1,2]/small.table[2,2])) {
            change = 'enriched'
            } else {
                change = 'depleted'
            }

        sc = as.character(supercluster.key[cluster,])      
        
        write(paste0(motif.name,'\t',cluster,'\t',change,'\t',sc),file="significant_motifs.txt",append=TRUE)
        sig=TRUE
    }
    if (sig == TRUE) {            
        pdf(file = paste0(motif.name,'.enrichment.barchart.pdf'),height=9)
        par(las=2)
        barplot(result, col = colors, cex.names= 1.2, legend.text = TRUE, args.legend = list(x=2.5,y=113), main = paste0(motif.name,' Enrichment'))
        dev.off()
    }
}
