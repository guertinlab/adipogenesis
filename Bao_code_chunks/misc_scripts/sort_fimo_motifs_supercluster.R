#Organize all the fimo identified motifs by supercluster and save them in new directory

setwd('/scratch/bhn9by/ATAC')
fimo = read.table('fimo_motif_enrichment/significant_motifs.txt',sep='\t',header=F)

colnames(fimo) = c('motif','cluster','change','supercluster')
enriched = fimo[fimo$change == 'enriched',]

#had to edit the jaspar motif names some during FIMO, changing them back here
enriched$motif = gsub('.var.2.','(var.2)',enriched$motif)
enriched$motif = gsub('\\.\\.','::',enriched$motif)

#enriched$new.motif = paste0(enriched$motif,'_',enriched$database)

#up.flat - 9,23,17,11
#grad.up - 5,8,10,1
#up.down - 6,2,12
#grad.down - 7,3,4,13,14
#down.up - 24

superclusters = unique(enriched$supercluster)

for (sc in superclusters) {
    print(sc)
    write((as.vector(unique(enriched[enriched$supercluster == sc,]$motif))),file=paste0('supercluster_fimo_motifs/',sc,'.fimo.motifs.txt'))
}
