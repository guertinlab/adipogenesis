setwd('/scratch/bhn9by/ATAC/SP_KLF_split')

x = read.table('/scratch/abd3x/ATAC/fimo_composites/PSWM_family_7_fimo.bed', stringsAsFactors=FALSE)
x = x[x[,6] != -1,]

a <- apply(x[,c(1:9) ], 1 ,paste0, collapse = ":" )
a = gsub(" ", "", a, fixed = TRUE)
b = x[,10]
for (i in 1:length(a)) {
    write.table(paste('>', a[i], sep =''), file = 'sp_fimo.txt', append=TRUE, col.names=FALSE, row.names=FALSE, sep ='', quote= FALSE)
    write.table(b[i], file = 'sp_fimo.txt', append=TRUE, col.names=FALSE, row.names=FALSE, sep ='', quote= FALSE)
}

x = read.table('/scratch/abd3x/ATAC/fimo_composites/PSWM_family_7_fimo_nondyn.bed', stringsAsFactors=FALSE)
x = x[x[,6] != -1,]

a <- apply(x[,c(1:9) ], 1 ,paste0, collapse = ":" )
a = gsub(" ", "", a, fixed = TRUE)
b = x[,10]
for (i in 1:length(a)) {
    write.table(paste('>', a[i], sep =''), file = 'sp_fimo_nondyn.txt', append=TRUE, col.names=FALSE, row.names=FALSE, sep ='', quote= FALSE)
    write.table(b[i], file = 'sp_fimo_nondyn.txt', append=TRUE, col.names=FALSE, row.names=FALSE, sep ='', quote= FALSE)
}

#takes a long time, should figure out a better way
x = read.table('/scratch/abd3x/ATAC/fimo_composites/PSWM_family_7_2M.txt', stringsAsFactors=FALSE)
a <- apply(x[,c(1:6) ], 1 ,paste0, collapse = ":" )
a = gsub(" ", "", a, fixed = TRUE)
b = x[,7]
for (i in 1:length(a)) {
    write.table(paste('>', a[i], sep =''), file = 'sp_klf_2M.txt', append=TRUE, col.names=FALSE, row.names=FALSE, sep ='', quote= FALSE)
    write.table(b[i], file = 'sp_klf_2M.txt', append=TRUE, col.names=FALSE, row.names=FALSE, sep ='', quote= FALSE)
}
