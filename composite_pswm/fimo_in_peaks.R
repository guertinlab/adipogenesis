                                        #how I got: dynamic_peaks.bed

chrom = vector(mode="character", length = length(unique(plot.df.atac$genes)))
start = vector(mode="integer", length = length(unique(plot.df.atac$genes)))
end = vector(mode="integer", length = length(unique(plot.df.atac$genes)))

count = 0
for (j in unique(plot.df.atac$genes)) {
    count = count +1
    chrom[count] = strsplit(j, ":")[[1]][1]
    start[count] = as.numeric(strsplit(strsplit(j, ":")[[1]][2], "-")[[1]][1])
    end[count] = as.numeric(strsplit(strsplit(j, ":")[[1]][2], "-")[[1]][2])
}













#get fimo info
res = data.frame(cbind(chrom, start, end))
colnames(res) = c('chr', 'start', 'end')
res$start = as.numeric(as.character(res$start))
res$end = as.numeric(as.character(res$end))

write.table(res, file='dynamic_peaks.bed', sep='\t', quote=F, row.names=F, col.names=F)
rownames(res) = paste0(res[,1], ':', res[,2], '-', res[,3])


res = head(res, 40)

count = -1
vec.names = c()
for(bed.file in Sys.glob(file.path('/Users/guertinlab/Desktop/fimo_stuff/*family_fimo.bed'))) {
    print(bed.file)
    count = count + 1
    factor.name = strsplit(bed.file, "/")[[1]]
    factor.name = strsplit(factor.name[length(factor.name)],
                           '_fimo.bed')[[1]][1]
    print(factor.name)
    x = read.table(bed.file, stringsAsFactors=FALSE)
    x = x[x[,6] != -1,]
    y = aggregate(as.numeric(V8)~V1+V2+V3, data=x, FUN=sum)
    colnames(y) = c('chr', 'start', 'end', factor.name)
    rownames(y) = paste0(y[,1], ':', y[,2], '-', y[,3])
    print(head(y))
    res = merge(res, y, by="row.names", all = TRUE)
    rnames = res[,1]
    print(head(res))
    vec.names = c(vec.names, factor.name)
    res = data.frame(res[,-c((ncol(res) -3 ):(ncol(res) -1))])
    res = data.frame(res[,(ncol(res) - count):ncol(res)])
    colnames(res) = vec.names
    rownames(res) = rnames
    print(head(res))
}

fimo.scores.atac = res
save(fimo.scores.atac, file = 'fimo.scores.atac.Rdata')


sapply(fimo.scores.atac, function(x) sum(!is.na(x)))

barplot(sapply(fimo.scores.atac, function(x) sum(!is.na(x))))
