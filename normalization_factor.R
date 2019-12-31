Args=commandArgs(TRUE)
arg1 = Args[1]
arg2 = Args[2]
arg3 = Args[3]
arg4 = Args[4]
arg5 = Args[5]
arg6 = Args[6]

library(bigWig)

print(arg1)
bigwig1 = load.bigWig(arg1)
bigwig2 = load.bigWig(arg2)
bigwig3 = load.bigWig(arg3)
bigwig4 = load.bigWig(arg4)
bigwig5 = load.bigWig(arg5)
bigwig6 = load.bigWig(arg6)

coverage.bigwig1 = bigwig1$basesCovered*bigwig1$mean
coverage.bigwig2 = bigwig2$basesCovered*bigwig2$mean
coverage.bigwig3 = bigwig3$basesCovered*bigwig3$mean
coverage.bigwig4 = bigwig4$basesCovered*bigwig4$mean
coverage.bigwig5 = bigwig5$basesCovered*bigwig5$mean
coverage.bigwig6 = bigwig6$basesCovered*bigwig6$mean

coverage = c(coverage.bigwig1, coverage.bigwig2, coverage.bigwig3, coverage.bigwig4, coverage.bigwig5, coverage.bigwig6)
file = c(arg1, arg2, arg3, arg4, arg5, arg6)

nm = strsplit(arg1, "_rep1_pro_plus.bigWig")[[1]][1]


write.table(data.frame(file, coverage), file = paste0(nm, '_normalization.txt'), sep = '\t', row.names = FALSE, col.names=FALSE, quote =FALSE)

unload.bigWig(bigwig1)
unload.bigWig(bigwig2)

unload.bigWig(bigwig3)
unload.bigWig(bigwig4)

unload.bigWig(bigwig5)
unload.bigWig(bigwig6)
