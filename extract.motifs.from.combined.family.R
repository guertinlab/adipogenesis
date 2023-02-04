library(ggseqlogo)
library(lattice)

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

setwd('/scratch/abd3x/ATAC/SP_KLF_split/')

load('/scratch/abd3x/ATAC/plot.df.atac.Rdata')
load('/scratch/abd3x/ATAC/fimo.scores.atac.Rdata')

plot.df = plot.df.atac
plot.df$status = 'Activated'
plot.df[plot.df$supercluster == 'gradual.down' | plot.df$supercluster == 'down.up',]$status = 'Repressed'

bed.file = '/scratch/abd3x/ATAC/fimo_composites/PSWM_family_7_fimo.bed'
factor.name = 'SP.KLF'
x = read.table(bed.file, stringsAsFactors=FALSE)
x = x[x[,6] != -1,]
x = x[,c(1,2,3,10)]
x$factor = factor.name
colnames(x) = c('chr', 'start', 'end', 'sequence', 'factor')
x$re = paste0(x[,1], ':', x[,2], '-', x[,3])
df.seq = x

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

#plot dynamic traces of all SP/KLF peaks
peaks = df.seq$re

#for(i in 1:4) {
#    temp.peaks = rownames(fimo.scores.atac[!is.na(fimo.scores.atac[,i]),])
#    peaks = setdiff(peaks,temp.peaks)
#}

plot.df = plot.df.atac[plot.df.atac$genes %in% peaks,]
plot.df$family = 'SP.KLF'

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

pdf(file = 'SP_KLF_dynamic_accessibility.pdf')

trellis.par.set(box.umbrella = list(lty = 1, col="black", lwd=1),
                box.rectangle = list( lwd=1.0, col="black", alpha = 1.0),
                plot.symbol = list(col="black", lwd=1.0, pch ='.'))

print(
    xyplot(value ~  sample.conditions, group = genes, data = plot.df, type = c('l'),#type = c('l','p'),
       scales=list(x=list(cex=1.5,relation = "free", rot = 45), y =list(cex=1.5, relation="free")),
       aspect=1.0,
       between=list(y=0.5, x=0.5),
       main = list(label = 'SP/KLF Motifs in Dynamic Peaks', cex = 2.0),
       ylab = list(label = 'Normalized ATAC signal', cex = 2.0),
       xlab = list(label = 'Time (minutes)', cex = 2.0),
       par.settings = list(superpose.symbol = list(pch = c(16),col=c('grey20'), cex =0.5),
                                   #strip.background=list(col="grey80"),
                                   superpose.line = list(col = as.character(cat.colours$colour), lwd=c(1),lty = c(1))),
       panel = function(x, y, ...) {
           panel.xyplot(x, y, ...)
           #panel.bwplot(x, y, pch = '|', horizontal = FALSE, box.width = 15, do.out = FALSE)
           #panel.spline(x, y, col ='blue', lwd =2.0, ...) 
           
       })

)
dev.off()
