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


                                        #from Arun: load('/Users/guertinlab/Desktop/fimo_stuff/plot.df.atac.Rdata')

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


#add this to supercluster plot:
#this takes a while to run
for (i in 1:nrow(plot.df.atac)) {
    re = plot.df.atac$genes[i] 
    plot.df.atac[i, 11:15] = fimo.scores.atac[re,]
}

#save(plot.df.atac, file = 'plot.df.atac.MJG.Rdata')

#organizing the data to plot and get the traces correct.
plot.df.atac$supercluster <- as.factor(plot.df.atac$supercluster)
plot.df.atac$genes <- as.factor(plot.df.atac$genes)

cat.colours = plot.df.atac[plot.df.atac$merge == 'one_groupt0', c(1,10)]
cat.colours$genes <- as.factor(cat.colours$genes)
cat.colours$supercluster <- as.factor(cat.colours$supercluster)


cat.colours$colour[cat.colours$supercluster == 'up.flat'] <- '#FF000008'
cat.colours$colour[cat.colours$supercluster == 'up.down'] <- '#FF000008'
cat.colours$colour[cat.colours$supercluster == 'gradual.up'] <- '#FF000008'
cat.colours$colour[cat.colours$supercluster == 'gradual.down'] <- '#0000FF08'
cat.colours$colour[cat.colours$supercluster == 'down.up'] <- '#0000FF08'

cat.colours$colour <- as.factor(cat.colours$colour)

cat.colours <- cat.colours[match(levels(plot.df.atac$genes), cat.colours$genes), ]

plot.df.atac$sc[plot.df.atac$supercluster == 'up.flat'] <- 'Immediate Increase'
plot.df.atac$sc[plot.df.atac$supercluster == 'up.down'] <- 'Transient Increase'
plot.df.atac$sc[plot.df.atac$supercluster == 'gradual.up'] <- 'Gradual Increase'
plot.df.atac$sc[plot.df.atac$supercluster == 'gradual.down'] <- 'Gradual Decrease'
plot.df.atac$sc[plot.df.atac$supercluster == 'down.up'] <- 'Transient Decrease'


pdf(file = paste0('supercluster.traces.transcol.pdf'),width=12,height=4)

trellis.par.set(box.umbrella = list(lty = 1, col="black", lwd=1),
                box.rectangle = list( lwd=1.0, col="black", alpha = 1.0),
                plot.symbol = list(col="black", lwd=1.0, pch ='.'))

print(
    xyplot(value ~  sample.conditions | sc, group = genes, data = plot.df.atac, type = c('l'),#type = c('l','p'),
           scales=list(x=list(cex=1.0,relation = "free", rot = 45), y = list(cex=1.0, relation="free")),
           aspect=1.0,
           layout = c(5,1),
           between=list(y=0.5, x=0.5),
           index.cond=list(c(3,5,2,4,1)),
           ylab = list(label = 'Normalized ATAC signal', cex =1.0),
           xlab = list(label = 'Time (minutes)', cex =1.0),
           par.settings = list(superpose.symbol = list(pch = c(16),col=c('grey20'), cex =0.5),
                               strip.background=list(col="grey80"),
                               superpose.line = list(col = as.character(cat.colours$colour), lwd=c(1),lty = c(1))),
           panel = function(x, y, ...) {
               panel.xyplot(x, y, ...)
               #panel.violin(x, y, horizontal = FALSE, col = "transparent",
               #         box.width = 20, do.out = FALSE)
               panel.bwplot(x, y, pch = '|', horizontal = FALSE, box.width = 10, do.out = FALSE)
               panel.spline(x, y, col = 'grey70', lwd =3.0, ...) 
           })
)

dev.off()

sum(plot.df.atac$sc == 'Immediate Increase')/7
sum(plot.df.atac$sc == 'Transient Increase')/7
sum(plot.df.atac$sc == 'Gradual Increase')/7
sum(plot.df.atac$sc == 'Gradual Decrease')/7
sum(plot.df.atac$sc == 'Transient Decrease')/7


plot.df.atac$AP1 <- ifelse(is.na(plot.df.atac$AP1_family), FALSE, TRUE)
plot.df.atac$GR <- ifelse(is.na(plot.df.atac$GR_family), FALSE, TRUE)
plot.df.atac$SP <- ifelse(is.na(plot.df.atac$SP_family), FALSE, TRUE)
plot.df.atac$TWIST <- ifelse(is.na(plot.df.atac$TWIST_family), FALSE, TRUE)
plot.df.atac$bHLH <- ifelse(is.na(plot.df.atac$bHLH_family), FALSE, TRUE)


plot.df.atac$status[plot.df.atac$supercluster == 'up.flat'] <- 'Activated'
plot.df.atac$status[plot.df.atac$supercluster == 'up.down'] <- 'Activated'
plot.df.atac$status[plot.df.atac$supercluster == 'gradual.up'] <- 'Activated'
plot.df.atac$status[plot.df.atac$supercluster == 'gradual.down'] <- 'Repressed'
plot.df.atac$status[plot.df.atac$supercluster == 'down.up'] <- 'Repressed'

save(plot.df.atac, file = 'plot.df.atac.MJG.Rdata')

#plot.df.atac now has all the column (I believe)


ap1 = plot.df.atac[plot.df.atac$AP1 == TRUE & plot.df.atac$GR == FALSE & plot.df.atac$SP == FALSE & plot.df.atac$TWIST == FALSE & plot.df.atac$bHLH == FALSE,]
ap1$motif = 'AP1 family'
gr = plot.df.atac[plot.df.atac$AP1 == FALSE & plot.df.atac$GR == TRUE & plot.df.atac$SP == FALSE & plot.df.atac$TWIST == FALSE & plot.df.atac$bHLH == FALSE,]
gr$motif = 'GR family'
sp = plot.df.atac[plot.df.atac$AP1 == FALSE & plot.df.atac$GR == FALSE & plot.df.atac$SP == TRUE & plot.df.atac$TWIST == FALSE & plot.df.atac$bHLH == FALSE,]
sp$motif = 'SP family'
twist = plot.df.atac[plot.df.atac$AP1 == FALSE & plot.df.atac$GR == FALSE & plot.df.atac$SP == FALSE & plot.df.atac$TWIST == TRUE & plot.df.atac$bHLH == FALSE,]
twist$motif = 'TWIST family'
bhlh = plot.df.atac[plot.df.atac$AP1 == FALSE & plot.df.atac$GR == FALSE & plot.df.atac$SP == FALSE & plot.df.atac$TWIST == FALSE & plot.df.atac$bHLH == TRUE,]
bhlh$motif = 'bHLH family'


fig.2j = rbind(ap1, gr, sp, twist, bhlh)

#added these value in AI
nrow(ap1[ap1$status == 'Activated',])/7
nrow(ap1[ap1$status == 'Repressed',])/7

nrow(gr[gr$status == 'Activated',])/7
nrow(gr[gr$status == 'Repressed',])/7

nrow(sp[sp$status == 'Activated',])/7
nrow(sp[sp$status == 'Repressed',])/7

nrow(twist[twist$status == 'Activated',])/7
nrow(twist[twist$status == 'Repressed',])/7

nrow(bhlh[bhlh$status == 'Activated',])/7
nrow(bhlh[bhlh$status == 'Repressed',])/7


pdf(file = paste0('2j.traces.pdf'),width=10,height=4)
trellis.par.set(box.umbrella = list(lty = 1, col="black", lwd=1),
                box.rectangle = list( lwd=1.0, col="black", alpha = 1.0),
                plot.symbol = list(col="black", lwd=1.0, pch ='.'))
print(
    xyplot(value ~  sample.conditions | motif, group = genes, data = fig.2j, type = c('l'),#type = c('l','p'),
           scales=list(x=list(cex=1.0,relation = "free", rot = 45), y = list(cex=1.0, relation="free")),
           aspect=1.0,
           between=list(y=0.5, x=0.5),
           ylab = list(label = 'Normalized ATAC signal', cex =1.0),
           index.cond=list(c(1,3,4,5,2)),
           xlab = list(label = 'Time (minutes)', cex =1.0),
           par.settings = list(superpose.symbol = list(pch = c(16),col=c('grey20'), cex =0.5),
                               #I want to change the background strip to the corresponding motif color
                               strip.background=list(col=c("grey80", "green", "purple", "yellow","teal")),
                               superpose.line = list(col = as.character(cat.colours$colour), lwd=c(1),lty = c(1))),
)
)
dev.off()


ap1.gr = plot.df.atac[(plot.df.atac$AP1 == TRUE | plot.df.atac$GR == TRUE) & plot.df.atac$SP == FALSE & plot.df.atac$TWIST == FALSE & plot.df.atac$bHLH == FALSE,]
ap1.gr$motif = 'AP1 or GR only'

sp.twist.bhlh = plot.df.atac[plot.df.atac$AP1 == FALSE & plot.df.atac$GR == FALSE & (plot.df.atac$SP == TRUE | plot.df.atac$TWIST == TRUE | plot.df.atac$bHLH == TRUE),]
sp.twist.bhlh$motif = 'SP or TWIST or bHLH only'

#twist.bhlh = plot.df.atac[plot.df.atac$AP1 == FALSE & plot.df.atac$GR == FALSE & plot.df.atac$SP == FALSE & (plot.df.atac$TWIST == TRUE | plot.df.atac$bHLH == TRUE),]
#twist.bhlh$motif = 'TWIST and/or bHLH '

nrow(ap1.gr[ap1.gr$status == 'Activated',])/7
nrow(ap1.gr[ap1.gr$status == 'Repressed',])/7

nrow(sp.twist.bhlh[sp.twist.bhlh$status == 'Activated',])/7
nrow(sp.twist.bhlh[sp.twist.bhlh$status == 'Repressed',])/7

#nrow(twist.bhlh[twist.bhlh$status == 'Activated',])/7
#nrow(twist.bhlh[twist.bhlh$status == 'Repressed',])/7


fig.2j.2 = rbind(ap1.gr, sp.twist.bhlh)


pdf(file = paste0('2j.traces.act.rep.pdf'),width=10,height=4.0)
trellis.par.set(box.umbrella = list(lty = 1, col="black", lwd=1),
                box.rectangle = list( lwd=1.0, col="black", alpha = 1.0),
                plot.symbol = list(col="black", lwd=1.0, pch ='.'))
print(
    xyplot(value ~  sample.conditions | motif, group = genes, data = fig.2j.2, type = c('l'),#type = c('l','p'),
           scales=list(x=list(cex=1.0,relation = "free", rot = 45), y = list(cex=1.0, relation="free")),
           aspect=1.0,
           between=list(y=0.5, x=0.5),
           ylab = list(label = 'Normalized ATAC signal', cex =1.0),
           xlab = list(label = 'Time (minutes)', cex =1.0),
           par.settings = list(superpose.symbol = list(pch = c(16),col=c('grey20'), cex =0.5),
                               #I want to change the background strip to the corresponding motif color
                               strip.background=list(col=c("grey80", "green", "purple", "yellow","teal")),
                               superpose.line = list(col = as.character(cat.colours$colour), lwd=c(1),lty = c(1))),
)
)
dev.off()

                                        #are GREs and AP1 motifs in the activated?
act.triple = unique(sp.twist.bhlh[sp.twist.bhlh$status == 'Activated',]$genes)

chrom = vector(mode="character", length = length(act.triple))
start = vector(mode="integer", length = length(act.triple))
end = vector(mode="integer", length = length(act.triple))

count = 0
for (j in act.triple) {
    count = count +1
    chrom[count] = strsplit(j, ":")[[1]][1]
    start[count] = as.numeric(strsplit(strsplit(j, ":")[[1]][2], "-")[[1]][1])
    end[count] = as.numeric(strsplit(strsplit(j, ":")[[1]][2], "-")[[1]][2])
}

res = data.frame(cbind(chrom, start, end))
colnames(res) = c('chr', 'start', 'end')
res$start = as.numeric(as.character(res$start))
res$end = as.numeric(as.character(res$end))

write.table(res, file='activated_triple.bed', sep='\t', quote=F, row.names=F, col.names=F)

#fastaFromBed -fi /Volumes/GUERTIN_2/adipogenesis/atac/mm10.fa -bed activated_triple.bed -fo activated_triple.fasta

#meme -oc activated_triple.meme_output -nmotifs 5 -objfun classic -evt 0.01 -searchsize 0 -minw 6 -maxw 15 -revcomp -dna -markov_order 3 -maxsize 100000000 activated_triple.fasta



#makes the box whiskey plot of motif scores.
fig.2j.uniq = fig.2j[!duplicated(paste0(fig.2j$genes, fig.2j$motif)),]

x = rbind(tail(head(fig.2j.uniq, 1000)), tail(head(fig.2j.uniq, 100)), tail(head(fig.2j.uniq, 7000)), tail(head(fig.2j.uniq, 6500)))
                                        #

fig.2j.uniq$new[fig.2j.uniq$AP1 == TRUE] <- fig.2j.uniq$AP1_family[fig.2j.uniq$AP1 == TRUE]
fig.2j.uniq$new[fig.2j.uniq$GR == TRUE] <- fig.2j.uniq$GR_family[fig.2j.uniq$GR == TRUE]
fig.2j.uniq$new[fig.2j.uniq$SP == TRUE] <- fig.2j.uniq$SP_family[fig.2j.uniq$SP == TRUE]
fig.2j.uniq$new[fig.2j.uniq$TWIST == TRUE] <- fig.2j.uniq$TWIST_family[fig.2j.uniq$TWIST == TRUE]
fig.2j.uniq$new[fig.2j.uniq$bHLH == TRUE] <- fig.2j.uniq$bHLH_family[fig.2j.uniq$bHLH == TRUE]

fig.2j.uniq$new <- as.numeric(as.character(fig.2j.uniq$new))



pdf(file = paste0('2j.bw.pdf'),width=8,height=4)
trellis.par.set(box.umbrella = list(lty = 1, col=c("red", "blue"), lwd=2),
                box.rectangle = list( lwd=2.0, col=c("red", "blue"), alpha = 1.0),
                plot.symbol = list(col=c("red", "blue"), lwd=2.0, pch ='.'))
print(
    
bwplot(log(as.numeric(as.character(new)), base =10) ~ status | motif, data = fig.2j.uniq, horizontal = FALSE, pch = '|', do.out = FALSE,
       scales=list(x=list(cex=1.0, relation = "free", rot = 45), y = list(cex=1.0, relation="free")),
       aspect=2.0,
       between=list(y=0.5, x=0.5),
       index.cond=list(c(1,3,4,5,2)),
       ylab = expression('log'[10]* paste(Sigma,'(motif scores)')),
       xlab = expression('ATAC Peak category'),
#manually settign to avoid outlier, since do.out = FALSE       
       ylim = list(c(0.9, 1.9),  c(0.9, 1.7), c(0.8, 1.8), c(1, 2.15),  c(0.96, 1.4)),
       par.settings = list(superpose.symbol = list(pch = c(16),col=c('grey20'), cex =0.5),
                               #I want to change the background strip to the corresponding motif color
                               strip.background=list(col=c("grey80"))))

)
dev.off()


                                        #need to plot seqLogo


count = -1
vec.names = c()

df.seq = as.data.frame(matrix(ncol=5, nrow=0),stringsAsFactors = FALSE)
colnames(df.seq) = c('chr', 'start', 'end', 'sequence', 'factor')

for(bed.file in Sys.glob(file.path('/Users/guertinlab/Desktop/fimo_stuff/*family_fimo.bed'))) {
    print(bed.file)
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

                                        #df.seq has the sequences

df.seq.ap1.rep = df.seq[df.seq$re %in% fig.2j.uniq[fig.2j.uniq$status == 'Repressed' & fig.2j.uniq$motif == 'AP1 family',]$genes,]
df.seq.ap1.act = df.seq[df.seq$re %in% fig.2j.uniq[fig.2j.uniq$status == 'Activated' & fig.2j.uniq$motif == 'AP1 family',]$genes,]

df.seq.gr.rep = df.seq[df.seq$re %in% fig.2j.uniq[fig.2j.uniq$status == 'Repressed' & fig.2j.uniq$motif == 'GR family',]$genes,]
df.seq.gr.act = df.seq[df.seq$re %in% fig.2j.uniq[fig.2j.uniq$status == 'Activated' & fig.2j.uniq$motif == 'GR family',]$genes,]

df.seq.sp.rep = df.seq[df.seq$re %in% fig.2j.uniq[fig.2j.uniq$status == 'Repressed' & fig.2j.uniq$motif == 'SP family',]$genes,]
df.seq.sp.act = df.seq[df.seq$re %in% fig.2j.uniq[fig.2j.uniq$status == 'Activated' & fig.2j.uniq$motif == 'SP family',]$genes,]

df.seq.twist.rep = df.seq[df.seq$re %in% fig.2j.uniq[fig.2j.uniq$status == 'Repressed' & fig.2j.uniq$motif == 'TWIST family',]$genes,]
df.seq.twist.act = df.seq[df.seq$re %in% fig.2j.uniq[fig.2j.uniq$status == 'Activated' & fig.2j.uniq$motif == 'TWIST family',]$genes,]

df.seq.bhlh.rep = df.seq[df.seq$re %in% fig.2j.uniq[fig.2j.uniq$status == 'Repressed' & fig.2j.uniq$motif == 'bHLH family',]$genes,]
df.seq.bhlh.act = df.seq[df.seq$re %in% fig.2j.uniq[fig.2j.uniq$status == 'Activated' & fig.2j.uniq$motif == 'bHLH family',]$genes,]


toupper(df.seq.bhlh.act$sequence)

library(ggseqlogo)

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



pswm.fimo(df.seq.ap1.rep , out = 'AP1_repressed')
pswm.fimo(df.seq.ap1.act , out = 'AP1_activated')
pswm.fimo(df.seq.gr.rep , out = 'GR_repressed')
pswm.fimo(df.seq.gr.act , out = 'GR_activated')
pswm.fimo(df.seq.sp.rep , out = 'SP_repressed')
pswm.fimo(df.seq.sp.act , out = 'SP_activated')
pswm.fimo(df.seq.twist.rep , out = 'TWIST_repressed')
pswm.fimo(df.seq.twist.act , out = 'TWIST_activated')
pswm.fimo(df.seq.bhlh.rep , out = 'bHLH_repressed')
pswm.fimo(df.seq.bhlh.act , out = 'bHLH_activated')



# in SHELL
ceqlogo -i AP1_activated.txt -m AP1_activated -o AP1_activated.eps
ceqlogo -i AP1_repressed.txt -m AP1_repressed -o AP1_repressed.eps
ceqlogo -i GR_activated.txt -m GR_activated -o GR_activated.eps
ceqlogo -i GR_repressed.txt -m GR_repressed -o GR_repressed.eps
ceqlogo -i SP_activated.txt -m SP_activated -o SP_activated.eps
ceqlogo -i SP_repressed.txt -m SP_repressed -o SP_repressed.eps
ceqlogo -i TWIST_activated.txt -m TWIST_activated -o TWIST_activated.eps
ceqlogo -i TWIST_repressed.txt -m TWIST_repressed -o TWIST_repressed.eps
ceqlogo -i bHLH_activated.txt -m bHLH_activated -o bHLH_activated.eps
ceqlogo -i bHLH_repressed.txt -m bHLH_repressed -o bHLH_repressed.eps


plot.seqlogo.func <- function(x, outfile = "fimo_out.pdf") {
  w =  0.663 + (ncol(x) + 1)*0.018 + (ncol(x)+2)* .336
  pdf(outfile, useDingbats=FALSE, width=w, height=2.695)
  print(ggseqlogo(x,  facet = "wrap", font = 'helvetica_bold'))
  dev.off()
}


plot.seqlogo.func(pswm.fimo(df.seq.bhlh.act, out = 'bHLH_activated'))


#


x = read.table('/Users/guertinlab/Desktop/fimo_stuff/SP_family_fimo.bed', stringsAsFactors=FALSE)
x = x[x[,6] != -1,]

a <- apply(x[,c(1:9) ], 1 ,paste0, collapse = ":" )
a = gsub(" ", "", a, fixed = TRUE)
b = x[,10]
for (i in 1:length(a)) {
    write.table(paste('>', a[i], sep =''), file = 'sp_fimo.txt', append=TRUE, col.names=FALSE, row.names=FALSE, sep ='', quote= FALSE)
    write.table(b[i], file = 'sp_fimo.txt', append=TRUE, col.names=FALSE, row.names=FALSE, sep ='', quote= FALSE)
}

                                        #shell
#mast -oc SP1.mast_output -ev 0.01 -hit_list -comp -best /Users/guertinlab/Desktop/fimo_stuff/SP_meme.txt /Users/guertinlab/Desktop/fimo_stuff/sp_fimo.txt > output_sp1.txt
#mast -oc KLF.mast_output -ev 0.01 -hit_list -comp -best /Users/guertinlab/Desktop/fimo_stuff/KLF_meme.txt /Users/guertinlab/Desktop/fimo_stuff/sp_fimo.txt > output_klf.txt

fimo --thresh 0.01 --text /Users/guertinlab/Desktop/fimo_stuff/SP_meme.txt /Users/guertinlab/Desktop/fimo_stuff/sp_fimo.txt > output_sp1.txt
fimo --thresh 0.01 --text /Users/guertinlab/Desktop/fimo_stuff/KLF_meme.txt /Users/guertinlab/Desktop/fimo_stuff/sp_fimo.txt > output_klf.txt


sp.mast = read.table('output_sp1.txt', stringsAsFactors=FALSE, sep = '\t', header = TRUE)
klf.mast = read.table('output_klf.txt', stringsAsFactors=FALSE, sep = '\t', header = TRUE)

sp.mast = sp.mast[sp.mast[,6] != -1,]
sp.mast = aggregate(as.numeric(score)~sequence_name, data=sp.mast, FUN=max)
rownames(sp.mast) = sp.mast[,1]

klf.mast = klf.mast[klf.mast[,6] != -1,]
klf.mast = aggregate(as.numeric(score)~sequence_name, data=klf.mast, FUN=max)
rownames(klf.mast) = klf.mast[,1]

sp.data = merge(sp.mast, klf.mast, by = 1)
colnames(sp.data) = c('seq', 'sp', 'klf')

xyplot(sp.data$sp ~ sp.data$klf, pch = 16, cex = 0.5, xlim = c(-5,25), ylim = c(-5,25), aspect = 1, panel = function(x, y, ...) {
         panel.xyplot(x, y, ...)
         panel.abline(0, 1)
})

sp.motifs.lda = aggregate(as.numeric(score)~X1+X2+X3, data=all.sp.motifs, FUN=max)
klf.motifs.lda = aggregate(as.numeric(score)~X1+X2+X3, data=all.klf.motifs, FUN=max)
colnames(sp.motifs.lda) = c('chr', 'start', 'end', 'sp')
colnames(klf.motifs.lda) = c('chr', 'start', 'end', 'klf')


rownames(klf.motifs.lda) = paste0(klf.motifs.lda[,1], ':', klf.motifs.lda[,2], '-', klf.motifs.lda[,3])
rownames(sp.motifs.lda) = paste0(sp.motifs.lda[,1], ':', sp.motifs.lda[,2], '-', sp.motifs.lda[,3])




#this could be the input to LDA
chrom = vector(mode="character", length = length(unique(plot.df.atac$genes)))
start = vector(mode="integer", length = length(unique(plot.df.atac$genes)))
end = vector(mode="integer", length = length(unique(plot.df.atac$genes)))
cls = vector(mode="integer", length = length(unique(plot.df.atac$genes)))

count = 0
for (j in unique(plot.df.atac$genes)) {
    count = count +1
    chrom[count] = strsplit(j, ":")[[1]][1]
    start[count] = as.numeric(strsplit(strsplit(j, ":")[[1]][2], "-")[[1]][1])
    end[count] = as.numeric(strsplit(strsplit(j, ":")[[1]][2], "-")[[1]][2])    
}

res.lda = data.frame(cbind(chrom, start, end))
colnames(res.lda) = c('chr', 'start', 'end')
res.lda$start = as.numeric(as.character(res.lda$start))
res.lda$end = as.numeric(as.character(res.lda$end))
rownames(res.lda) = paste0(res.lda[,1], ':', res.lda[,2], '-', res.lda[,3])

res.lda = merge(res.lda, klf.motifs.lda, by="row.names", all = TRUE)
rnames = res.lda[,1]

res.lda = data.frame(res.lda[,c(2:4, 8)])
rownames(res.lda) = rnames

                                        #stopped here:
res.lda = merge(res.lda, sp.motifs.lda, by="row.names", all = TRUE)
rnames = res.lda[,1]

res.lda = data.frame(res.lda[,c(2:5, 9)])

rownames(res.lda) = rnames

save(res.lda, file = 'res.lda.Rdata')

compare.lda = plot.df.atac[!duplicated(plot.df.atac$genes),c(1,22)]
rownames(compare.lda) = compare.lda[,1]


res.lda = merge(res.lda, compare.lda, by="row.names", all = TRUE)


rownames(res.lda) = res.lda[,1]
res.lda = res.lda[,c(5,6,8)]

res.lda = res.lda[!is.na(res.lda$sp) & !is.na(res.lda$sp),]

xyplot(sp ~ klf, groups = status, data= res.lda, pch = 16, col = c('red', 'blue'),
       panel = function(x, y, ...) {
         panel.xyplot(x, y, ...)
       })

library(MASS)

fit <- lda(status ~ sp + klf, data= res.lda)

plot(fit)

ct <- table(res.lda$status, fit$class)
diag(prop.table(ct, 1))
sum(diag(prop.table(ct)))


gmean <- fit$prior%*%fit$means
const <- drop(gmean%*%fit$scaling)

slope <- -fit$scaling[1]/fit$scaling[2]
intercept <- const/fit$scaling[2]

slope
#7.794952
intercept
#-76.05604

res.lda$status = as.factor(res.lda$status)


relevel(res.lda$status, "Repressed")

pdf(file = paste0('LDA_analysis.pdf'),width=3,height=3)

print(
    
    xyplot(klf ~ sp, groups = status, data= res.lda, pch = 16, cex = 0.6,
           col = c('#FF000066', '#0000FF40'),
#       ylim = c(-5, 25),
#       xlim = c(-5, 25),
       xlab = "SP score",
       ylab = "KLF score",
       scales=list(alternating=c(1,0)),
       panel = function(x, y, ...) {
         panel.xyplot(x, y, ...)
         panel.abline(-76.05604, 7.794952, lty = 2, col = 'grey45', lwd = 2.3)
       })

)
dev.off()






                                        #get_unique_peaks


out.klf <- strsplit(as.character(klf.mast[,1]), ':') 
all.klf.motifs <- data.frame(do.call(rbind, out.klf, quote=FALSE), klf.mast[,c(2)])
colnames(all.klf.motifs) = c(colnames(all.klf.motifs)[1:9], 'score')

out.sp <- strsplit(as.character(sp.mast[,1]), ':') 
all.sp.motifs <- data.frame(do.call(rbind, out.sp, quote=FALSE), sp.mast[,c(2)])
colnames(all.sp.motifs) = c(colnames(all.sp.motifs)[1:9], 'score')


klf.motifs = all.klf.motifs[(all.klf.motifs$score > (all.sp.motifs$score*slope + intercept)),]
sp.motifs = all.sp.motifs[(all.klf.motifs$score <= (all.sp.motifs$score*slope + intercept)),]


plot(all.sp.motifs[(all.klf.motifs$score > (all.sp.motifs$score*slope + intercept)),]$score,
     all.klf.motifs[(all.klf.motifs$score > (all.sp.motifs$score*slope + intercept)),]$score)

plot(all.sp.motifs[(all.klf.motifs$score <= (all.sp.motifs$score*slope + intercept)),]$score,
     all.klf.motifs[(all.klf.motifs$score <= (all.sp.motifs$score*slope + intercept)),]$score)


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

res = data.frame(cbind(chrom, start, end))
colnames(res) = c('chr', 'start', 'end')
res$start = as.numeric(as.character(res$start))
res$end = as.numeric(as.character(res$end))
rownames(res) = paste0(res[,1], ':', res[,2], '-', res[,3])

sp.motifs.1 = aggregate(as.numeric(score)~X1+X2+X3, data=sp.motifs, FUN=sum)
klf.motifs.1 = aggregate(as.numeric(score)~X1+X2+X3, data=klf.motifs, FUN=sum)


rownames(sp.motifs.1) = paste0(sp.motifs.1[,1], ':', sp.motifs.1[,2], '-', sp.motifs.1[,3])
rownames(klf.motifs.1) = paste0(klf.motifs.1[,1], ':', klf.motifs.1[,2], '-', klf.motifs.1[,3])

colnames(sp.motifs.1) = c('chr', 'start', 'end', 'sp')
colnames(klf.motifs.1) = c('chr', 'start', 'end', 'klf')

res = merge(res, sp.motifs.1, by="row.names", all = TRUE)
rnames = res[,1]

res = data.frame(res[,c(2:4, 8)])

colnames(res) = c('chr', 'start','end', 'SP_alone')
rownames(res) = rnames

print(head(klf.motifs.1))

res = merge(res, klf.motifs.1, by="row.names", all = TRUE)
rnames = res[,1]


res = data.frame(res[,c(2:5, 9)])

colnames(res) = c('chr', 'start','end', 'SP_alone', 'KLF_alone')
rownames(res) = rnames

sp.klf.scores.atac = res
save(sp.klf.scores.atac, file = 'sp.klf.scores.atac.Rdata')

x = plot.df.atac

for (i in 1:nrow(plot.df.atac)) {
    re = plot.df.atac$genes[i] 
    plot.df.atac[i, 23:24] = sp.klf.scores.atac[re,c(4:5)]
}

plot.df.atac$SP_alone_T <- ifelse(is.na(plot.df.atac$SP_alone), FALSE, TRUE)
plot.df.atac$KLF_alone_T <- ifelse(is.na(plot.df.atac$KLF_alone), FALSE, TRUE)


save(plot.df.atac, file = 'plot.df.atac.MJG.Rdata')


sp1.plot = plot.df.atac[plot.df.atac$AP1 == FALSE & plot.df.atac$GR == FALSE & plot.df.atac$SP == TRUE & plot.df.atac$TWIST == FALSE & plot.df.atac$bHLH == FALSE & plot.df.atac$SP_alone_T == TRUE & plot.df.atac$KLF_alone_T == FALSE,]
sp1.plot$motif = 'SP family'
klf.plot = plot.df.atac[plot.df.atac$AP1 == FALSE & plot.df.atac$GR == FALSE & plot.df.atac$SP == TRUE & plot.df.atac$TWIST == FALSE & plot.df.atac$bHLH == FALSE & plot.df.atac$SP_alone_T == FALSE & plot.df.atac$KLF_alone_T == TRUE,]
klf.plot$motif = 'KLF family'

fig.2k = rbind(sp1.plot, klf.plot)

#added these value in AI
nrow(sp1.plot[sp1.plot$status == 'Activated',])/7
nrow(sp1.plot[sp1.plot$status == 'Repressed',])/7

nrow(klf.plot[klf.plot$status == 'Activated',])/7
nrow(klf.plot[klf.plot$status == 'Repressed',])/7


#pdf(file = paste0('partition_sp_klf.pdf'),width=3,height=3)
#
#print(xyplot(sp.data$sp ~ sp.data$klf, pch = 16, cex = 0.5, xlim = c(-5,25),
#             ylim = c(-5,25), aspect = 1, col = '#8e8e8e4D',
#             ylab = 'SP1 Score', xlab = 'KLF Score',  
#             panel = function(x, y, ...) {
#         panel.xyplot(x, y, ...)
#         panel.abline(0, 1, col = 'red')
#})
#)
#dev.off()


pdf(file = paste0('2k.traces.pdf'),width=4,height=4)
trellis.par.set(box.umbrella = list(lty = 1, col="black", lwd=1),
                box.rectangle = list( lwd=1.0, col="black", alpha = 1.0),
                plot.symbol = list(col="black", lwd=1.0, pch ='.'))
print(
    xyplot(value ~  sample.conditions | motif, group = genes, data = fig.2k, type = c('l'),#type = c('l','p'),
           scales=list(x=list(cex=1.0,relation = "free", rot = 45), y = list(cex=1.0, relation="free")),
           aspect=1.0,
           between=list(y=0.5, x=0.5),
           ylab = list(label = 'Normalized ATAC signal', cex =1.0),
#           index.cond=list(c(1,3,4,5,2)),
           xlab = list(label = 'Time (minutes)', cex =1.0),
           par.settings = list(superpose.symbol = list(pch = c(16),col=c('grey20'), cex =0.5),
                               #I want to change the background strip to the corresponding motif color
                               strip.background=list(col=c("grey80", "green", "purple", "yellow","teal")),
                               superpose.line = list(col = as.character(cat.colours$colour), lwd=c(1),lty = c(1))),
)
)
dev.off()


fig.2k.uniq = fig.2k[!duplicated(paste0(fig.2k$genes, fig.2k$motif)),]

                                        #

fig.2k.uniq$new[fig.2k.uniq$KLF_alone_T == TRUE] <- fig.2k.uniq$KLF_alone[fig.2k.uniq$KLF_alone_T == TRUE]
fig.2k.uniq$new[fig.2k.uniq$SP_alone_T == TRUE] <- fig.2k.uniq$SP_alone[fig.2k.uniq$SP_alone_T == TRUE]

fig.2k.uniq$new <- as.numeric(as.character(fig.2k.uniq$new))


pdf(file = paste0('2k.bw.pdf'),width=4,height=4)
trellis.par.set(box.umbrella = list(lty = 1, col=c("red", "blue"), lwd=2),
                box.rectangle = list( lwd=2.0, col=c("red", "blue"), alpha = 1.0),
                plot.symbol = list(col=c("red", "blue"), lwd=2.0, pch ='.'))
print(
    
bwplot(log(as.numeric(as.character(new)), base =10) ~ status | motif, data = fig.2k.uniq, horizontal = FALSE, pch = '|', do.out = FALSE,
       scales=list(x=list(cex=1.0, relation = "free", rot = 45), y = list(cex=1.0, relation="free")),
       aspect=2.0,
       between=list(y=0.5, x=0.5),
#       index.cond=list(c(1,3,4,5,2)),
       ylab = expression('log'[10]* paste(Sigma,'(motif scores)')),
       xlab = expression('ATAC Peak category'),
#manually settign to avoid outlier, since do.out = FALSE       
       ylim = list(c(0.8, 1.4),  c(1, 2.1)),
       par.settings = list(superpose.symbol = list(pch = c(16),col=c('grey20'), cex =0.5),
                               #I want to change the background strip to the corresponding motif color
                               strip.background=list(col=c("grey80"))))

)
dev.off()



ap1 = plot.df.atac[plot.df.atac$AP1 == TRUE & plot.df.atac$GR == FALSE & plot.df.atac$SP_alone_T== FALSE & plot.df.atac$TWIST == FALSE & plot.df.atac$bHLH == FALSE & plot.df.atac$KLF_alone_T == FALSE,]
ap1$motif = 'AP1 family'
gr = plot.df.atac[plot.df.atac$AP1 == FALSE & plot.df.atac$GR == TRUE & plot.df.atac$SP_alone_T == FALSE & plot.df.atac$TWIST == FALSE & plot.df.atac$bHLH == FALSE & plot.df.atac$KLF_alone_T == FALSE,]
gr$motif = 'GR family'
sp = plot.df.atac[plot.df.atac$AP1 == FALSE & plot.df.atac$GR == FALSE & plot.df.atac$SP_alone_T == TRUE & plot.df.atac$TWIST == FALSE & plot.df.atac$bHLH == FALSE & plot.df.atac$KLF_alone_T == FALSE,]
sp$motif = 'SP family'
klf = plot.df.atac[plot.df.atac$AP1 == FALSE & plot.df.atac$GR == FALSE & plot.df.atac$SP_alone_T == FALSE & plot.df.atac$TWIST == FALSE & plot.df.atac$bHLH == FALSE & plot.df.atac$KLF_alone_T == TRUE,]
klf$motif = 'KLF family'

twist = plot.df.atac[plot.df.atac$AP1 == FALSE & plot.df.atac$GR == FALSE & plot.df.atac$SP_alone_T == FALSE & plot.df.atac$TWIST == TRUE & plot.df.atac$bHLH == FALSE & plot.df.atac$KLF_alone_T == FALSE,]
twist$motif = 'TWIST family'
bhlh = plot.df.atac[plot.df.atac$AP1 == FALSE & plot.df.atac$GR == FALSE & plot.df.atac$SP_alone_T == FALSE & plot.df.atac$TWIST == FALSE & plot.df.atac$bHLH == TRUE & plot.df.atac$KLF_alone_T == FALSE,]
bhlh$motif = 'bHLH family'



fig.2j.concise = rbind(ap1, gr, sp, klf, twist, bhlh)

#added these value in AI
nrow(ap1[ap1$status == 'Activated',])/7
nrow(ap1[ap1$status == 'Repressed',])/7

nrow(gr[gr$status == 'Activated',])/7
nrow(gr[gr$status == 'Repressed',])/7

nrow(sp[sp$status == 'Activated',])/7
nrow(sp[sp$status == 'Repressed',])/7

nrow(klf[klf$status == 'Activated',])/7
nrow(klf[klf$status == 'Repressed',])/7


nrow(twist[twist$status == 'Activated',])/7
nrow(twist[twist$status == 'Repressed',])/7

nrow(bhlh[bhlh$status == 'Activated',])/7
nrow(bhlh[bhlh$status == 'Repressed',])/7


pdf(file = paste0('2j.traces.concise.pdf'),width=12,height=4)
trellis.par.set(box.umbrella = list(lty = 1, col="black", lwd=1),
                box.rectangle = list( lwd=1.0, col="black", alpha = 1.0),
                plot.symbol = list(col="black", lwd=1.0, pch ='.'))
print(
    xyplot(value ~  sample.conditions | motif, group = genes, data = fig.2j.concise, type = c('l'),#type = c('l','p'),
           scales=list(x=list(cex=1.0,relation = "free", rot = 45), y = list(cex=1.0, relation="free")),
           aspect=1.0,
           between=list(y=0.5, x=0.5),
           ylab = list(label = 'Normalized ATAC signal', cex =1.0),
           index.cond=list(c(1,3,4,5,2, 6)),
           xlab = list(label = 'Time (minutes)', cex =1.0),
           par.settings = list(superpose.symbol = list(pch = c(16),col=c('grey20'), cex =0.5),
                               #I want to change the background strip to the corresponding motif color
                               strip.background=list(col=c("grey80", "green", "purple", "yellow","teal")),
                               superpose.line = list(col = as.character(cat.colours$colour), lwd=c(1),lty = c(1))),
)
)
dev.off()



fig.2j.uniq = fig.2j.concise[!duplicated(paste0(fig.2j.concise$genes, fig.2j.concise$motif)),]

                                        #

fig.2j.uniq$new[fig.2j.uniq$AP1 == TRUE] <- fig.2j.uniq$AP1_family[fig.2j.uniq$AP1 == TRUE]
fig.2j.uniq$new[fig.2j.uniq$GR == TRUE] <- fig.2j.uniq$GR_family[fig.2j.uniq$GR == TRUE]
fig.2j.uniq$new[fig.2j.uniq$SP_alone_T == TRUE] <- fig.2j.uniq$SP_alone[fig.2j.uniq$SP_alone_T  == TRUE]
fig.2j.uniq$new[fig.2j.uniq$KLF_alone_T == TRUE] <- fig.2j.uniq$KLF_alone[fig.2j.uniq$KLF_alone_T  == TRUE]

fig.2j.uniq$new[fig.2j.uniq$TWIST == TRUE] <- fig.2j.uniq$TWIST_family[fig.2j.uniq$TWIST == TRUE]
fig.2j.uniq$new[fig.2j.uniq$bHLH == TRUE] <- fig.2j.uniq$bHLH_family[fig.2j.uniq$bHLH == TRUE]

fig.2j.uniq$new <- as.numeric(as.character(fig.2j.uniq$new))



pdf(file = paste0('2j.bw.concise.pdf'),width=12,height=4)
trellis.par.set(box.umbrella = list(lty = 1, col=c("red", "blue"), lwd=2),
                box.rectangle = list( lwd=2.0, col=c("red", "blue"), alpha = 1.0),
                plot.symbol = list(col=c("red", "blue"), lwd=2.0, pch ='.'))
print(
    
bwplot(log(as.numeric(as.character(new)), base =10) ~ status | motif, data = fig.2j.uniq, horizontal = FALSE, pch = '|', do.out = FALSE,
       scales=list(x=list(cex=1.0, relation = "free", rot = 45), y = list(cex=1.0, relation="free")),
       aspect=2.0,
       between=list(y=0.5, x=0.5),
       index.cond=list(c(1,3,4,5,2, 6)),
       ylab = expression('log'[10]* paste(Sigma,'(motif scores)')),
       xlab = expression('ATAC Peak category'),
#manually settign to avoid outlier, since do.out = FALSE       
#       ylim = list(c(0.9, 1.9),  c(0.9, 1.7), c(0.8, 1.8), c(1, 2.15),  c(0.96, 1.4)),
       par.settings = list(superpose.symbol = list(pch = c(16),col=c('grey20'), cex =0.5),
                               #I want to change the background strip to the corresponding motif color
                               strip.background=list(col=c("grey80"))))

)
dev.off()


                                        #need to plot seqLogo


count = -1
vec.names = c()

df.seq = as.data.frame(matrix(ncol=5, nrow=0),stringsAsFactors = FALSE)
colnames(df.seq) = c('chr', 'start', 'end', 'sequence', 'factor')

for(bed.file in Sys.glob(file.path('/Users/guertinlab/Desktop/fimo_stuff/*family_fimo.bed'))) {
    print(bed.file)
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

                                        #df.seq has the sequences

df.seq.ap1.rep = df.seq[df.seq$re %in% fig.2j.uniq[fig.2j.uniq$status == 'Repressed' & fig.2j.uniq$motif == 'AP1 family',]$genes,]
df.seq.ap1.act = df.seq[df.seq$re %in% fig.2j.uniq[fig.2j.uniq$status == 'Activated' & fig.2j.uniq$motif == 'AP1 family',]$genes,]

df.seq.gr.rep = df.seq[df.seq$re %in% fig.2j.uniq[fig.2j.uniq$status == 'Repressed' & fig.2j.uniq$motif == 'GR family',]$genes,]
df.seq.gr.act = df.seq[df.seq$re %in% fig.2j.uniq[fig.2j.uniq$status == 'Activated' & fig.2j.uniq$motif == 'GR family',]$genes,]

df.seq.sp.rep = df.seq[df.seq$re %in% fig.2j.uniq[fig.2j.uniq$status == 'Repressed' & fig.2j.uniq$motif == 'SP family',]$genes,]
df.seq.sp.act = df.seq[df.seq$re %in% fig.2j.uniq[fig.2j.uniq$status == 'Activated' & fig.2j.uniq$motif == 'SP family',]$genes,]


df.seq.klf.rep = df.seq[df.seq$re %in% fig.2j.uniq[fig.2j.uniq$status == 'Repressed' & fig.2j.uniq$motif == 'KLF family',]$genes,]
df.seq.klf.act = df.seq[df.seq$re %in% fig.2j.uniq[fig.2j.uniq$status == 'Activated' & fig.2j.uniq$motif == 'KLF family',]$genes,]



df.seq.twist.rep = df.seq[df.seq$re %in% fig.2j.uniq[fig.2j.uniq$status == 'Repressed' & fig.2j.uniq$motif == 'TWIST family',]$genes,]
df.seq.twist.act = df.seq[df.seq$re %in% fig.2j.uniq[fig.2j.uniq$status == 'Activated' & fig.2j.uniq$motif == 'TWIST family',]$genes,]

df.seq.bhlh.rep = df.seq[df.seq$re %in% fig.2j.uniq[fig.2j.uniq$status == 'Repressed' & fig.2j.uniq$motif == 'bHLH family',]$genes,]
df.seq.bhlh.act = df.seq[df.seq$re %in% fig.2j.uniq[fig.2j.uniq$status == 'Activated' & fig.2j.uniq$motif == 'bHLH family',]$genes,]


toupper(df.seq.bhlh.act$sequence)


pswm.fimo(df.seq.ap1.rep , out = 'AP1_repressed')
pswm.fimo(df.seq.ap1.act , out = 'AP1_activated')
pswm.fimo(df.seq.gr.rep , out = 'GR_repressed')
pswm.fimo(df.seq.gr.act , out = 'GR_activated')
pswm.fimo(df.seq.sp.rep , out = 'SP_repressed')
pswm.fimo(df.seq.sp.act , out = 'SP_activated')
pswm.fimo(df.seq.twist.rep , out = 'TWIST_repressed')
pswm.fimo(df.seq.twist.act , out = 'TWIST_activated')
pswm.fimo(df.seq.bhlh.rep , out = 'bHLH_repressed')
pswm.fimo(df.seq.bhlh.act , out = 'bHLH_activated')
pswm.fimo(df.seq.klf.rep , out = 'KLF_repressed')
pswm.fimo(df.seq.klf.act , out = 'KLF_activated')


# in SHELL
ceqlogo -i AP1_activated.txt -m AP1_activated -o AP1_activated.eps
ceqlogo -i AP1_repressed.txt -m AP1_repressed -o AP1_repressed.eps
ceqlogo -i GR_activated.txt -m GR_activated -o GR_activated.eps
ceqlogo -i GR_repressed.txt -m GR_repressed -o GR_repressed.eps
ceqlogo -i SP_activated.txt -m SP_activated -o SP_activated.eps
ceqlogo -i SP_repressed.txt -m SP_repressed -o SP_repressed.eps
ceqlogo -i TWIST_activated.txt -m TWIST_activated -o TWIST_activated.eps
ceqlogo -i TWIST_repressed.txt -m TWIST_repressed -o TWIST_repressed.eps
ceqlogo -i bHLH_activated.txt -m bHLH_activated -o bHLH_activated.eps
ceqlogo -i bHLH_repressed.txt -m bHLH_repressed -o bHLH_repressed.eps
ceqlogo -i KLF_activated.txt -m KLF_activated -o KLF_activated.eps
ceqlogo -i KLF_repressed.txt -m KLF_repressed -o KLF_repressed.eps
