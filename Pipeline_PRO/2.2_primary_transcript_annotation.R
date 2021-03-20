library(devtools)
library(NMF)
library(dplyr)
library(bigWig)
library(pracma)
library(RColorBrewer)
#install_github("WarrenDavidAnderson/genomicsRpackage/primaryTranscriptAnnotation")
library(primaryTranscriptAnnotation)

setwd('/scratch/bhn9by/PRO/primary_transcript_annotation')

source('pTA.functions.R')

# import data for first exons, annotate, and remove duplicate transcripts
fname = "../gencode.mm10.firstExon.bed"
dat0 = read.table(fname,header=F,stringsAsFactors=F)
names(dat0) = c('chr', 'start', 'end', 'gene', 'xy', 'strand')
dat0 = unique(dat0)
gencode.firstExon = dat0
# import data for all transcripts, annotate, and remove duplicate transcripts
fname = "../gencode.mm10.transcript.bed"
dat0 = read.table(fname,header=F,stringsAsFactors=F)
names(dat0) = c('chr', 'start', 'end', 'gene', 'xy', 'strand')
dat0 = unique(dat0)
gencode.transcript = dat0
# chromosome sizes
chrom.sizes = read.table("../mm10.chrom.sizes",stringsAsFactors=F,header=F)
names(chrom.sizes) = c("chr","size")

plus.file = "../pro_plus_merged.bigWig"
minus.file = "../pro_minus_merged.bigWig"
bw.plus = load.bigWig(plus.file)
bw.minus = load.bigWig(minus.file)

#get intervals for furthest TSS and TTS +/- interval
largest.interval.bed = get.largest.interval(bed=gencode.transcript)
# filtering which read count and read density might be consider as insignificant
transcript.reads = read.count.transcript(bed=gencode.transcript, bw.plus=bw.plus, bw.minus=bw.minus)

# evaluate count and density distributions
pdf("Histogram_density_distibution.pdf", useDingbats = FALSE, width=6.4, height=5.4)
par(mfrow=c(1,2))
hist(log(transcript.reads$density), breaks=200,
col="black",xlab="log10 read density",main="")
hist(log(transcript.reads$counts), breaks=200,
col="black",xlab="log10 read count",main="")
dev.off()

# specify which genes to cut based on low expression, visualize cutoffs
pdf("Histogram_density_distibution_treshold.pdf", useDingbats = FALSE, width=6.4, height=5.4)
den.cut = -4
cnt.cut = 2
ind.cut.den = which(log(transcript.reads$density) < den.cut)
ind.cut.cnt = which(log(transcript.reads$counts) < cnt.cut)
ind.cut = union(ind.cut.den, ind.cut.cnt)
par(mfrow=c(1,2))
hist(log(transcript.reads$density), breaks=200,
col="black",xlab="log"[10]~" read density",main="")
abline(v=den.cut, col="red")
hist(log(transcript.reads$counts), breaks=200,
col="black",xlab="log"[10]~" read count",main="")
abline(v=cnt.cut, col="red")
dev.off()

# remove "unexpressed" genes
unexp = names(transcript.reads$counts)[ind.cut]
largest.interval.expr.bed = largest.interval.bed[
    !(largest.interval.bed$gene %in% unexp),]
# select the TSS for each gene and incorporate these TSSs
# into the largest interval coordinates
bp.range = c(20,120)
cnt.thresh = 2
bed.out = largest.interval.expr.bed
bed.in = gencode.firstExon[gencode.firstExon$gene %in% bed.out$gene,]
TSS.gene = get.TSS(bed.in=bed.in, bed.out=bed.out,
                   bw.plus=bw.plus, bw.minus=bw.minus,
                   bp.range=bp.range, cnt.thresh=cnt.thresh)
TSS.gene = TSS.gene$bed

# parameters and analysis
window = 1000
bp.bin = 10
bed1 = TSS.gene
bed2 = largest.interval.expr.bed
bed2 = bed2[bed2$gene %in% bed1$gene,]
tss.eval = eval.tss(bed1=bed1, bed2=bed2,
                    bw.plus=bw.plus, bw.minus=bw.minus,
window=window, bp.bin=bp.bin, fname="TSSres.pdf")

                                        #Inferred TSSs showed as distribution of paused RNA polymerase density
pdf("TS.pdf", useDingbats = FALSE, width=6.4, height=5.4)
bk = seq(-window/2,window/2, bp.bin)
hist(tss.eval$tss.dists.inf$dist, main="overlay",
     xlab="dist from TSS to read max (bp)",
     breaks=bk, col='black')
bk = seq(-window/2,window/2, bp.bin)
hist( tss.eval$tss.dists.lng$dist,
breaks=bk, col='red', add=T)
dev.off()

# identify duplicates
dups = get.dups(bed = TSS.gene)
max(dups$cases)
head(dups)
write.table(dups,"duplicateCoords.bed",quote=F,sep="\t",col.names=F,row.names=F)

#put *merged.bigWig files on cyverse and look at duplicate examples in browser
#annotations that should be entirely removed go into 'duplicate_remove.csv'
#annotations that should be kept should be ordered and put into 'duplicate_keepadjacent.csv'
remove.genes.id <- read.csv(file="duplicate_remove.csv",
header=TRUE, sep=",", stringsAsFactors=F)
fix.genes.id <- read.csv("duplicate_keepadjacent.csv",
header= TRUE, sep=",", stringsAsFactors=F)
# remove genes based on manual analysis, genes cosidered as adjecent will be adressed below
genes.remove1 = c(remove.genes.id$remove, fix.genes.id$upstream, fix.genes.id$downstream)
TSS.gene.filtered1 = TSS.gene[!(TSS.gene$gene %in% genes.remove1),]
# confirm that there are no any overlaps
#this object should have no data
dups1 =get.dups(bed = TSS.gene.filtered1)

# overlap analysis
overlap.data = gene.overlaps(bed = TSS.gene.filtered1)
has.start.inside = overlap.data$has.start.inside
is.a.start.inside = overlap.data$is.a.start.inside
head(has.start.inside)
dim(has.start.inside) #271
head(is.a.start.inside)
dim(is.a.start.inside) #323
head(overlap.data$cases)
dim(overlap.data$cases) #503

# Checking those genes in overlap.gene lists
overlap.genes = overlap.data$cases$gene %>% unique
length(grep("^AC[0-9][0-9]",overlap.genes)) #7
length(grep("^AL[0-9][0-9]",overlap.genes)) #3
length(grep("^AP[0-9][0-9]",overlap.genes)) #0
# set to remove usually poorly annotated genes in overlapping regions
genes.remove2 = overlap.genes[grep("^AC[0-9][0-9]",overlap.genes)]
genes.remove2 = c(genes.remove2, overlap.genes[grep("^AL[0-9][0-9]",overlap.genes)])
genes.remove2 = c(genes.remove2, overlap.genes[grep("^AP[0-9][0-9]",overlap.genes)])
# identify genes with multiple starts inside (i.e. 'big' genes)
mult.inside.starts = inside.starts(vec = is.a.start.inside$xy)
length(mult.inside.starts) #31
# add genes with multiple starts inside to remove
genes.remove2 = c(genes.remove2, mult.inside.starts %>% unique) #41
# remove filtered genes and re-run overlap analysis
in.dat = TSS.gene.filtered1[!(TSS.gene.filtered1$gene %in% genes.remove2),]
overlap.data = gene.overlaps( bed = in.dat )
has.start.inside = overlap.data$has.start.inside
is.a.start.inside = overlap.data$is.a.start.inside
case.dat = overlap.data$cases
length(unique(case.dat$cases)) #171
# data for manual analysis and curration
fname = "nonid_overlaps.txt"
write.table(case.dat,fname,col.names=T,row.names=F,sep="\t",quote=F)
overlaps = case.dat %>% select(chr,start,end,gene,cases)
write.table(overlaps,"nonid_overlaps.bed",quote=F,sep="\t",col.names=F,row.names=F)

write.table(TSS.gene.filtered1,"TSS.gene.filtered1.bed",quote=F,sep="\t",col.names=F,row.names=F)

remove.genes.ov = read.csv("overlaps_remove.csv", 
header=TRUE, stringsAsFactors=F)
fix.genes.ov =read.csv("overlaps_keepadjacent.csv", 
header=TRUE, stringsAsFactors=F)
genes.remove3 = c(genes.remove2, remove.genes.ov$remove,
fix.genes.ov$upstream, fix.genes.ov$downstream)
# read corrected start/end for 10 genes in bed file

TSS.gene.filtered2 = TSS.gene.filtered1[
!(TSS.gene.filtered1$gene %in% genes.remove3),]
# verify the abscence of overlaps
overlap.data = gene.overlaps( bed = TSS.gene.filtered2 )
overlap.data$cases

fix.genes.ov =read.csv("overlaps_keepadjacent.csv", 
header=TRUE, stringsAsFactors=F)
# get the coordinates for adjacent gene pairs
fix.genes = rbind(fix.genes.id, fix.genes.ov)
bp.bin = 5
knot.div = 40
shift.up = 100
delta.tss = 50
diff.tss = 1000
dist.from.start = 50
adjacent.coords = adjacent.gene.coords(fix.genes=fix.genes, bed.long=TSS.gene,
exon1=gencode.firstExon,
bw.plus=bw.plus, bw.minus=bw.minus,
knot.div=knot.div, bp.bin=bp.bin,
shift.up=shift.up, delta.tss=delta.tss,
dist.from.start=dist.from.start,
diff.tss=diff.tss, fname="adjacentSplines.pdf")

# visualize coordinates for a specific pair
adjacent.coords.plot(adjacent.coords=adjacent.coords,
pair=fix.genes[1,],
bw.plus=bw.plus,
bw.minus=bw.minus)

# aggregate downstream adjacent genes with main data
TSS.gene.filtered3 = rbind(TSS.gene.filtered2, adjacent.coords)
overlap.data = gene.overlaps( bed = TSS.gene.filtered3 )
overlap.data$cases

#get intervals for TTS evaluation
add.to.end = 100000
fraction.end = 0.2
dist.from.start = 50
bed.for.tts.eval = get.end.intervals(bed=TSS.gene.filtered3,
add.to.end=add.to.end,
fraction.end=fraction.end,
dist.from.start=dist.from.start)
# distribution of clip distances
pdf("Histogram_TTS_clipdistance.pdf", useDingbats = FALSE, width=6.4, height=5.4)
hist(bed.for.tts.eval$xy, xlab="clip distance (bp)",col="black",main="")
dev.off()

# identify gene ends
add.to.end = max(bed.for.tts.eval$xy)
knot.div = 40
pk.thresh = 0.05
bp.bin = 50
knot.thresh = 5
cnt.thresh = 5
tau.dist = 50000
frac.max = 1
frac.min = 0.3
inferred.coords = get.TTS(bed=bed.for.tts.eval, tss=TSS.gene.filtered3,
bw.plus=bw.plus, bw.minus=bw.minus,
bp.bin=bp.bin, add.to.end=add.to.end,
pk.thresh=pk.thresh, knot.thresh=knot.thresh,
cnt.thresh=cnt.thresh, tau.dist=tau.dist,
frac.max=frac.max, frac.min=frac.min, knot.div=knot.div)

coords = inferred.coords$bed
# check for the percentage of identified TTSs that match the search region boundry
TTS.boundary.match(coords=coords, bed.for.tts.eval=bed.for.tts.eval) #0.1968006
# look at whether TTSs identified at the search boundry had clipped boundries
frac.bound.clip = TTS.boundary.clip(coords=coords, bed.for.tts.eval=bed.for.tts.eval)
frac.bound.clip$frac.clip #0.919607
frac.bound.clip$frac.noclip #0.08039303
# plot didtribution of clip distances for genes with boundary TTSs
pdf("bp_clipped_for_boundry.pdf", useDingbats = FALSE, width=6.4, height=5.4)
hist(frac.bound.clip$clip.tts,main="",col="black",
xlab="bp clipped for boundry genes")
dev.off()

# plot new coord id
gene.end = bed.for.tts.eval
# ploting NR3C1
long.gene = largest.interval.bed
pdf("Genenr3c1.pdf", onefile=FALSE)
tts.plot(coords=coords, gene.end=gene.end, long.gene=long.gene,
gene = "Nr3c1", xper=0.1, yper=0.2,
bw.plus=bw.plus, bw.minus=bw.minus, bp.bin=5,
frac.min=frac.min, frac.max=frac.max,
add.to.end=add.to.end, tau.dist=tau.dist)
dev.off()

# plot tts curve
pdf("Gene_curvenr3c1.pdf", onefile=FALSE)
gene.end.plot(bed=bed.for.tts.eval, gene="Nr3c1",
bw.plus=bw.plus, bw.minus=bw.minus,
bp.bin=bp.bin, add.to.end=add.to.end, knot.div=knot.div,
pk.thresh=pk.thresh, knot.thresh=knot.thresh,
cnt.thresh=cnt.thresh, tau.dist=tau.dist,
frac.max=frac.max, frac.min=frac.min)
dev.off()

write.table(coords,'primary_transcript_annotation.bed',sep='\t',quote=F,col.names=F,row.names=F)
