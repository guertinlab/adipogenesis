library(bigWig)
library(DESeq2)
library(DEGreport)
library(tibble)
library(lattice)
source('https://raw.githubusercontent.com/mjg54/znf143_pro_seq_analysis/master/docs/ZNF143_functions.R')

directory = '/Volumes/GUERTIN_2/adipogenesis/atac/'
setwd(directory)
preadipo.file = read.table("4hour/atac4h_0.05/3T3_atac_summit_200window.bed")

get.raw.counts.interval <- function(df, path.to.bigWig, file.prefix = 'H') {
    df = df[,1:5]
    vec.names = c()
    inten.df=data.frame(matrix(ncol = 0, nrow = nrow(df)))
    for (mod.bigWig in Sys.glob(file.path(path.to.bigWig,
                                          paste(file.prefix, "*.bigWig", sep ='')))) {
        factor.name = strsplit(strsplit(mod.bigWig, "/")[[1]][length(strsplit(mod.bigWig, "/")[[1]])], '\\.')[[1]][1]
        print(factor.name)
        vec.names = c(vec.names, factor.name)
        loaded.bw = load.bigWig(mod.bigWig)
        mod.inten = bed.region.bpQuery.bigWig(loaded.bw, df[,1:3])
        inten.df = cbind(inten.df, mod.inten)
    }
    colnames(inten.df) = vec.names
    r.names = paste(df[,1], ':', df[,2], '-', df[,3], sep='')
    row.names(inten.df) = r.names
    return(inten.df)
}

df.preadipo = get.raw.counts.interval(preadipo.file, directory, file.prefix = '3')
                                        #follow this: https://hbctraining.github.io/DGE_workshop/lessons/08_DGE_LRT.html
sample.conditions = factor(sapply(strsplit(as.character(colnames(df.preadipo)), '_'), '[', 2))

sample.conditions = factor(sample.conditions, levels=c("t0","20min","40min","60min","2hr","3hr","4hr"))

deseq.counts.table = DESeqDataSetFromMatrix(df.preadipo, as.data.frame(sample.conditions), ~ sample.conditions)

dds = DESeq(deseq.counts.table)

                                        #exploratory
rld = rlog(dds, blind=TRUE)

                                        #I think two ATAC-seq samples were labeled incorrectly.
#
x = plotPCA(rld, intgroup="sample.conditions", returnData=TRUE)
plotPCAlattice(x, file = 'PCA_adipogenesis_lattice_guertin.pdf')

                                        #lrt

dds.lrt = DESeq(dds, test="LRT", reduced = ~ 1)

#use for outside R
normalized.counts = counts(dds, normalized=TRUE)

res.lrt = results(dds_lrt)

padj.cutoff = 0.000001

siglrt.re = res.lrt[res.lrt$padj < padj.cutoff & !is.na(res.lrt$padj),]

dim(siglrt.re)                                    

rld_mat <- assay(rld)
cluster_rlog = rld_mat[rownames(siglrt.re),]
meta = as.data.frame(sample.conditions)
rownames(meta) = colnames(cluster_rlog)

save.image('191115_adipogenesis.Rdata')

#need to play around with this:
clusters <- degPatterns(cluster_rlog[1:1000,], metadata = meta, minc = 200, time = "sample.conditions", col=NULL, eachStep = TRUE)


class(clusters)
head(clusters$df)
cluster_groups <- clusters$df
group1 <- clusters$df %>%
          filter(cluster == 1)








#useful for network inference section 5.6 two wave network
                                        #https://kateto.net/netscix2016.html
#layered graph bipartite
#made a cyclical https://rdrr.io/cran/igraph/man/layout_with_sugiyama.html
#term sugiyama 
#https://stackoverflow.com/questions/51194653/r-is-not-taking-the-parameter-hgap-in-layout-with-sugiyama

sample.conditions = factor(c("untreated","untreated","untreated","untreated","treated","treated"), levels=c("t0","20min","40min","60min","2hr","3hr","4hr"))        
  deseq.counts.table = DESeqDataSetFromMatrix(mat, DataFrame(sample.conditions), ~ sample.conditions);
  colData(deseq.counts.table)$condition<-factor(colData(deseq.counts.table)$sample.conditions, levels=c('untreated','treated'));
dds = DESeq(deseq.counts.table);



conditions = sapply(names(reads0),function(x){strsplit(x,"_")[[1]][1]})
conditions = colnames(df.preadipo)
conditions = factor(conditions, levels=c("t0","20min","40min","60min","2hr","3hr","4hr"))


dds = DESeqDataSetFromMatrix(countData = df.preadipo, colData = meta, design = ~ sampletype)
#follow atac_time_cluster.R
# estimate size factors
size.factors.df.preadipo = estimateSizeFactorsForMatrix(df.preadipo)

## generate DESeqDataSet object 
conditions = sapply(names(df.preadipo),function(x){strsplit(x,"_")[[1]][1]})
colData = as.data.frame(conditions)

deseq.df.preadipo = DESeqDataSetFromMatrix(countData=df.preadipo, 
                                           colData=colData, design=~conditions)
sizeFactors(deseq.df.preadipo) = sizefactors.deseq.df.preadipo



# basic counts, size factors, and annotation
counts = counts(deseq.df.preadipo)
sizefac = sizeFactors(deseq.df.preadipo)
times = colData(deseq.df.preadipo)$conditions
t.order = as.data.frame(cbind(c("t0","20min","40min","60min","2hr","3hr","4hr"), 
                              c(0,0.33,0.67,1,2,3,4)),stringsAsFactors=FALSE)

names(t.order) = c("condition","time")

# specify new conditions as times (factor)
conditions = sapply(as.character(times),function(x){
  t.order$time[which(t.order$condition==x)]})
conditions = factor(conditions, levels=c(0,0.33,0.67,1,2,3,4))

# set new deseq object
colData = as.data.frame(conditions)
deseq_obj = DESeqDataSetFromMatrix(countData=counts, 
                  colData=colData, design=~conditions)
sizeFactors(deseq_obj) = estimateSizeFactorsForMatrix(counts)

# remove the 6day data
counts_preadip = as.data.frame(counts(deseq_obj))
#counts_preadip = counts[,-grep("6d",names(counts))]
sizefac_preadip = sizeFactors(counts_preadip)
conditions_preadip = colData(counts_preadip)

# set new deseq object
colData_preadip = as.data.frame(conditions_preadip)
deseq_obj_preadip = DESeqDataSetFromMatrix(countData=counts_preadip, 
                                           colData=colData_preadip, 
                                           design=~conditions_preadip)
sizeFactors(deseq_obj_preadip) = sizefac_preadip



############################################################
## implement differential peak analysis
############################################################

# differential expreassion analysis - temporal dynamics based on LRT
deg_preadip = DESeq(deseq_obj_preadip, test="LRT", 
                    full=~conditions_preadip,  reduced=~1)
res_preadip = results(deg_preadip)

############################################################
## isolate dynamic peaks for meme
############################################################

# set significant and unsignificant peaks
ind.sig = which(res_preadip$padj<sig.thresh)
ind.un = which(res_preadip$padj>sig.un)
length(ind.sig)
length(ind.un)

# function to find peak info from rownames of the results frame
get.map.from.res = function(res,map){
  namen = rownames(res)
  chrs = sapply(namen,function(x){strsplit(x,":")[[1]][1]})
  ends = sapply(namen,function(x){strsplit(x,"-")[[1]][2]})
  strs = sapply(namen,function(x){y=strsplit(x,"-")[[1]][1]; 
    return(strsplit(y,":")[[1]][2]) })
  coords = cbind(chrs,strs,ends) %>% as.data.frame(stringsAsFactors=F)
  out = row.match(coords, map[,1:3])
  return( cbind(res, map[out,]) )
} # get.map.from.res

# get data for sig and unsig peaks
sig.pks = get.map.from.res(res=res_preadip[ind.sig,], map=bed0)
uns.pks = get.map.from.res(res=res_preadip[ind.un,], map=bed0)

# set coordinates around each summit based on pk.dist
pk.coords.for.meme = function(res=NULL, pk.dist=NULL){
  out = c()
  for(ii in 1:nrow(res)){
    pks = res$summits[ii]
    for(jj in 1:length(pks)){
      new = c(res$chr[ii], pks[jj]-pk.dist, pks[jj]+pk.dist, rownames(res)[ii],
              signif(res$log2FoldChange[ii],3), signif(res$padj[ii],3) )
      out = rbind(out,new) 
    } # jj
  } # ii
  out = out %>% as.data.frame(stringsAsFactors=FALSE)
  names(out) = c("chr","start","end","peakID","log2fc","fdr")
  rownames(out) = c(1:nrow(out))
  return(out)
} # pk.coords.for.meme

out.sig = pk.coords.for.meme(res=sig.pks, pk.dist=pk.dist)
out.uns = pk.coords.for.meme(res=uns.pks, pk.dist=pk.dist)

# save insignificant peaks for enrichment analysis
# see Motif_in_peak_diffDyn.R and enrichSummary_diffDyn.sh
save(out.uns, file="out.uns.RData")
















############################################################
## normalize raw counts to size factors, generate DESeq object
############################################################

# estimate size factors
size_factors = estimateSizeFactorsForMatrix(reads0)

## generate DESeqDataSet object and transform count data to log2 scale
## note. enter size factors before implementing log transform
conditions = sapply(names(reads0),function(x){strsplit(x,"_")[[1]][1]})
colData = as.data.frame(conditions)
deseq_obj_sizefac = DESeqDataSetFromMatrix(countData=reads0, 
                                           colData=colData, design=~conditions)
sizeFactors(deseq_obj_sizefac) = size_factors


############################################################
## adjust data organization
## exclude 6day data
############################################################

# basic counts, size factors, and annotation
counts = counts(deseq_obj_sizefac)
sizefac = sizeFactors(deseq_obj_sizefac)
times = colData(deseq_obj_sizefac)$conditions
t.order = cbind(c("t0","20min","40min","60min","2hr","3hr","4hr","6d"),
                c(0,0.33,0.67,1,2,3,4,144)) %>% 
  as.data.frame(stringsAsFactors=FALSE)
names(t.order) = c("condition","time")

# specify new conditions as times (factor)
conditions = sapply(as.character(times),function(x){
  t.order$time[which(t.order$condition==x)]})
conditions = factor(conditions, levels=c(0,0.33,0.67,1,2,3,4,144))

# set new deseq object
colData = as.data.frame(conditions)
deseq_obj = DESeqDataSetFromMatrix(countData=counts, 
                  colData=colData, design=~conditions)
sizeFactors(deseq_obj) = estimateSizeFactorsForMatrix(counts)

# remove the 6day data
counts = counts(deseq_obj) %>% as.data.frame
counts_preadip = counts[,-grep("6d",names(counts))]
sizefac_preadip = sizeFactors(deseq_obj)[-grep("6d",names(counts))]
conditions_preadip = colData(deseq_obj)$conditions[-grep("6d",names(counts))]

# set new deseq object
colData_preadip = as.data.frame(conditions_preadip)
deseq_obj_preadip = DESeqDataSetFromMatrix(countData=counts_preadip, 
                                           colData=colData_preadip, 
                                           design=~conditions_preadip)
sizeFactors(deseq_obj_preadip) = sizefac_preadip

############################################################
## implement differential peak analysis
############################################################

# differential expreassion analysis - temporal dynamics based on LRT
deg_preadip = DESeq(deseq_obj_preadip, test="LRT", 
                    full=~conditions_preadip,  reduced=~1)
res_preadip = results(deg_preadip)

############################################################
## isolate dynamic peaks for meme
############################################################

# set significant and unsignificant peaks
ind.sig = which(res_preadip$padj<sig.thresh)
ind.un = which(res_preadip$padj>sig.un)
length(ind.sig)
length(ind.un)

# function to find peak info from rownames of the results frame
get.map.from.res = function(res,map){
  namen = rownames(res)
  chrs = sapply(namen,function(x){strsplit(x,":")[[1]][1]})
  ends = sapply(namen,function(x){strsplit(x,"-")[[1]][2]})
  strs = sapply(namen,function(x){y=strsplit(x,"-")[[1]][1]; 
    return(strsplit(y,":")[[1]][2]) })
  coords = cbind(chrs,strs,ends) %>% as.data.frame(stringsAsFactors=F)
  out = row.match(coords, map[,1:3])
  return( cbind(res, map[out,]) )
} # get.map.from.res

# get data for sig and unsig peaks
sig.pks = get.map.from.res(res=res_preadip[ind.sig,], map=bed0)
uns.pks = get.map.from.res(res=res_preadip[ind.un,], map=bed0)

# set coordinates around each summit based on pk.dist
pk.coords.for.meme = function(res=NULL, pk.dist=NULL){
  out = c()
  for(ii in 1:nrow(res)){
    pks = res$summits[ii]
    for(jj in 1:length(pks)){
      new = c(res$chr[ii], pks[jj]-pk.dist, pks[jj]+pk.dist, rownames(res)[ii],
              signif(res$log2FoldChange[ii],3), signif(res$padj[ii],3) )
      out = rbind(out,new) 
    } # jj
  } # ii
  out = out %>% as.data.frame(stringsAsFactors=FALSE)
  names(out) = c("chr","start","end","peakID","log2fc","fdr")
  rownames(out) = c(1:nrow(out))
  return(out)
} # pk.coords.for.meme

out.sig = pk.coords.for.meme(res=sig.pks, pk.dist=pk.dist)
out.uns = pk.coords.for.meme(res=uns.pks, pk.dist=pk.dist)

# save insignificant peaks for enrichment analysis
# see Motif_in_peak_diffDyn.R and enrichSummary_diffDyn.sh
save(out.uns, file="out.uns.RData")









































#each column of the count matrix can be divided by the respective size factor to get a matrix of relative accessibility per replicate/condition
#check to see if df.preadipo is updated (or changed) when this is run
head(sapply(c(1:24),function(x){df.preadipo[,x]/sizefactors.deseq.df.preadipo[x]})[,1:5])
head(counts(deseq.df.preadipo, normalized=TRUE)[,1:5])


# basic counts, size factors, and annotation
counts = counts(deseq.df.preadipo)
sizefac = sizeFactors(deseq.df.preadipo)
times = colData(deseq.df.preadipo)$conditions
t.order = cbind(c("t0","20min","40min","60min","2hr","3hr","4hr","6d"),
                c(0,0.33,0.67,1,2,3,4,144)) %>% 
  as.data.frame(stringsAsFactors=FALSE)
names(t.order) = c("condition","time")

# specify new conditions as times (factor)
conditions = sapply(as.character(times),function(x){
  t.order$time[which(t.order$condition==x)]})
conditions = factor(conditions, levels=c(0,0.33,0.67,1,2,3,4,144))

# set new deseq object
colData = as.data.frame(conditions)
deseq_obj = DESeqDataSetFromMatrix(countData=counts, 
                                   colData=colData, design=~conditions)
sizeFactors(deseq_obj) = estimateSizeFactorsForMatrix(counts)


############################################################
## exclude 6day data
############################################################

# remove the 6day data
counts = counts(deseq_obj) %>% as.data.frame
counts_preadip = counts[,-grep("6d",names(counts))]
sizefac_preadip = sizeFactors(deseq_obj)[-grep("6d",names(counts))]
conditions_preadip = colData(deseq_obj)$conditions[-grep("6d",names(counts))]
condit.map = cbind(as.character(conditions_preadip), names(counts_preadip)) %>% 
  as.data.frame(stringsAsFactors=F)
names(condit.map) = c("time","rep")

# set new deseq object
colData_preadip = as.data.frame(conditions_preadip)
deseq_obj_preadip = DESeqDataSetFromMatrix(countData=counts_preadip, 
                                           colData=colData_preadip, 
                                           design=~conditions_preadip)
sizeFactors(deseq_obj_preadip) = sizefac_preadip

# save(deseq_obj_preadip, file="atacDeseq.RData")

############################################################
## generate all pairwise comparisons
############################################################

# loop through all pairwise combinations and run DESeq 
condits = c(0,0.33,0.67,1,2,3,4) %>% as.character
res.pairs = list()
for(ii in 1:(length(condits)-1)){
  for(jj in (ii+1):length(condits)){
    
    print(paste0("compare ",condits[ii]," to ",condits[jj]))
    
    # basic data annotation
    ind_ii = which(conditions_preadip == condits[ii])
    ind_jj = which(conditions_preadip == condits[jj])
    
    # set new deseq object
    des = conditions_preadip[c(ind_ii,ind_jj)]
    colData_ij = as.data.frame(des)
    cnt_ij = counts_preadip[,c(ind_ii,ind_jj)]
    deseq_obj_preadip = DESeqDataSetFromMatrix(countData=cnt_ij, 
                                               colData=colData_ij, 
                                               design=~des)
    sizeFactors(deseq_obj_preadip) = sizefac_preadip[c(ind_ii,ind_jj)]
    colData(deseq_obj_preadip)$condition = des
    
    # pairwise DEG analysis
    dds = DESeq(deseq_obj_preadip)
    res = results(dds)[order(results(dds)$padj),]
    res.pairs[[paste0(condits[ii],"_",condits[jj])]] = res
    
  } # jj
} # ii

# save(res.pairs, file="pairwise.deg.RData")
# load("pairwise.deg.RData") 






















library(dplyr)
library(DESeq2)
library(ggplot2)
library(bigWig)
library(prodlim)

source("atac_norm_functions.R")

dir = paste0("/media/wa3j/Seagate2/Documents/PRO/",
             "adipogenesis/July2018/atac_time_deg")
setwd(dir)


############################################################
## key analysis parameters
############################################################

# set number of base pairs around peaks for motif analysis
pk.dist = 50

# parameters for designating dynamic peaks
sig.thresh = 0.001
fc.thresh = 1

# parameters for designating non-dynamic peaks
sig.un = 0.5
fc.un = 0.25


############################################################
## import macs2 peak data
############################################################

# peak and summit coordinates (summits.bed, from empiricalSummits.R)
load("ATACsummits_20191013.RData")
bed0 = summits.bed

# chrom sizes
chrm.size = read.table("mm10.chrom.sizes",header=F,stringsAsFactors=F,sep="\t")
chrm.size = cbind(chrm.size[,1], 1, chrm.size[,2])
chrm.size = as.data.frame(chrm.size, stringsAsFactors=F)
names(chrm.size) = c("chr","start","end")
chrm.size[,2:3] = apply(chrm.size[,2:3],2,function(x){data.matrix(x) %>% as.numeric})

# write peaks
# write.table(bed0[,1:3],"atacPeaks.bed0",col.names=F,row.names=F,sep="\t",quote=F)

############################################################
## process peak/summit data 
############################################################

# filter chromosomes
unique(bed0$chr)
dim(bed0) # 157049
chr.keep = paste0("chr",c(1:19))
bed0 = bed0[bed0$chr %in% chr.keep,]
dim(bed0) # 152729

# look at peak distances
d = bed0$end - bed0$start
hist(d,col="black",xlab="atac peak dist (bp)")
median(d)
quantile(d)

# save coords
save(bed0, file="bed.map20191014.RData")


############################################################
## import atac-seq data and set directory for data output
############################################################

# atac data mapped to peak regions
bigWig.path = paste0("/media/wa3j/Seagate2/Documents/PRO/adipogenesis",
                     "/July2018/atac_bw_20191014/atac_all")
file.prefix = "3T3_"
file.suffix = ".bigWig"
min.length = 1
reads0 = get.counts.interval(bed=bed0[,1:3], bigWig.path=bigWig.path, 
                             file.prefix=file.prefix,
                             file.suffix=file.suffix,
                             min.length=min.length)


############################################################
## normalize raw counts to size factors, generate DESeq object
############################################################

# estimate size factors
size_factors = estimateSizeFactorsForMatrix(reads0)

## generate DESeqDataSet object and transform count data to log2 scale
conditions = sapply(names(reads0),function(x){strsplit(x,"_")[[1]][1]})
colData = as.data.frame(conditions)
deseq_obj_sizefac = DESeqDataSetFromMatrix(countData=reads0, 
                                           colData=colData, design=~conditions)
sizeFactors(deseq_obj_sizefac) = size_factors

# note that each column of the count matrix was divided 
# by the respective size factor
head(sapply(c(1:24),function(x){reads0[,x]/size_factors[x]})[,1:5])
head(counts(deseq_obj_sizefac, normalized=TRUE)[,1:5])

############################################################
## adjust data organization
############################################################

# basic counts, size factors, and annotation
counts = counts(deseq_obj_sizefac)
sizefac = sizeFactors(deseq_obj_sizefac)
times = colData(deseq_obj_sizefac)$conditions
t.order = cbind(c("t0","20min","40min","60min","2hr","3hr","4hr","6d"),
                c(0,0.33,0.67,1,2,3,4,144)) %>% 
  as.data.frame(stringsAsFactors=FALSE)
names(t.order) = c("condition","time")

# specify new conditions as times (factor)
conditions = sapply(as.character(times),function(x){
  t.order$time[which(t.order$condition==x)]})
conditions = factor(conditions, levels=c(0,0.33,0.67,1,2,3,4,144))

# set new deseq object
colData = as.data.frame(conditions)
deseq_obj = DESeqDataSetFromMatrix(countData=counts, 
                                   colData=colData, design=~conditions)
sizeFactors(deseq_obj) = estimateSizeFactorsForMatrix(counts)


############################################################
## exclude 6day data
############################################################

# remove the 6day data
counts = counts(deseq_obj) %>% as.data.frame
counts_preadip = counts[,-grep("6d",names(counts))]
sizefac_preadip = sizeFactors(deseq_obj)[-grep("6d",names(counts))]
conditions_preadip = colData(deseq_obj)$conditions[-grep("6d",names(counts))]
condit.map = cbind(as.character(conditions_preadip), names(counts_preadip)) %>% 
  as.data.frame(stringsAsFactors=F)
names(condit.map) = c("time","rep")

# set new deseq object
colData_preadip = as.data.frame(conditions_preadip)
deseq_obj_preadip = DESeqDataSetFromMatrix(countData=counts_preadip, 
                                           colData=colData_preadip, 
                                           design=~conditions_preadip)
sizeFactors(deseq_obj_preadip) = sizefac_preadip

# save(deseq_obj_preadip, file="atacDeseq.RData")

############################################################
## generate all pairwise comparisons
############################################################

# loop through all pairwise combinations and run DESeq 
condits = c(0,0.33,0.67,1,2,3,4) %>% as.character
res.pairs = list()
for(ii in 1:(length(condits)-1)){
  for(jj in (ii+1):length(condits)){
    
    print(paste0("compare ",condits[ii]," to ",condits[jj]))
    
    # basic data annotation
    ind_ii = which(conditions_preadip == condits[ii])
    ind_jj = which(conditions_preadip == condits[jj])
    
    # set new deseq object
    des = conditions_preadip[c(ind_ii,ind_jj)]
    colData_ij = as.data.frame(des)
    cnt_ij = counts_preadip[,c(ind_ii,ind_jj)]
    deseq_obj_preadip = DESeqDataSetFromMatrix(countData=cnt_ij, 
                                               colData=colData_ij, 
                                               design=~des)
    sizeFactors(deseq_obj_preadip) = sizefac_preadip[c(ind_ii,ind_jj)]
    colData(deseq_obj_preadip)$condition = des
    
    # pairwise DEG analysis
    dds = DESeq(deseq_obj_preadip)
    res = results(dds)[order(results(dds)$padj),]
    res.pairs[[paste0(condits[ii],"_",condits[jj])]] = res
    
  } # jj
} # ii

# save(res.pairs, file="pairwise.deg.RData")
# load("pairwise.deg.RData") 

############################################################
## generate output for meme
############################################################

# function to find peak info from rownames of the results frame
get.map.from.res = function(res,map){
  namen = rownames(res)
  chrs = sapply(namen,function(x){strsplit(x,":")[[1]][1]})
  ends = sapply(namen,function(x){strsplit(x,"-")[[1]][2]})
  strs = sapply(namen,function(x){y=strsplit(x,"-")[[1]][1]; 
        return(strsplit(y,":")[[1]][2]) })
  coords = cbind(chrs,strs,ends) %>% as.data.frame(stringsAsFactors=F)
  out = row.match(coords, map[,1:3])
  return( cbind(res, map[out,]) )
} # get.map.from.res

# set coordinates around each summit based on pk.dist
# generate output format for meme
pk.coords.for.meme = function(res=NULL, pk.dist=NULL){
  out = c()
  for(ii in 1:nrow(res)){
    pks = res$summits[ii]
    new = c(res$chr[ii], pks-pk.dist, pks+pk.dist, rownames(res)[ii],
              signif(res$log2FoldChange[ii],3), "+", signif(res$padj[ii],3) )
    out = rbind(out,new) 
  } # ii
  out = out %>% as.data.frame(stringsAsFactors=FALSE)
  names(out) = c("chr","start","end","peak_name","log2fc","str","fdr")
  rownames(out) = c(1:nrow(out))
  return(out)
} # pk.coords.for.meme

# loop through all pairwise data and write output for meme
nosig = c()
for(ii in 1:length(res.pairs)){
  
  # basic annotation
  comp_ii = names(res.pairs)[[ii]]
  res_ii = res.pairs[[ii]]
  indsig = which(res_ii$padj<sig.thresh & abs(res_ii$log2FoldChange)>fc.thresh)
  ind.up = which(res_ii$padj<sig.thresh & res_ii$log2FoldChange>fc.thresh)
  ind.dn = which(res_ii$padj<sig.thresh & res_ii$log2FoldChange<(-1)*fc.thresh)
  ind.un = which(res_ii$padj>sig.un & abs(res_ii$log2FoldChange)<fc.un)
  
  # if there are more insignificant peaks, use the number of sig peaks
  # select those with the highest FDRs
  if(length(indsig)==0){nosig = c(nosig, comp_ii); next}
  if(length(ind.un) > length(indsig)){ind.un = rev(ind.un)[1:length(indsig)]}
  
  # get macs2 data for sig and unsig peaks
  sig.pks.up = get.map.from.res(res=res_ii[ind.up,], map=bed0)
  sig.pks.dn = get.map.from.res(res=res_ii[ind.dn,], map=bed0)
  uns.pks = get.map.from.res(res=res_ii[ind.un,], map=bed0)
  out.sig.up = pk.coords.for.meme(res=sig.pks.up, pk.dist=pk.dist)
  out.sig.dn = pk.coords.for.meme(res=sig.pks.dn, pk.dist=pk.dist)
  out.uns = pk.coords.for.meme(res=uns.pks, pk.dist=pk.dist)
  
  # set directory and output data for meme analysis
  dir_ii = paste0("meme_",comp_ii)
  system(paste0("mkdir ",dir_ii))
  setwd(dir_ii)
  fname = paste0("upsig_",comp_ii,".bed")
  write.table(out.sig.up,fname,sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE)
  fname = paste0("downsig_",comp_ii,".bed")
  write.table(out.sig.dn,fname,sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE)
  fname = paste0("unsig_",comp_ii,".bed")
  write.table(out.uns,fname,sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE)
  setwd(dir)
  
} # ii










library(dplyr)
library(DESeq2)
library(ggplot2)
library(bigWig)
library(prodlim)

source("atac_norm_functions.R")

dir = paste0("/media/wa3j/Seagate2/Documents/PRO/",
             "adipogenesis/July2018/atac_time_deg")
setwd(dir)


############################################################
## key analysis parameters
############################################################

# set number of base pairs around peaks for motif analysis
pk.dist = 50

# parameters for designating dynamic peaks
sig.thresh = 0.001

# parameters for designating non-dynamic peaks
sig.un = 0.5


############################################################
## import macs2 peak data
############################################################

# load bed0 - filtered ATAC peak coords, see atac_time_deg_meme.R
load("bed.map20191014.RData")

############################################################
## import atac-seq data and set directory for data output
############################################################

# atac data mapped to peak regions
bigWig.path = paste0("/media/wa3j/Seagate2/Documents/PRO/adipogenesis",
                     "/July2018/atac_bw_20191014/atac_all")
file.prefix = "3T3_"
file.suffix = ".bigWig"
min.length = 1
reads0 = get.counts.interval(bed=bed0[,1:3], bigWig.path=bigWig.path, 
                             file.prefix=file.prefix,
                             file.suffix=file.suffix,
                             min.length=min.length)


############################################################
## normalize raw counts to size factors, generate DESeq object
############################################################

# estimate size factors
size_factors = estimateSizeFactorsForMatrix(reads0)

## generate DESeqDataSet object and transform count data to log2 scale
## note. enter size factors before implementing log transform
conditions = sapply(names(reads0),function(x){strsplit(x,"_")[[1]][1]})
colData = as.data.frame(conditions)
deseq_obj_sizefac = DESeqDataSetFromMatrix(countData=reads0, 
                                           colData=colData, design=~conditions)
sizeFactors(deseq_obj_sizefac) = size_factors


############################################################
## adjust data organization
## exclude 6day data
############################################################

# basic counts, size factors, and annotation
counts = counts(deseq_obj_sizefac)
sizefac = sizeFactors(deseq_obj_sizefac)
times = colData(deseq_obj_sizefac)$conditions
t.order = cbind(c("t0","20min","40min","60min","2hr","3hr","4hr","6d"),
                c(0,0.33,0.67,1,2,3,4,144)) %>% 
  as.data.frame(stringsAsFactors=FALSE)
names(t.order) = c("condition","time")

# specify new conditions as times (factor)
conditions = sapply(as.character(times),function(x){
  t.order$time[which(t.order$condition==x)]})
conditions = factor(conditions, levels=c(0,0.33,0.67,1,2,3,4,144))

# set new deseq object
colData = as.data.frame(conditions)
deseq_obj = DESeqDataSetFromMatrix(countData=counts, 
                  colData=colData, design=~conditions)
sizeFactors(deseq_obj) = estimateSizeFactorsForMatrix(counts)

# remove the 6day data
counts = counts(deseq_obj) %>% as.data.frame
counts_preadip = counts[,-grep("6d",names(counts))]
sizefac_preadip = sizeFactors(deseq_obj)[-grep("6d",names(counts))]
conditions_preadip = colData(deseq_obj)$conditions[-grep("6d",names(counts))]

# set new deseq object
colData_preadip = as.data.frame(conditions_preadip)
deseq_obj_preadip = DESeqDataSetFromMatrix(countData=counts_preadip, 
                                           colData=colData_preadip, 
                                           design=~conditions_preadip)
sizeFactors(deseq_obj_preadip) = sizefac_preadip

############################################################
## implement differential peak analysis
############################################################

# differential expreassion analysis - temporal dynamics based on LRT
deg_preadip = DESeq(deseq_obj_preadip, test="LRT", 
                    full=~conditions_preadip,  reduced=~1)
res_preadip = results(deg_preadip)

############################################################
## isolate dynamic peaks for meme
############################################################

# set significant and unsignificant peaks
ind.sig = which(res_preadip$padj<sig.thresh)
ind.un = which(res_preadip$padj>sig.un)
length(ind.sig)
length(ind.un)

# function to find peak info from rownames of the results frame
get.map.from.res = function(res,map){
  namen = rownames(res)
  chrs = sapply(namen,function(x){strsplit(x,":")[[1]][1]})
  ends = sapply(namen,function(x){strsplit(x,"-")[[1]][2]})
  strs = sapply(namen,function(x){y=strsplit(x,"-")[[1]][1]; 
    return(strsplit(y,":")[[1]][2]) })
  coords = cbind(chrs,strs,ends) %>% as.data.frame(stringsAsFactors=F)
  out = row.match(coords, map[,1:3])
  return( cbind(res, map[out,]) )
} # get.map.from.res

# get data for sig and unsig peaks
sig.pks = get.map.from.res(res=res_preadip[ind.sig,], map=bed0)
uns.pks = get.map.from.res(res=res_preadip[ind.un,], map=bed0)

# set coordinates around each summit based on pk.dist
pk.coords.for.meme = function(res=NULL, pk.dist=NULL){
  out = c()
  for(ii in 1:nrow(res)){
    pks = res$summits[ii]
    for(jj in 1:length(pks)){
      new = c(res$chr[ii], pks[jj]-pk.dist, pks[jj]+pk.dist, rownames(res)[ii],
              signif(res$log2FoldChange[ii],3), signif(res$padj[ii],3) )
      out = rbind(out,new) 
    } # jj
  } # ii
  out = out %>% as.data.frame(stringsAsFactors=FALSE)
  names(out) = c("chr","start","end","peakID","log2fc","fdr")
  rownames(out) = c(1:nrow(out))
  return(out)
} # pk.coords.for.meme

out.sig = pk.coords.for.meme(res=sig.pks, pk.dist=pk.dist)
out.uns = pk.coords.for.meme(res=uns.pks, pk.dist=pk.dist)

# save insignificant peaks for enrichment analysis
# see Motif_in_peak_diffDyn.R and enrichSummary_diffDyn.sh
save(out.uns, file="out.uns.RData")

############################################################
## data transformations for STEM
## /media/wa3j/Seagate2/Documents/software/stem
## java -mx1024M -jar /media/wa3j/Seagate2/Documents/software/stem/stem.jar
############################################################

# significance filter
ind.sig = sapply(out.sig$peakID,function(x){which(rownames(deseq_obj_preadip)==x)})

# log transformation
log2_exp = rlogTransformation(deseq_obj_preadip[ind.sig,], blind=TRUE)
log.pk.cnt = assay(log2_exp)

# z-score the data for each peak
pk.z.scores = scale(t(log.pk.cnt)) %>% t

# get mean z-scores for each time point
times = c("t0","20min","40min","60min","2h","3h","4h")
t.profiles = sapply(times,function(x){ 
  apply(pk.z.scores[,grep(x,colnames(pk.z.scores))],1,mean) 
})
colnames(t.profiles) = c(0,0.33,0.67,1,2,3,4)

# revise names for timeseries clustering
gene_symbol = rownames(t.profiles)
out = cbind(gene_symbol, t.profiles)
write.table(out,"pkdata_20191014.txt",sep="\t",quote=F,col.names=T,row.names=F)

############################################################
## process/graph STEM clustering results
############################################################

# function to capitalize first letter
simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
        sep="", collapse=" ")
}

# import STEM results, adjust gene names
stem.res0 = read.table("time_cluster_results_20191014.txt",sep="\t",
                       header=T,fill=T,stringsAsFactors=F,skip=1)
stem.res = stem.res0 %>% select(gene_symbol, Profile)
stem.res$gene_symbol = tolower(stem.res$gene_symbol)

# write out coordinate information for enrichment analysis
# see Motif_in_peak_diffDyn.R and enrichSummary_diffDyn.sh
save(stem.res, file="stem.res.RData")

# threshold number of profiles for plotting
n.profiles = 0

# generate plots
library(graphics)
tt = c(0,0.33,0.67,1,2,3,4)
pdf("ATACpeak_timeseries_clusters.pdf"); par(mfrow=c(3,3))
for(ii in 1:length(unique(stem.res$Profile))){
  pr = unique(stem.res$Profile)[ii]
  pk.pr = stem.res$gene_symbol[which(stem.res$Profile==pr)]
  if(length(pk.pr) < n.profiles){next}
  ex.pr = t.profiles[rownames(t.profiles) %in% pk.pr,]
  av.pr = apply(ex.pr,2,mean)
  tp = rep("l",nrow(ex.pr))
  ly = rep(1,nrow(ex.pr))
  matplot(tt,t(ex.pr),type=tp,lty=ly,col="black",
          ylab=paste0("profile ",pr),main="",xlab="time (hrs)")
  lines(tt,av.pr,col="red",cex=4)
}
dev.off()


############################################################
## output profiles of interest for meme analysis
############################################################

# clusters
profiles = unique(stem.res$Profile)

# loop through profiles and write data for meme
for(ii in profiles){
  pk.pr = stem.res$gene_symbol[which(stem.res$Profile==ii)] %>%
    unlist
  pk.un = out.uns$peakID[sample.int( nrow(out.uns), min(nrow(out.uns), length(pk.pr)) )]
  fname = paste0("preadip_sigDynamics_profile",ii,".bed")
  out = out.sig[out.sig$peakID %in% pk.pr,]
  write.table(out,fname,sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE)
  fname = paste0("preadip_unDynamics_profile",ii,".bed")
  out = out.uns[out.uns$peakID %in% pk.un,]
  write.table(out,fname,sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE)
} ## ii, loop profiles



