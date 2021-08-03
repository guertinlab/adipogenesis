library(igraph)
library(dichromat)

setwd('~/Desktop/fasta')

threecol=read.table("test_output.txt",header=F,stringsAsFactors = F,sep='\t')
threecol = threecol[,c(1,3, 13, 14)]
colnames(threecol)=c('from','to','eval','sw')

                                        #filter based on evalue
threecol= threecol[threecol$eval < 1e-5,]
threecol$eval=abs(log(threecol$eval, base = 10))



threecol$to = sapply(strsplit(sapply(strsplit(as.character(threecol$to), "Mus_musculus_"), "[[", 2),"-annot"), "[[", 1)
threecol$from = sapply(strsplit(sapply(strsplit(as.character(threecol$from), "Mus_musculus_"), "[[", 2),"-annot"), "[[", 1)

threecol$to = sapply(strsplit(sapply(strsplit(as.character(threecol$to), "annot"), "[[", 1),"_"), "[[", 1)
threecol$from = sapply(strsplit(sapply(strsplit(as.character(threecol$from), "annot"), "[[", 1),"_"), "[[", 1)


threecol = aggregate(eval~from + to, data=threecol, FUN=max)


#create the graph variable
g=graph.data.frame(threecol,directed=F)
g=simplify(g)

cluster=clusters(g)

for(i in 1:length(groups(cluster))) {
write.table(groups(cluster)[i],file=paste0('FASTA_family_',i,'.txt'),col.names = F, row.names = F, quote = F, sep = '\t')
}

l=layout.fruchterman.reingold(g)
l=layout.norm(l,-1,1,-1,1)

colfunc<-colorRampPalette(c("red","yellow","springgreen","royalblue","purple"))
#pick a distinct color from the palette for each disease and save the list of colors as a vector
mycol = colfunc(length(groups(cluster)))

pdf(file='dbd_families.pdf',width=10,height=10)
plot(g,layout=l,rescale=F,vertex.label.cex=.5,xlim=range(l[,1]),  ylim=range(l[,2]),
edge.width=E(g)$weight/20,vertex.size=degree(g,mode='out')/5,
edge.curved=T,vertex.label=NA,vertex.color=mycol[cluster$membership],
margin=0,asp=0)
dev.off()




df = as.data.frame(matrix(ncol=0, nrow=0),stringsAsFactors = FALSE)

vec.names = c()
for (i in unique(threecol[,1])) {
	vec.names = c(vec.names,i)
	y = threecol[threecol[,1] == i,]
        #print(y)
	rown = y[,2]
	#print(rown)
	y = as.data.frame(y[,3], stringsAsFactors = FALSE)
	rownames(y) = rown
	colnames(y) = i
	#print(y)
	df = merge(df, y, by="row.names", all = TRUE)
	rownam = df[,1]
	df = as.data.frame(df[,2:ncol(df)])
	colnames(df) = vec.names
	rownames(df) = rownam
}



#makes any Inf values the maximum of the df
#mx = max(threecol[,3][which(threecol[,3] < Inf)])
#z = do.call(data.frame, lapply(df, function(x) replace(x, is.infinite(x),mx)))

#rownames(z) = rownames(df)
#colnames(z) = colnames(df)

   
order.class = read.table('order_file.txt')
order.class$V1 = sapply(strsplit(sapply(strsplit(as.character(order.class$V1), "Mus_musculus_"), "[[", 2),"-annot"), "[[", 1)
order.class$V1 = sapply(strsplit(sapply(strsplit(as.character(order.class$V1), "annot"), "[[", 1),"_"), "[[", 1)
order.class$V1 %in% rownames(df)

order.class.vec = unique(order.class$V1)
d = df[match(order.class.vec, rownames(df)), match(order.class.vec, colnames(df))]



#d = z[rownames(z), colnames(z)]
###d = z[match(ord, rownames(z)), match(ord, colnames(z))]
#d = z[match(colnames(z), rownames(z)),]

td <- t(d)

#this makes any NA values with diagonal mirrored NAs convert to the
                                        #diagonal mirror value (reciprocal comparisons have differnt values)

for (i in 1:nrow(d)) {
    for (j in 1:ncol(d)) {
        if (is.na(d[i,j]) & !is.na(d[j,i])) {
            d[i,j] = d[j,i]
            }}}
for (i in 1:nrow(td)) {
    for (j in 1:ncol(td)) {
        if (is.na(td[i,j]) & !is.na(td[j,i])) {
            td[i,j] = td[j,i]
            }}}


#reciprocal comparisons have slightly differnt e-values, average them:
avg.df = (d + td)/2

#only one triangle
avg.df[upper.tri(avg.df, diag = FALSE)] <- NA 



library(gplots)
#plot, the manipulate in AI:
pdf('pdf_test.pdf', width=10.4, height=9.5)
heatmap.2(as.matrix(avg.df),col=colorRampPalette(c("light blue", "blue", "dark blue"))(100),
          symbreaks=FALSE,scale="none",na.rm=TRUE, dendrogram = 'none', symm =TRUE, density.info=c("none"), key.title= '',trace=c("none"), Colv = FALSE, 
key.xlab = expression('-log'[10]*'(e-value)'), Rowv =FALSE, lhei=c(0.75,4), lwid = c(1.2, 4))
dev.off()















