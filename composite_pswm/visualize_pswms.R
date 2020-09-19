library(igraph)
library(dichromat)

setwd('~/Desktop')

threecol=read.table("3_col_combined_motif_db.txt",header=F,stringsAsFactors = F,sep='\t')
colnames(threecol)=c('from','to','e_value')
threecol$weight=abs(log(threecol$e_value, base = 10))



#create the graph variable
g=graph.data.frame(threecol,directed=F)
g=simplify(g)

cluster=clusters(g)

for(i in 1:length(groups(cluster))) {
write.table(groups(cluster)[i],file=paste0('PSWM_family_',i,'.txt'),col.names = F, row.names = F, quote = F, sep = '\t')
}

l=layout.fruchterman.reingold(g)
l=layout.norm(l,-1,1,-1,1)

colfunc<-colorRampPalette(c("red","yellow","springgreen","royalblue","purple"))
#pick a distinct color from the palette for each disease and save the list of colors as a vector
mycol = colfunc(length(groups(cluster)))

pdf(file='families.pdf',width=10,height=10)
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
	rown = y[,2]
	print(rown)
	y = as.data.frame(y[,4], stringsAsFactors = FALSE)
	rownames(y) = rown
	colnames(y) = i
	print(y)
	df = merge(df, y, by="row.names", all = TRUE)
	rownam = df[,1]
	df = as.data.frame(df[,2:ncol(df)])
	colnames(df) = vec.names
	rownames(df) = rownam
}



#makes any Inf values the maximum of the df
mx = max(threecol[,4][which(threecol[,4] < Inf)])
z = do.call(data.frame, lapply(df, function(x) replace(x, is.infinite(x),mx)))

rownames(z) = rownames(df)
colnames(z) = colnames(df)

                                       
#need to order based on conformity to consensus
#ord = names(cluster$membership)[order(cluster$membership)]

#doing it manually:
family1 = read.table('mike_generate_composites/PSWM_family_1_ranks.tomtom_output/tomtom.tsv', header = TRUE, skip = '#', sep = '\t', stringsAsFactors = FALSE)
family1 = family1[family1$Query_ID %in% names(cluster$membership[cluster$membership == 1]),]
ord1 = family1[order(family1$E.value),]$Query_ID

family2 = read.table('mike_generate_composites/PSWM_family_2_ranks.tomtom_output/tomtom.tsv', header = TRUE, skip = '#', sep = '\t', stringsAsFactors = FALSE)
family2 = family2[family2$Query_ID %in% names(cluster$membership[cluster$membership == 2]),]
ord2 = family2[order(family2$E.value),]$Query_ID

family3 = read.table('mike_generate_composites/PSWM_family_3_ranks.tomtom_output/tomtom.tsv', header = TRUE, skip = '#', sep = '\t', stringsAsFactors = FALSE)
family3 = family3[family3$Query_ID %in% names(cluster$membership[cluster$membership == 3]),]
ord3 = family3[order(family3$E.value),]$Query_ID

family4 = read.table('mike_generate_composites/PSWM_family_4_ranks.tomtom_output/tomtom.tsv', header = TRUE, skip = '#', sep = '\t', stringsAsFactors = FALSE)
family4 = family4[family4$Query_ID %in% names(cluster$membership[cluster$membership == 4]),]
ord4 = family4[order(family4$E.value),]$Query_ID

family5 = read.table('mike_generate_composites/PSWM_family_5_ranks.tomtom_output/tomtom.tsv', header = TRUE, skip = '#', sep = '\t', stringsAsFactors = FALSE)
family5 = family5[family5$Query_ID %in% names(cluster$membership[cluster$membership == 5]),]
ord5 = family5[order(family5$E.value),]$Query_ID

ord = c(ord1, ord2, ord3, ord4, ord5)



d = z[match(ord, rownames(z)), match(ord, colnames(z))]
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


#plot, the manipulate in AI:
pdf('pdf_test.pdf', width=10.4, height=9.5)
heatmap.2(as.matrix(avg.df),col=colorRampPalette(c("light blue", "blue", "dark blue"))(100),
          symbreaks=FALSE,scale="none",na.rm=TRUE, dendrogram = 'none', symm =TRUE, density.info=c("none"), key.title= '',trace=c("none"), Colv = FALSE, 
key.xlab = expression('-log'[10]*'(e-value)'), Rowv =FALSE, lhei=c(0.75,4), lwid = c(1.2, 4))
dev.off()



















