Args=commandArgs(TRUE)
arg.dir = Args[1]

library(igraph)

setwd(arg.dir)

threecol=read.csv(paste0(arg.dir,"/3_col_combined_motif_db.txt"),
header=F,stringsAsFactors = F,sep='\t')
colnames(threecol)=c('from','to','e_value')
threecol$weight=abs(log(threecol$e_value))

#create the graph variable
g=graph.data.frame(threecol,directed=F)
g=simplify(g)

cluster=clusters(g)

for(i in 1:length(groups(cluster))) {
write.table(groups(cluster)[i],file=paste0('PSWM_family_',i,'.txt'),col.names = F, row.names = F, quote = F, sep = '\t')
}

l=layout.fruchterman.reingold(g)
l=layout.norm(l,-1,1,-1,1)

#note: please change the number of colors depending on how many of communities you got
#spare colors: , "#66C2A5", "#BC80BD", "#8DA0CB","#FC8D62"
mycol=c("#7570B3", "#E5C494", "#FB8072", "#B15928", "red", "blue", "#FFFF99", "#A65628", "#FFD92F")

pdf(file='families.pdf',width=10,height=10)
plot(g,layout=l,rescale=F,vertex.label.cex=.5,xlim=range(l[,1]),  ylim=range(l[,2]),
edge.width=E(g)$weight/20,vertex.size=degree(g,mode='out')/5,
edge.curved=T,vertex.label=NA,vertex.color=mycol[cluster$membership],
margin=0,asp=0)
dev.off()

######## BELOW NOT NEEDED? #########

l=layout.fruchterman.reingold(g)
l=layout.norm(l,-1,1,-1,1)
comm.g = fastgreedy.community(g)

length(comm.g)#returns the number of communities
sizes(comm.g)#returns the community sizes, in the order of their ids
membership(comm.g)#gives the division of the vertices, into communities. 

#note: please change the number of colors depending on how many of communities you got
#spare colors: , "#66C2A5", "#BC80BD", "#8DA0CB","#FC8D62"
mycol=c("#7570B3", "#E5C494", "#FB8072", "#B15928", "red", "blue", "#FFFF99", "#A65628", "#FFD92F")

#save the motif communites to a text file
length=1
for (i in 1:length(comm.g)) {
    length=max(length,length(as.character(unlist(comm.g[i]))))
}
print(length)

cluster.df = data.frame(matrix('a', ncol = length, nrow = 1),stringsAsFactors=F)
for (i in 1:length(comm.g))
{
community=as.character(unlist(comm.g[i]))
length(community)=length
cluster.df=rbind(cluster.df,community)
}

cluster.df=cluster.df[-1,]

write.table(cluster.df,file="cluster.df.txt",quote=F,sep="\t",row.names=F,col.names=F,na='')

pdf(paste0('family','_all_no_comm','.pdf'),width=10,height=10)
plot(g,layout=l,rescale=F,vertex.label.cex=.5,xlim=range(l[,1]),  ylim=range(l[,2]),
edge.width=E(g)$weight/20,vertex.size=degree(g,mode='out')/5,
edge.curved=T,vertex.label=NA,vertex.color=mycol[cluster$membership],
margin=0,asp=0)
dev.off()

pdf(paste0('family','_all_commu','.pdf'),width=10,height=10)
plot(g,layout=l,rescale=F,vertex.label.cex=.5,xlim=range(l[,1]),  mark.groups = communities(comm.g), ylim=range(l[,2]),
edge.width=E(g)$weight/20,vertex.size=degree(g,mode='out')/5,
edge.curved=T,vertex.label=NA,vertex.color=mycol[cluster$membership],
margin=0,asp=0)
dev.off()


for (i in 1:length(comm.g)) {
    write.table(as.data.frame(comm.g[[i]]), file = paste0('PSWM_community_',i,'.txt'), col.names = FALSE, row.names = FALSE, quote = FALSE, sep = '\t')
}




