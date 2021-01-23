setwd('/scratch/bhn9by/ATAC/SP_KLF_split/KLF')

#read in values object
values = read.table('composite.values.txt',sep='\t')
values$name <- values$V1
colnames(values) = c('count','afreq','cfreq','gfreq','tfreq','name')

#read in index object
index = read.table('composite.index.txt',sep='\t')
colnames(index) = c('count','rc','offset')
index = rbind(index,c('query','n',0))
index$offset = as.numeric(index$offset)

#add PSWM length and eventual starting row to index object
pswm.length = c()
for (count in index$count) {
    df = values[values$count == count,]
    pswm.length = append(pswm.length,nrow(df))
}
index$pswm.length = pswm.length

index$start.row = (-1*index$offset) + max(index$offset) + 1
index$end.row = index$start.row + index$pswm.length - 1 

#define how many bases to the left and right of the query PSWM
#for the left, it's just the maximum offset
#left.shift = max(index$offset)
#for the right, it's the minimum offset (i.e. the most negative value) plus the length of the longest PSWM associated with that offset minus the query length  
#right.shift = abs(min(index$offset)) + max(index[index$offset == min(index$offset),]$pswm.length) - index$pswm.length[nrow(index)]

#the final length of the composite is the length of the query plus the extra positions on the right and left
composite.length = max(index$end.row)

                                        #list of all the members of the family
#changed this
factors=as.vector(unique(values$name))



#define a dataframe for each base
#each row will be a position in the composite PSWM, and each column is the contribution from each family member
df.a = data.frame(matrix(nrow=composite.length,ncol=length(factors)))
colnames(df.a) = factors
df.c = df.g = df.t = df.a

#iterate through each factor
for (count in index$count) {

    #isolate values from each factor for convenience
    df = values[values$count == count,]

    #deal with RC
    if (index[index$count == count,]$rc == 'y') {
        df = df[order(nrow(df):1),]
        colnames(df) = c('count','name','tfreq','gfreq','cfreq','afreq')
    }

    #determine appropriate starting row for each PSWM in the final df
    start.row = index[index$count == count,]$start.row
    end.row = index[index$count == count,]$end.row
    #determine range of row values for the final df
    range = start.row:end.row
    #determine corresponding factor
    col = unique(as.vector(df$name))
    #iterate through relevant rows
    for (i in 1:length(range)) {
        #convert row numbers in final df to row numbers in original df
        j = range[i]
        #assign values in final dfs from orginial df
        df.a[j,col] = df$afreq[i]       
        df.c[j,col] = df$cfreq[i]
        df.g[j,col] = df$gfreq[i]
        df.t[j,col] = df$tfreq[i]
    }   
}

#initalize final composite PSWM
composite = data.frame(matrix(nrow=nrow(df.a),ncol=4))

#take medians of all the rows
med.a = c()
med.c = c()
med.g = c()
med.t = c()
for (i in 1:nrow(df.a)) {
    if (sum(is.na(df.a[i,]))==0) {
        med.a = append(med.a,median(as.numeric(df.a[i,]),na.rm=T))
    } else {
        x.a = df.a[i,]
        x.a[is.na(x.a)] <- 0.25
        med.a = append(med.a, mean(as.numeric(x.a)))
    }
    if (sum(is.na(df.c[i,]))==0) {
        med.c = append(med.c,median(as.numeric(df.c[i,]),na.rm=T))
    } else {
        x.c = df.c[i,]
        x.c[is.na(x.c)] <- 0.25
        med.c = append(med.c, mean(as.numeric(x.c)))
    }
    if (sum(is.na(df.g[i,]))==0) {
        med.g = append(med.g,median(as.numeric(df.g[i,]),na.rm=T))
    } else {
        x.g = df.g[i,]
        x.g[is.na(x.g)] <- 0.25
        med.g = append(med.g, mean(as.numeric(x.g)))
    }
    if (sum(is.na(df.t[i,]))==0) {
        med.t = append(med.t,median(as.numeric(df.t[i,]),na.rm=T))
    } else {
        x.t = df.t[i,]
        x.t[is.na(x.t)] <- 0.25
        med.t = append(med.t, mean(as.numeric(x.t)))
    }
}

composite[,1] = med.a
composite[,2] = med.c
composite[,3] = med.g
composite[,4] = med.t

#some rows don't add up to 1 so normalize to row sum (?)
composite = composite/rowSums(composite)

write.table(composite, file = paste0('KLF_composite_PSWM.txt'), sep='\t' ,row.names = F, col.names = F, quote = F)
