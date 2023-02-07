
library(GSVA)
data=read.delim(file.choose(),row.names = 1,header = T)###The first line is ID, the first column is gene symbol
genelist=read.delim(file.choose(),row.names = 1,header=F)###The first line is genelist
genelist<-as.matrix(genelist) 
nrow=nrow(genelist)
data=data.matrix(data)
genelist[genelist==""]=NA  
genesets=c()

for (i in 1:nrow)
  
{
  a=genelist[i,]
  a=a[!is.na(a)]  
  a=list(as.character(as.matrix(a)))  
  genesets=c(genesets,a)
}

overlap_num<-c()
for(i in 1:nrow(genelist)){
  b<-as.character(unlist(genelist[i,]))
  o<-intersect(rownames(data),b)
  overlap_num<-c(overlap_num,length(o))
}

index<-overlap_num>0
d=rownames(genelist)[index]

result=gsva(data,genesets,mx.diff=FALSE, verbose=FALSE)
result=data.matrix(result)
rownames(result)=d
colnames(result)=colnames(data)

write.table(result,file="GSVA_result.txt",quote=F,sep="\t")
