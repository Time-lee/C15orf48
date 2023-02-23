#library(utils)
#rforge <- "http://r-forge.r-project.org"
#install.packages("estimate", repos=rforge, dependencies=TRUE)

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")


library(limma)
library(estimate)
setwd("D:\\Lcl\\pancancer\\14.estimate")        

files=dir()
files=grep("^symbol.",files,value=T)

outTab=data.frame()
for(i in files){
  CancerType=gsub("symbol\\.|\\.txt","",i)
  rt=read.table(i,sep="\t",header=T,check.names=F)
  
  rt=as.matrix(rt)
  rownames(rt)=rt[,1]
  exp=rt[,2:ncol(rt)]
  dimnames=list(rownames(exp),colnames(exp))
  data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
  data=avereps(data)
  
  group=sapply(strsplit(colnames(data),"\\-"),"[",4)
  group=sapply(strsplit(group,""),"[",1)
  group=gsub("2","1",group)
  data=data[,group==0]
  out=data[rowMeans(data)>0,]
  out=rbind(ID=colnames(out),out)
  write.table(out,file="uniq.symbol.txt",sep="\t",quote=F,col.names=F)
  
  #Run the estimate package
  uniqFile="uniq.symbol.txt"
  inputDs="commonGenes.gct"
  outputDs="estimateScore.gct"
  filterCommonGenes(input.f=uniqFile, output.f=inputDs, id="GeneSymbol")
  estimateScore(input.ds =inputDs ,output.ds=outputDs)
  
  #Output the score for each sample
  scores=read.table("estimateScore.gct",skip = 2,header = T)
  rownames(scores)=scores[,1]
  scores=t(scores[,3:ncol(scores)])
  rownames(scores)=gsub("\\.","\\-",rownames(scores))
  outTab=rbind(outTab,cbind(scores,CancerType))
  file.remove(uniqFile)
  file.remove(inputDs)
  file.remove(outputDs)
}
out=cbind(ID=row.names(outTab),outTab)
write.table(out,file="estimateScores.txt",sep="\t",quote=F,row.names=F)


#Estimate Correlation
#install.packages("corrplot")

library(corrplot)                    
expFile="panGeneExp.txt"             
scoreFile="estimateScores.txt"       
scoreType="StromalScore"             
setwd("D:\\Lcl\\pancancer\\15.estimateCor")     

exp=read.table(expFile, header=T,sep="\t",row.names=1,check.names=F)
exp=exp[(exp[,"Type"]=="Tumor"),]

TME=read.table(scoreFile, header=T,sep="\t",row.names=1,check.names=F)

sameSample=intersect(row.names(TME),row.names(exp))
TME=TME[sameSample,]
exp=exp[sameSample,]

outTab=data.frame()
pTab=data.frame()
for(i in levels(factor(exp[,"CancerType"]))){
  exp1=exp[(exp[,"CancerType"]==i),]
  TME1=TME[(TME[,"CancerType"]==i),]
  x=as.numeric(TME1[,scoreType])
  pVector=data.frame(i)
  outVector=data.frame(i)
  genes=colnames(exp1)[1:(ncol(exp1)-2)]
  for(j in genes){
    y=as.numeric(exp1[,j])
    corT=cor.test(x,y,method="spearman")
    cor=corT$estimate
    pValue=corT$p.value
    pVector=cbind(pVector,pValue)
    outVector=cbind(outVector,cor)
  }
  pTab=rbind(pTab,pVector)
  outTab=rbind(outTab,outVector)
}
colNames=c("CancerType",colnames(exp1)[1:(ncol(exp1)-2)])
colnames(outTab)=colNames
write.table(outTab,file="estimateCor.cor.txt",sep="\t",row.names=F,quote=F)

colnames(pTab)=colNames
write.table(pTab,file="estimateCor.pval.txt",sep="\t",row.names=F,quote=F)