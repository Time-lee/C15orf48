#TMB
setwd("D:\\Lcl\\pancancer\\15.TMBcor")          

#read expression file
exp=read.table("singleGeneExp.txt", header=T,sep="\t",row.names=1,check.names=F)
#read TMB file
TMB=read.table("TMB.txt", header=T,sep="\t",row.names=1,check.names=F)
#Remove normal samples
group=sapply(strsplit(row.names(exp),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
exp=exp[group==0,]
#Sample intersection
sameSample=intersect(row.names(TMB),row.names(exp))
TMB=TMB[sameSample,]
exp=exp[sameSample,]

#correlation test
outTab=data.frame()
fmsbTab=data.frame()
for(i in levels(factor(exp[,"CancerType"]))){
  exp1=exp[(exp[,"CancerType"]==i),]
  TMB1=TMB[(TMB[,"CancerType"]==i),]
  x=as.numeric(TMB1[,1])
  y=as.numeric(exp1[,1])
  corT=cor.test(x,y,method="spearman")
  cor=corT$estimate
  pValue=corT$p.value
  sig=ifelse(pValue<0.001,"***",ifelse(pValue<0.01,"**",ifelse(pValue<0.05,"*"," ")))
  outTab=rbind(outTab,cbind(CancerType=i,cor=cor,pValue=pValue,sig))
  fmsbTab=rbind(fmsbTab,cbind(CancerType=i,cor=cor))
}
write.table(outTab,file="corStat.txt",sep="\t",row.names=F,quote=F)
write.table(t(fmsbTab),file="fmsbInput.txt",sep="\t",col.names=F,quote=F)

#install.packages("fmsb")

library(fmsb) 
setwd("D:\\Lcl\\pancancer\\16.TMBradar")        
data=read.table("fmsbInput.txt",header=T,sep="\t",row.names=1,check.names=F)   
maxValue=ceiling(max(abs(data))*10)/10
data=rbind(rep(maxValue,ncol(data)),rep(-maxValue,ncol(data)),data)
#define color
colors="red"
#define significance
corStat=read.table("corStat.txt",header=T,sep="\t",row.names=1,check.names=F)
colnames(data)=paste0(colnames(data),corStat$sig)

#output result
pdf(file="radar.pdf",height=7,width=7)
radarchart( data, axistype=1 , 
            pcol=colors,                 
            plwd=2 ,                     
            plty=1,                      
            cglcol="grey",              
            cglty=1,                      
            caxislabels=seq(-maxValue,maxValue,maxValue/2),   
            cglwd=1.2,                   
            axislabcol="blue",          
            vlcex=0.8                   
)
dev.off()


#MSI
setwd("D:\\Lcl\\pancancer\\17.MSIcor")            

#read expression file
exp=read.table("singleGeneExp.txt", header=T,sep="\t",row.names=1,check.names=F)
#Read MSI file
MSI=read.table("MSI.txt", header=T,sep="\t",row.names=1,check.names=F)
#Remove normal samples
group=sapply(strsplit(row.names(exp),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
exp=exp[group==0,]
#Sample intersection
sameSample=intersect(row.names(MSI),row.names(exp))
MSI=MSI[sameSample,]
exp=exp[sameSample,]

#correlation test
outTab=data.frame()
fmsbTab=data.frame()
for(i in levels(factor(exp[,"CancerType"]))){
  exp1=exp[(exp[,"CancerType"]==i),]
  MSI1=MSI[(MSI[,"CancerType"]==i),]
  x=as.numeric(MSI1[,1])
  y=as.numeric(exp1[,1])
  corT=cor.test(x,y,method="spearman")
  cor=corT$estimate
  pValue=corT$p.value
  sig=ifelse(pValue<0.001,"***",ifelse(pValue<0.01,"**",ifelse(pValue<0.05,"*"," ")))
  outTab=rbind(outTab,cbind(CancerType=i,cor=cor,pValue=pValue,sig))
  fmsbTab=rbind(fmsbTab,cbind(CancerType=i,cor=cor))
}
write.table(outTab,file="corStat.txt",sep="\t",row.names=F,quote=F)
write.table(t(fmsbTab),file="fmsbInput.txt",sep="\t",col.names=F,quote=F)

#install.packages("fmsb")

library(fmsb) 
setwd("D:\\Lcl\\pancancer\\18.MSIradar")        
data=read.table("fmsbInput.txt",header=T,sep="\t",row.names=1,check.names=F)  
maxValue=ceiling(max(abs(data))*10)/10
data=rbind(rep(maxValue,ncol(data)),rep(-maxValue,ncol(data)),data)
#define color
colors="blue"
#define significance
corStat=read.table("corStat.txt",header=T,sep="\t",row.names=1,check.names=F)
colnames(data)=paste0(colnames(data),corStat$sig)

#output result
pdf(file="radar.pdf",height=7,width=7)
radarchart( data, axistype=1 , 
            pcol=colors,                 
            plwd=2 ,                    
            plty=1,                      
            cglcol="grey",               
            cglty=1,                     
            caxislabels=seq(-maxValue,maxValue,maxValue/2),    
            cglwd=1.2,                   
            axislabcol="green",           
            vlcex=0.8                   
)
dev.off()
