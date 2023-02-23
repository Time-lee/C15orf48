#preSurvial

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")
library(limma)

expFile="GeneExp.txt"                   #expression data file
setwd("D:\\Lcl\\pancancer\\07.preSurvival")     

#Get all survival files in the directory
files=dir()
files=grep(".survival.tsv$",files,value=T)

#Read the survival data file for each tumor
surTime=data.frame()
for(i in 1:length(files)){
  inputFile=files[i]
  rt=read.table(inputFile,header=T,sep="\t",check.names=F,row.names=1)
  rt=rt[,c(3,1)]
  surTime=rbind(surTime,rt)
}
colnames(surTime)=c("futime","fustat")

#Read the expression file and organize the input file
exp=read.table(expFile,sep="\t",header=T,check.names=F,row.names=1)
exp=exp[(exp[,"Type"]=="Tumor"),]

#Combine data and output results
sameSample=intersect(row.names(surTime),row.names(exp))
surTime=surTime[sameSample,]
exp=exp[sameSample,]
surData=cbind(surTime,exp)
surData=cbind(id=row.names(surData),surData)
write.table(surData,file="expTime.txt",sep="\t",quote=F,row.names=F)

#Survival

#install.packages("survival")
#install.packages("survminer")

library(survival)
library(survminer)

setwd("D:\\Lcl\\pancancer\\07.preSurvival\\08.survival")                     #设置工作目录
rt=read.table("expTime.txt",header=T,sep="\t",check.names=F,row.names=1)       #读取输入文件
rt$futime=rt$futime/365
gene=colnames(rt)[3]
pFilter=0.05           

#Cycle through tumor types
for(i in levels(factor(rt[,"CancerType"]))){
  rt1=rt[(rt[,"CancerType"]==i),]
  group=ifelse(rt1[,gene]>median(rt1[,gene]),"high","low")
  diff=survdiff(Surv(futime, fustat) ~group,data = rt1)
  pValue=1-pchisq(diff$chisq,df=1)
  if(pValue<pFilter){
    if(pValue<0.001){
      pValue="p<0.001"
    }else{
      pValue=paste0("p=",sprintf("%.03f",pValue))
    }
    fit <- survfit(Surv(futime, fustat) ~ group, data = rt1)
    #plot survival curve
    surPlot=ggsurvplot(fit, 
                       data=rt1,
                       title=paste0("Cancer: ",i),
                       pval=pValue,
                       pval.size=6,
                       legend.labs=c("high","low"),
                       legend.title=paste0(gene," levels"),
                       font.legend=12,
                       break.time.by = 1,
                       palette=c("red","blue"),
                       conf.int=F,
                       fontsize=4,
                       risk.table=TRUE,
                       risk.table.title="",
                       ylab="Overall survival",
                       xlab="Time(years)",
                       risk.table.height=.25)
    pdf(file=paste0("survival.",i,".pdf"),onefile = FALSE,
        width = 6,             
        height =5)             
    print(surPlot)
    dev.off()
  }
}

