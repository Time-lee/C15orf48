library(stringr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(data.table)


#################======= step1: clean GTEx pheno data   =======#################
gtex <- read.table("samplepair.txt",header=T,sep='\t')

tcga_ref <- gtex[,1:2]

gtex$type <- paste0(gtex$TCGA,"_normal_GTEx")
gtex$sample_type <-"normal"
gtex <- gtex[,c("TCGA","GTEx","type","sample_type")]
names(gtex)[1:2] <- c("tissue","X_primary_site")


gp <- read.delim(file="GTEX_phenotype.gz",header=T,as.is = T)
gtex2tcga <- merge(gtex,gp,by="X_primary_site")
gtex_data <- gtex2tcga[,c(5,2:4)]
names(gtex_data)[1] <- "sample"
#write.table(gtex_data,"GTEx_pheno.txt",row.names=F,quote=F,sep='\t')

#################======= step2: clean a TCGA pheno data   =======#################
tcga <- read.delim(file="TCGA_phenotype_denseDataOnlyDownload.tsv.gz",header=T,as.is = T)
tcga <- merge(tcga_ref,tcga,by.y="X_primary_disease",by.x="Detail",all.y = T)
tcga <- tcga[tcga$sample_type %in% c("Primary Tumor","Solid Tissue Normal"),]

tcga$type <- ifelse(tcga$sample_type=='Solid Tissue Normal',
                    paste(tcga$TCGA,"normal_TCGA",sep="_"),paste(tcga$TCGA,"tumor_TCGA",sep="_"))
tcga$sample_type <- ifelse(tcga$sample_type=='Solid Tissue Normal',"normal","tumor")
tcga<-tcga[,c(3,2,6,5)]
names(tcga)[2] <- "tissue"
#write.table(tcga,"tcga_pheno.txt",row.names = F,quote=F,sep='\t')


#################======= step3: remove samples without tpm data =======############
gtex_exp <-  fread("gtex_RSEM_gene_tpm.gz",data.table = F)
gtexS <- gtex_data[ gtex_data$sample%in%colnames(gtex_exp)[-1],]

tcga_exp <- fread("tcga_RSEM_gene_tpm.gz",data.table = F)
tcgaS <- tcga[tcga$sample %in%colnames(tcga_exp)[-1],]
tcga_gtex <- rbind(tcgaS,gtexS)
write.table(tcga_gtex,"tcga_gtex_sample.txt",row.names = F,quote=F,sep='\t')

library(stringr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(data.table)
library(tibble)

rm(list = ls())
options(stringsAsFactors = FALSE) 

target <- "C15orf48"

idmap <- read.delim("gencode.v23.annotation.gene.probemap",as.is=T)
tcga_exp <- fread("tcga_RSEM_gene_tpm.gz",data.table = F)
gtex_exp <- fread("gtex_RSEM_gene_tpm.gz",data.table=F)
tcga_gtex <- read.table("tcga_gtex_sample.txt",sep='\t',header = T)

id <- idmap$id[which(idmap$gene==target)]
tcga_data <- t(tcga_exp[tcga_exp$sample==id,colnames(tcga_exp)%in%c("sample",tcga_gtex$sample)])
tcga_data <- data.frame(tcga_data[-1,])
tcga_data <- rownames_to_column(tcga_data,"sample")
names(tcga_data)[2] <- "tpm"

gtex_data <- t(gtex_exp[gtex_exp$sample==id,colnames(gtex_exp)%in%c("sample",tcga_gtex$sample)])
gtex_data <- data.frame(gtex_data[-1,])
gtex_data <- rownames_to_column(gtex_data,"sample")
names(gtex_data)[2] <- "tpm"

tmp  <- rbind(tcga_data,gtex_data)
exp <- merge(tmp,tcga_gtex,by="sample",all.x=T)
exp <- exp[,c("tissue","sample_type","tpm")]
exp <- arrange(exp,tissue)
write.table(exp,"C15orf48 expression.txt",row.names = F,quote=F,sep='\t')