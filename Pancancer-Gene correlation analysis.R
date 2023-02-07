##===========
library(TCGAbiolinks)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(reshape2)
library(ggplot2)
library(tibble)
library(ggpubr)
library(pheatmap)
library(ggrepel)
source("00_function.R")
gene = "C15orf48"

outdir = paste0("outputResult/",gene)
ifelse(dir.exists(outdir),print("folder already exists"),dir.create(outdir))



#Read the input file and process the input file
TPMFilePath <- dir("./01-data/02-TCGA-RNASeq-TPM/","RNASeq_TPM.RData$",full.names = T)
project <- getGDCprojects()$project_id
project <- project[grep("TCGA-",project)]

genelist = readLines("02-base_files/inhibitory immune checkpoints.txt")



###Calculate correlation
all_cancer_gene_cor = data.frame()
for(proj in project){
  if(proj == project[1]){
    message("========================================================")
  }
  cancer =  unlist(strsplit(proj,"-"))[2]
  load(TPMFilePath[grep(proj,TPMFilePath)])#expDataTPM
  
  SamN <- TCGAquery_SampleTypes(barcode = colnames(expDataTPM),typesample = "NT")
  SamT <- setdiff(colnames(expDataTPM),SamN)
  tur_expDataTPM = expDataTPM[,SamT]
  
  congene = intersect(genelist,rownames(tur_expDataTPM))
  if(proj == project[1]){
    message(paste0("============================\nNo mapped gene(s):",
                   setdiff(genelist,congene)),"\n============================")
  }
  message(proj)
  cordf = data.frame()
  for(g in congene){
    cor.spearman = cor.test(as.numeric(tur_expDataTPM[g,]),
                            as.numeric( tur_expDataTPM[gene,]),
                            method ="spearman")
    cor.pearson = cor.test(as.numeric(tur_expDataTPM[g,]),
                           as.numeric( tur_expDataTPM[gene,]),
                           method ="pearson")
    df = data.frame(cancer = cancer,
                    gene1 = gene,
                    gene2 = g,
                    spearman = cor.spearman[["estimate"]][["rho"]],
                    spearman.p.value = cor.spearman[["p.value"]],
                    pearson = cor.pearson[["estimate"]][["cor"]],
                    pearson.p.value = cor.pearson[["p.value"]])
    cordf = rbind(cordf,df)
  }
  all_cancer_gene_cor = rbind(all_cancer_gene_cor,cordf)
  if(proj == project[33]){
    message("========================================================")
  }
}



pdf(paste0(outdir,"/",gene,"correlation analysis.pdf"),height=5,width=5)
for(g2 in unique(all_cancer_gene_cor$gene2)){
  cod = all_cancer_gene_cor[all_cancer_gene_cor$gene2 == g2,]
  
  rho = cod[,c(1,4,5)]
  rhocancer = rho$cancer[abs(rho$spearman)> 0.2 & rho$spearman.p.value < 0.01]
  rho$label = ""
  rho$label[match(rhocancer,rho$cancer)] = rhocancer
  rho$Group = ifelse(abs(rho$spearman)> 0.2 & rho$spearman.p.value < 0.01,
                     ifelse(rho$spearman > 0.2,"positive","negative"),
                     "No")
  
  cor = cod[,c(1,6,7)]
  corcancer = cor$cancer[abs(cor$pearson)> 0.2 & cor$pearson.p.value < 0.01]
  cor$label = ""
  cor$label[match(corcancer,cor$cancer)] = corcancer
  cor$Group = ifelse(abs(cor$pearson)> 0.2 & cor$pearson.p.value < 0.01,
                     ifelse(cor$pearson > 0.2,"positive","negative"),
                     "No")
  
  ###=======pearson
  p1 = ggplot(data = cor, aes(x = pearson, 
                              y = -log10(pearson.p.value),
                              colour = Group)) + 
    geom_point(alpha = 0.9,shape = 19,size=3) +
    scale_color_manual(values = c("#00008B", "#808080","#DC143C")) + 
    theme_bw() +
    geom_text_repel(label = rho$label,
                    size = 5,
                    segment.color = "black", 
                    show.legend = FALSE) +
    labs(title = paste0())+
    theme(axis.title=element_text(size=15,face="plain",color="black"),
          axis.text = element_text(size=12,face="plain",color="black"),
          legend.position =  "none",
          panel.background = element_rect(fill = "transparent",colour = "black"),
          plot.background = element_blank(),
          plot.title = element_text(size=15, lineheight=.8,hjust=0.5, face="plain"),
          legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))+
    ylab(expression(-log[10]("p.value"))) +
    xlab(paste0("Cor of ",gene," and ",g2)) +
    geom_vline(xintercept = c(-0.2,0.2),
               lty = 2,
               col = "black",
               lwd = 0.3) +
    geom_hline(yintercept = -log10(0.01),
               lty = 2,
               col = "black",
               lwd = 0.3) 
  
  ##==========spearman
  p2 =  ggplot(data = rho, aes(x = spearman, 
                               y = -log10(spearman.p.value),
                               colour = Group)) + 
    geom_point(alpha = 0.9,shape = 19,size=3) +
    scale_color_manual(values = c("#00008B", "#808080","#DC143C")) + 
    theme_bw() +
    geom_text_repel(label = rho$label,
                    size = 5,
                    segment.color = "black", 
                    show.legend = FALSE) +
    labs(title = paste0())+
    theme(axis.title=element_text(size=15,face="plain",color="black"),
          axis.text = element_text(size=12,face="plain",color="black"),
          legend.position =  "none",
          panel.background = element_rect(fill = "transparent",colour = "black"),
          plot.background = element_blank(),
          plot.title = element_text(size=15, lineheight=.8,hjust=0.5, face="plain"),
          legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))+
    ylab(expression(-log[10]("p.value"))) +
    xlab(paste0("Rho of ",gene," and ",g2)) +
    geom_vline(xintercept = c(-0.2,0.2),
               lty = 2,
               col = "black",
               lwd = 0.3) +
    geom_hline(yintercept = -log10(0.01),
               lty = 2,
               col = "black",
               lwd = 0.3)
  print(p1)
  print(p2)
}
dev.off()

color.key<-c("#3300CC","#3399FF","white","#FF3333","#CC0000")
f1 = paste0(outdir,"/",gene,"-heatmap.pdf")
pdf(f1,height=6,width=9)
{
  rho = all_cancer_gene_cor[,c(1,3,4)]
  rhop = all_cancer_gene_cor[,c(1,3,5)]
  rho_wdata = spread(rho,cancer,spearman,fill = NA)
  rownames(rho_wdata) = rho_wdata$gene2
  rho_wdata = rho_wdata[,-1]
  
  rho_pdata = spread(rhop,cancer, spearman.p.value,fill = NA)
  rownames(rho_pdata) = rho_pdata$gene2
  rho_pdata = rho_pdata[,-1]
  display_numbers = matrix(ifelse(rho_pdata > 0.01 | is.na(rho_pdata), "×", ""), nrow(rho_pdata))
  pheatmap(rho_wdata,
           #annotation_col =annotation_col,
           color =colorRampPalette(color.key)(50),
           cluster_cols =T,
           fontsize=8,
           cluster_rows = F,
           fontsize_row=8,
           cellwidth =12,
           cellheight =12,
           fontsize_number = 16,
           display_numbers= display_numbers,
           #scale="row",
           show_colnames=T,
           fontsize_col=8,
           main = paste0("Spearman correlation analysis of ",
                         gene," and ","immune checkpoint",
                         "\nfrom the The Cancer Genome Atlas database")) %>% print()
  cor_wdata = spread(all_cancer_gene_cor[,c(1,3,6)],cancer,pearson,fill = NA)
  rownames(cor_wdata) = cor_wdata$gene2
  cor_wdata = cor_wdata[,-1]
  
  cor_pdata = spread(all_cancer_gene_cor[,c(1,3,7)],cancer, pearson.p.value,fill = NA)
  rownames(cor_pdata) = cor_pdata$gene2
  cor_pdata = cor_pdata[,-1]
  display_numbers = matrix(ifelse(cor_pdata > 0.01 | is.na(cor_pdata), "×", ""), nrow(cor_pdata))
  pheatmap(cor_wdata,
           #annotation_col =annotation_col,
           color =colorRampPalette(color.key)(50),
           cluster_cols =T,
           fontsize=8,
           cluster_rows = F,
           fontsize_row=8,
           cellwidth =12,
           cellheight =12,
           fontsize_number = 16,
           display_numbers= display_numbers,
           #scale="row",
           show_colnames=T,
           fontsize_col=8,
           main = paste0("Pearson correlation analysis of ",
                         gene," and ","immune checkpoint",
                         "\nfrom the The Cancer Genome Atlas database")) %>% print()
}
dev.off()
