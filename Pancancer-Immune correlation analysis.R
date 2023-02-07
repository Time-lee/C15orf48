
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

timer = read.table("01-data/04-infiltrationTIMER2/infiltration_estimation_for_tcga.csv",
                   sep = ",",
                   header = T,
                   check.names = F)
#Read the input file and process the input file
TPMFilePath <- dir("./01-data/02-TCGA-RNASeq-TPM/","RNASeq_TPM.RData$",full.names = T)
project <- getGDCprojects()$project_id
project <- project[grep("TCGA-",project)]




all_immue_eatimate = data.frame()
allcancer_ssGSEA_Score = data.frame()
gene_immCell_cor = data.frame()#ssGSEA
gene_immCell_cor2 = data.frame()#TIMER2
#proj = "TCGA-LAML"


###====================ssGSEA，TIMER
outdir = paste0("outputResult/",gene)
ifelse(dir.exists(outdir),print("folder already exists"),dir.create(outdir))


color.key<-c("#3300CC","#3399FF","white","#FF3333","#CC0000")
#proj = "TCGA-LUAD"
for(proj in project){
  message(proj)
  cancer =  unlist(strsplit(proj,"-"))[2]
  load(TPMFilePath[grep(proj,TPMFilePath)])#expDataTPM
  
  SamN <- TCGAquery_SampleTypes(barcode = colnames(expDataTPM),typesample = "NT")
  ##tumor tissue sample
  SamT <- setdiff(colnames(expDataTPM),SamN)
  tur_expDataTPM = expDataTPM[,SamT]

  ssGSEA_Score = read.table(paste0("01-data/03-TCGA_turmor_ssGSEA/",proj,"-norm_turmor_ssGSEA_Score.txt"),
                            header = T,sep = "\t",check.names = F)
  rownames(ssGSEA_Score) = ssGSEA_Score$id
  med = median(tur_expDataTPM[gene,])
  gene_meta = data.frame(tpm = tur_expDataTPM[gene,],
                         group = ifelse(tur_expDataTPM[gene,] > med,"High","Low"))
  
  coid = intersect(rownames(gene_meta), colnames(ssGSEA_Score))
  gene_meta = gene_meta[coid,]
  gene_meta = arrange(gene_meta,group)
  ssGSEA_Score = ssGSEA_Score[,rownames(gene_meta)]
  outdirp = paste0(outdir,"/",proj)
  ifelse(dir.exists(outdirp),
         print("folder already exists"),
         dir.create(outdirp))
  
  samples_meta = data.frame(gene =  gene_meta[,2])
  rownames(samples_meta) = rownames(gene_meta)
  ann_colors = list(gene = c(High = "#E41A1C",Low = "#4DAF4A"))
  
  
  f1 = paste0(outdirp,"/",proj,"-ssGSEA-heatmap.pdf")
  pdf(f1,height=3.5,width=9)
  pheatmap(ssGSEA_Score,
           annotation_col = samples_meta,
           color =colorRampPalette(color.key)(50),
           cluster_cols =F,
           fontsize=8,
           fontsize_row=8,
           scale="row",
           show_colnames=F,
           annotation_colors = ann_colors,
           fontsize_col=3)
  dev.off()
  norm_ssGSEA_Score2 =  ssGSEA_Score %>% t() %>% as.data.frame()
  norm_ssGSEA_Score2 = cbind(samples_meta,norm_ssGSEA_Score2)
  
  ldata <- melt(norm_ssGSEA_Score2,id.vars = "gene")
  head(ldata)
  
  p = ggplot(ldata, aes(x=variable, y=value,color = gene)) +
    #geom_jitter(size = 0.5,aes(x=variable, y=value,color = Group))+
    geom_boxplot(aes(color = gene),alpha =1,
                 lwd = 0.5, outlier.size = 1,
                 outlier.colour = "white")+ #color = c("red", "blue"),
    theme_bw()+
    stat_compare_means(label = "p.signif") +
    #labs(title = 'ImmuneCellAI') +
    theme(legend.position = "top",
          plot.title = element_text(color = 'black', size = 12, hjust = 0.5),
          axis.text.x = element_text(angle = 45,face = "bold",colour = "#1A1A1A",hjust=1,vjust=1),
          axis.text.y = element_text(face = "bold",colour = "#1A1A1A"),
          axis.title = element_text(size = 12,face = "bold", colour = "#1A1A1A"),
          legend.title =element_blank(),
          legend.text = element_text(size = 12, face = "bold",colour = "#1A1A1A")
    ) +
    ylab("ssGSEA score") +
    xlab("") +
    #ylim(c(0.45,1.6))+
    scale_color_manual(values = c("#E41A1C","#3300CC")) #+ ,
  pdf(paste0(outdirp,"/",proj,"-ssGSEA-boxplot.pdf"),height=5.5,width=9)
  print(p)
  dev.off()
  
  ssGSEA_Score = ssGSEA_Score %>% t() %>% as.data.frame()
  ssGSEA_Score = rownames_to_column(ssGSEA_Score,var = "id")
  ssGSEA_Score = add_column(ssGSEA_Score,cancer = proj,.before = 1)
  gene_meta = rownames_to_column(gene_meta,var = "id")
  ssGSEA_Score = merge(gene_meta,ssGSEA_Score,by = "id")
  allcancer_ssGSEA_Score = rbind(allcancer_ssGSEA_Score, ssGSEA_Score)
  
  tpm = ssGSEA_Score$tpm
  one_cell_df = data.frame()
  for(cell in colnames(ssGSEA_Score)[-c(1:4)]){
    immScore = ssGSEA_Score[,cell]
    cor.spearman = cor.test(tpm,immScore,method ="spearman")
    df = data.frame(cancer = cancer,cell = cell,
                    spearman= cor.spearman[["estimate"]],
                    p.value = cor.spearman[["p.value"]])
    one_cell_df = rbind(one_cell_df,df)
  }
  
  gene_immCell_cor = rbind(gene_immCell_cor,one_cell_df)
  
  
  
  ##-----------------------------
  samples_meta$gene = factor(samples_meta$gene,levels = c("High","Low"))
  
  id = gsub("(.*?)-(.*?)-(.*?)-([0-9][1-9])[A-Z]-.*","\\1-\\2-\\3-\\4",rownames(samples_meta))
  samples_meta$cell_type = id
  samples_meta = samples_meta[!duplicated(samples_meta$cell_type),]
  
  conp = intersect(samples_meta$cell_type,timer$cell_type)
  onetimer = timer[timer$cell_type %in% conp,]
  onetimer = merge(samples_meta,onetimer,by = "cell_type")
  rownames(onetimer) = onetimer$cell_type
  immue_eatimate = add_column(onetimer,cancer = proj,.before = 1)
  colnames(immue_eatimate)[2] = "id"
  
  gene_meta$id = gsub("(.*?)-(.*?)-(.*?)-([0-9][1-9])[A-Z]-.*","\\1-\\2-\\3-\\4",gene_meta$id)
  immue_eatimate = merge(gene_meta,immue_eatimate,by = "id")
  immue_eatimate = immue_eatimate[,-5]
  all_immue_eatimate = rbind(all_immue_eatimate,immue_eatimate)
  
  iedf = data.frame()
  immue_eatimate = immue_eatimate %>% t() %>% na.omit() %>% t() %>% as.data.frame()
  for(cell_m in colnames(immue_eatimate)[-c(1:4)]){
    cor.spearman = cor.test(as.numeric(immue_eatimate[,"tpm"]),
                            as.numeric( immue_eatimate[,cell_m]),
                            method ="spearman")
    na.omit(c(NA,12,431))
    df = data.frame(cancer = cancer,
                    cell = unlist(strsplit(cell_m,"_" ))[1],
                    method = unlist(strsplit(cell_m,"_" ))[2],
                    spearman= cor.spearman[["estimate"]],
                    p.value = cor.spearman[["p.value"]])
    iedf = rbind(iedf,df)
  }
  gene_immCell_cor2 = rbind(gene_immCell_cor2,iedf )
  
  onetimer = onetimer[,-1]
  pdf(paste0(outdirp,"/",proj,"-TIMER2-CIBERSORT-boxplot.pdf"),height=5.5,width=9)
  print(plotbox_TIMER(data = onetimer,meth = "CIBERSORT",proj = proj)) %>% print()
  dev.off()
  pdf(paste0(outdirp,"/",proj,"-TIMER2-CIBERSORT-ABS-boxplot.pdf"),height=5.5,width=9)
  plotbox_TIMER(data = onetimer,meth = "CIBERSORT-ABS",proj = proj) %>% print()
  dev.off()
  pdf(paste0(outdirp,"/",proj,"-TIMER2-XCELL-boxplot.pdf"),height=5.5,width=15)
  plotbox_TIMER(data = onetimer,meth = "XCELL",proj = proj)%>% print()
  dev.off()
  
  pdf(paste0(outdirp,"/",proj,"-TIMER2-MCPCOUNTERboxplot.pdf"),height=5.5,width=7)
  plotbox_TIMER(data = onetimer,meth = "MCPCOUNTER",proj = proj)%>% print()
  dev.off()
  
  pdf(paste0(outdirp,"/",proj,"-TIMER2-QUANTISEQboxplot.pdf"),height=5.5,width=7)
  plotbox_TIMER(data = onetimer,meth = "QUANTISEQ",proj = proj)%>% print()
  dev.off()
  pdf(paste0(outdirp,"/",proj,"-TIMER2-EPICboxplot.pdf.pdf"),height=5.5,width=6)
  plotbox_TIMER(data = onetimer,meth = "EPIC",proj = proj)%>% print()
  dev.off()
  pdf(paste0(outdirp,"/",proj,"-TIMER2-TIMERboxplot.pdf.pdf"),height=5.5,width=6)
  plotbox_TIMER(data = onetimer,meth = "TIMER",proj = proj)%>% print()
  dev.off()
}

write.csv(all_immue_eatimate,
          file = paste0(outdir,"/all_immue_eatimate.csv"))
write.csv(allcancer_ssGSEA_Score,
          file = paste0(outdir,"/allcancer_ssGSEA_Score.csv"))
write.csv(gene_immCell_cor,
          file = paste0(outdir,"/gene_immCell_ssGSEA_cor.csv"))
write.csv(gene_immCell_cor2,
          file = paste0(outdir,"/gene_immCell_TIMER_cor.csv"))



####=============ssGSEA
library(ggplot2)
library(ggpubr)
library(ggrepel)
unique(gene_immCell_cor$cell)

pdf(paste0(outdir,"/",gene,"-ssGSEA_Scatterplot.pdf"),
    height=3.5,width=3.5)
#cell = "Plasmacytoid dendritic cell"
for(cell in unique(gene_immCell_cor$cell)){
  message(cell)
  onedata = gene_immCell_cor[gene_immCell_cor$cell == cell,]
  onedata$label = ""
  lab = onedata$cancer[abs(onedata$spearman)> 0.3 & onedata$p.value < 0.01]
  onedata$label[match(lab,onedata$cancer)] = lab
  onedata$Group = ifelse(abs(onedata$spearman)> 0.3 & onedata$p.value < 0.01,
                         ifelse(onedata$spearman > 0.3,"positive","negative"),
                         "No")
  onedata = arrange(onedata,Group)
  onedata = onedata[onedata$p.value!= "Inf",]
  p = ggplot(data = onedata, aes(x = spearman, y = -log10(p.value + 0.0001),colour = Group)) + 
    geom_point(alpha = 0.9,shape = 19,size=3) +
    scale_color_manual(values = c("#00008B", "#808080","#DC143C")) + 
    theme_bw() +#
    geom_text_repel(label = onedata$label,
                    size = 5,
                    segment.color = "black", 
                    show.legend = FALSE) +
    labs(title = cell)+
    theme(axis.title=element_text(size=15,face="plain",color="black"),
          axis.text = element_text(size=12,face="plain",color="black"),
          legend.position =  "none",
          panel.background = element_rect(fill = "transparent",colour = "black"),
          plot.background = element_blank(),
          plot.title = element_text(size=15, lineheight=.8,hjust=0.5, face="plain"),
          legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))+
    ylab(expression(-log[10]("p.value"))) +
    xlab("spearman correlation (r)") +
    geom_vline(xintercept = c(-0.3,0.3),
               lty = 2,
               col = "black",
               lwd = 0.3) +
    geom_hline(yintercept = -log10(0.01),
               lty = 2,
               col = "black",
               lwd = 0.3)
  print(p)
}
dev.off()

###================================
head(gene_immCell_cor)
imm_wdata = spread(gene_immCell_cor[,1:3],cancer,spearman,fill = NA)
rownames(imm_wdata) = imm_wdata$cell
imm_wdata = imm_wdata[,-1]

imm_pdata = spread(gene_immCell_cor[,-3],cancer, p.value,fill = NA)
rownames(imm_pdata) = imm_pdata$cell
imm_pdata = imm_pdata[,-1]
display_numbers = matrix(ifelse(imm_pdata > 0.05 | is.na(imm_pdata), "×", ""), nrow(imm_pdata))
pdf(paste0(outdir,"/",gene,"-ssGSEA_heatmap.pdf"),height=10,width=10)
pheatmap(imm_wdata,
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
         fontsize_col=8)
dev.off()


###= 6 algorithms to calculate immune cell infiltration
for(method in unique(gene_immCell_cor2$method)){
  message(method)
  methoddata = gene_immCell_cor2[gene_immCell_cor2$method == method,]
  file = paste0(outdir,"/",gene,"-",method,"-cor0.3.pdf")
  pdf(file,height=3.5,width=3.5)
  for(cell in unique(methoddata$cell)){
    onedata = methoddata[methoddata$cell == cell,]
    onedata = na.omit(onedata)
    onedata$label = ""
    lab = onedata$cancer[abs(onedata$spearman)> 0.3 & onedata$p.value < 0.01]
    onedata$label[match(lab,onedata$cancer)] = lab
    onedata$Group = ifelse(abs(onedata$spearman)> 0.3 & onedata$p.value < 0.01,
                           ifelse(onedata$spearman > 0.3,"positive","negative"),
                           "No")
    onedata = arrange(onedata,Group)
    p = ggplot(data = onedata, aes(x = spearman, y = -log10(p.value),colour = Group)) + 
      geom_point(alpha = 0.9,shape = 19,size=3) +
      scale_color_manual(values = c("#00008B", "#808080","#DC143C")) + 
      theme_bw() +
      geom_text_repel(label = onedata$label,
                      size = 5,
                      segment.color = "black", 
                      show.legend = FALSE) +
      labs(title = cell)+
      theme(axis.title=element_text(size=15,face="plain",color="black"),
            axis.text = element_text(size=12,face="plain",color="black"),
            legend.position =  "none",
            panel.background = element_rect(fill = "transparent",colour = "black"),
            plot.background = element_blank(),
            plot.title = element_text(size=15, lineheight=.8,hjust=0.5, face="plain"),
            legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))+
      ylab(expression(-log[10]("p.value"))) +
      xlab("spearman correlation (r)") +
      geom_vline(xintercept = c(-0.3,0.3),
                 lty = 2,
                 col = "black",
                 lwd = 0.3) +
      geom_hline(yintercept = -log10(0.01),
                 lty = 2,
                 col = "black",
                 lwd = 0.3)
    print(p)
  }
  dev.off()
}


for(method in unique(gene_immCell_cor2$method)){
  message(method)
  methoddata = gene_immCell_cor2[gene_immCell_cor2$method == method,]
  head(methoddata)
  meth_wdata = spread(methoddata[,c(1,2,4)],cancer,spearman,fill = NA)
  rownames(meth_wdata) = meth_wdata$cell
  meth_wdata = meth_wdata[,-1]
  
  meth_pdata = spread(methoddata[,c(1,2,5)],cancer, p.value,fill = NA)
  rownames(meth_pdata) = meth_pdata$cell
  meth_pdata = meth_pdata[,-1]
  dn = matrix(ifelse(meth_pdata > 0.05 | is.na(meth_pdata), "×", ""), nrow(meth_pdata))
  file = paste0(outdir,"/",gene,"-",method,"-heatmap.pdf")
  pdf(file,height=10,width=10)
  pheatmap(meth_wdata,
           #annotation_col =annotation_col,
           color =colorRampPalette(color.key)(50),
           cluster_cols =T,
           fontsize=8,
           cluster_rows = F,
           fontsize_row=8,
           cellwidth =12,
           cellheight =12,
           fontsize_number = 16,
           display_numbers= dn,
           #scale="row",
           show_colnames=T,
           fontsize_col=8)
  dev.off()
}


###================TIP================


activation = dir("01-data/05-TIP/Immune activity scores",".txt$",full.names = T)
infiltration = dir("01-data/05-TIP/Immune cell infiltration",".txt$",full.names = T)


all_cancer_infi_Score = data.frame()
all_cancer_inf_gene_cor = data.frame()
all_cancer_acti_Score = data.frame()
for(proj in project){
  message(proj)
  load(TPMFilePath[grep(proj,TPMFilePath)])#expDataTPM
  SamN <- TCGAquery_SampleTypes(barcode = colnames(expDataTPM),typesample = "NT")
  
  SamT <- setdiff(colnames(expDataTPM),SamN)
  tur_expDataTPM = expDataTPM[,SamT]
  
  med = median(tur_expDataTPM[gene,])
  gene_meta = data.frame(tpm = tur_expDataTPM[gene,],
                         group = ifelse(tur_expDataTPM[gene,] > med,"High","Low"))
  infi_Score = read.table(paste0("01-data/05-TIP/Immune cell infiltration/",proj,".txt"),
                          header = T,sep = "\t",check.names = F)
  infi_Score = infi_Score[,-c(16:18)]
  # rownames(infi_Score) = infi_Score$Mixture
  # infi_Score = infi_Score[,-1] %>% t() %>% as.data.frame()
  
  coid = intersect(rownames(gene_meta), infi_Score$Mixture)
  gene_meta = gene_meta[coid,]
  gene_meta = arrange(gene_meta,group)
  gene_meta = add_column(gene_meta,cancer = unlist(strsplit(proj,"-"))[2],.before = 1)
  gene_meta = add_column(gene_meta,Mixture = rownames(gene_meta),.before = 1)
  
  infi_Score = merge(gene_meta,infi_Score,by = "Mixture")
  colnames(infi_Score)[1] = "id"
  
  all_cancer_infi_Score = rbind(all_cancer_infi_Score,infi_Score)
  
  tpm = infi_Score$tpm
  ifdf = data.frame()
  for(cell in colnames(infi_Score)[-c(1:4)]){
    isco = infi_Score[,cell]
    if(var(isco)!=0){
      cor.spearman = cor.test(tpm,
                              infi_Score[,cell],
                              method ="spearman")
      
      df = data.frame(cancer = unlist(strsplit(proj,"-"))[2],
                      cell = cell,
                      spearman= cor.spearman[["estimate"]],
                      p.value = cor.spearman[["p.value"]])
      ifdf = rbind(ifdf,df)
    }
    
  }
  all_cancer_inf_gene_cor = rbind( all_cancer_inf_gene_cor,ifdf)
  
  
  ### activation score==========
  acti_Score = read.table(paste0("01-data/05-TIP/Immune activity scores/",proj,".txt"),
                          header = T,sep = "\t",check.names = F)
  rownames(acti_Score) = acti_Score$Steps
  acti_Score = acti_Score[,-1] %>% t() %>% as.data.frame()
  acti_Score = add_column(acti_Score,Mixture = rownames(acti_Score),.before = 1)
  acti_Score = merge(gene_meta,acti_Score,by = "Mixture")
  colnames(acti_Score)[1] = "id"
  all_cancer_acti_Score = rbind(all_cancer_acti_Score,acti_Score)
}
write.csv(all_cancer_infi_Score,
          file = paste0(outdir,"/",gene,"-all_cancer_infi_Score.csv"))
write.csv(all_cancer_inf_gene_cor,
          file = paste0(outdir,"/",gene,"-all_cancer_inf_gene_cor.csv"))
write.csv(all_cancer_acti_Score,
          file = paste0(outdir,"/",gene,"-all_cancer_acti_Score.csv"))

head(all_cancer_acti_Score)
ldata = all_cancer_acti_Score[,-c(1,3)]
ldata = gather(ldata,Step,activity_scores,Step1:Step7)
colnames(ldata)

head(ldata)
color.key<-c("#3300CC","#3399FF","white","#FF3333","#CC0000")
for(cancer in unique(ldata$cancer)){
  message(cancer)
  data= ldata[ldata$cancer == cancer,]
  p = ggplot(data, aes(x=Step, y=activity_scores,color = group)) +
    #geom_jitter(size = 0.5,aes(x=variable, y=value,color = Group))+
    geom_boxplot(aes(color = group),alpha =1,
                 lwd = 0.5, outlier.size = 1,
                 outlier.colour = "white")+ #color = c("red", "blue"),
    theme_bw()+
    stat_compare_means(label = "p.signif") +
    #labs(title = 'ImmuneCellAI') +
    theme(legend.position = "top",
          plot.title = element_text(color = 'black', size = 12, hjust = 0.5),
          axis.text.x = element_text(angle = 45,face = "bold",colour = "#1A1A1A",hjust=1,vjust=1),
          axis.text.y = element_text(face = "bold",colour = "#1A1A1A"),
          axis.title = element_text(size = 12,face = "bold", colour = "#1A1A1A"),
          legend.text = element_text(size = 12, face = "bold",colour = "#1A1A1A")
    ) +
    ylab("Immune activity scores") +
    xlab("") +
    #ylim(c(0.45,1.6))+
    scale_color_manual(values = c("#CC0000","#3300CC")) #+ ,"#EAE838"
  pdf(paste0(outdir,"/TCGA-",cancer,"/",gene,"-Immune activity scores-boxplot.pdf"),
      height=5.5,width=9)
  print(p)
  dev.off()
}

###TIP-heatmap
library(tidyr)

all_cancer_inf_gene_cor = read.csv(paste0(outdir,"/",gene,"-all_cancer_inf_gene_cor.csv"),
                                   header = T,row.names = 1)

unique(all_cancer_inf_gene_cor$cell)
inf_wdata = spread(all_cancer_inf_gene_cor[,1:3],cell,spearman)
rownames(inf_wdata) = inf_wdata$cancer
inf_wdata = inf_wdata[,-1]

inf_pval = spread(all_cancer_inf_gene_cor[,-3],cell,p.value)
rownames(inf_pval) = inf_pval$cancer
inf_pval = inf_pval[,-1]
inf_wdata = t(inf_wdata)

color.key<-c("#3300CC","#3399FF","white","#FF3333","#CC0000")
dn = matrix(ifelse(inf_pval > 0.05 | is.na(inf_pval), "×", ""), nrow(inf_pval))
dn = t(dn)
file = paste0(outdir,"/",gene,"-Immune cell infiltration-TIP-heatmap.pdf")
pdf(file,height=10,width=10)
pheatmap(inf_wdata,
         #annotation_col =annotation_col,
         color =colorRampPalette(color.key)(50),
         cluster_cols =T,
         fontsize=8,
         cluster_rows = F,
         fontsize_row=8,
         cellwidth =12,
         cellheight =12,
         fontsize_number = 16,
         display_numbers= dn,
         #scale="row",
         show_colnames=T,
         fontsize_col=8)
dev.off()

##===========immune modulator
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


geneType = read.table("02-base_files/immunomodulatorType.txt",header = T)


###Calculate correlation
all_cancer_gene_cor = data.frame()
for(proj in project){
  
  cancer =  unlist(strsplit(proj,"-"))[2]
  load(TPMFilePath[grep(proj,TPMFilePath)])#expDataTPM
  
  SamN <- TCGAquery_SampleTypes(barcode = colnames(expDataTPM),typesample = "NT")
  SamT <- setdiff(colnames(expDataTPM),SamN)
  tur_expDataTPM = expDataTPM[,SamT]
  
  congene = intersect(geneType$Gene,rownames(tur_expDataTPM))
  if(proj == project[1]){
    message(paste0("============================\nNo mapped gene(s):",
                   setdiff(geneType$Gene,congene)),"\n============================")
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
}


color.key<-c("#3300CC","#3399FF","white","#FF3333","#CC0000")
annotation_row = data.frame(Type = geneType$Type)
rownames(annotation_row) = geneType$Gene
annotation_row = arrange(annotation_row,Type)


f1 = paste0(outdir,"/",gene,"-heatmap.pdf")
pdf(f1,height=19,width=8)
{
  rho = all_cancer_gene_cor[,c(1,3,4)]
  rhop = all_cancer_gene_cor[,c(1,3,5)]
  rho_wdata = spread(rho,cancer,spearman,fill = NA)
  rownames(rho_wdata) = rho_wdata$gene2
  rho_wdata = rho_wdata[,-1]
  rho_wdata = rho_wdata[rownames(annotation_row),]
  
  rho_pdata = spread(rhop,cancer, spearman.p.value,fill = NA)
  rownames(rho_pdata) = rho_pdata$gene2
  rho_pdata = rho_pdata[,-1]
  display_numbers = matrix(ifelse(rho_pdata > 0.01 | is.na(rho_pdata), "×", ""), 
                           nrow(rho_pdata))
  rho_pdata= rho_pdata[rownames(annotation_row),]
  pheatmap(rho_wdata,
           #annotation_col =annotation_col,
           annotation_row = annotation_row,
           color =colorRampPalette(color.key)(50),
           cluster_cols =T,
           fontsize=6,
           cluster_rows = F,
           fontsize_row=6,
           cellwidth =10,
           cellheight =10,
           fontsize_number = 16,
           display_numbers= display_numbers,
           #scale="row",
           show_colnames=T,
           fontsize_col=8,
           main = paste0("Spearman correlation analysis of ",
                         gene," and ","immunomodulator",
                         "\nfrom the The Cancer Genome Atlas database")) %>% print()
  cor_wdata = spread(all_cancer_gene_cor[,c(1,3,6)],cancer,pearson,fill = NA)
  rownames(cor_wdata) = cor_wdata$gene2
  cor_wdata = cor_wdata[,-1]
  cor_wdata = cor_wdata[rownames(annotation_row),]
  
  cor_pdata = spread(all_cancer_gene_cor[,c(1,3,7)],cancer, pearson.p.value,fill = NA)
  rownames(cor_pdata) = cor_pdata$gene2
  cor_pdata = cor_pdata[,-1]
  display_numbers = matrix(ifelse(cor_pdata > 0.01 | is.na(cor_pdata), "×", ""), nrow(cor_pdata))
  cor_pdata = cor_pdata[rownames(annotation_row),]
  
  pheatmap(cor_wdata,
           #annotation_col =annotation_col,
           annotation_row = annotation_row,
           color =colorRampPalette(color.key)(50),
           cluster_cols =T,
           fontsize=6,
           cluster_rows = F,
           fontsize_row=6,
           cellwidth =10,
           cellheight =10,
           fontsize_number = 16,
           display_numbers= display_numbers,
           #scale="row",
           show_colnames=T,
           fontsize_col=8,
           main = paste0("Pearson correlation analysis of ",
                         gene," and ","immunomodulator",
                         "\nfrom the The Cancer Genome Atlas database")) %>% print()
}
dev.off()

