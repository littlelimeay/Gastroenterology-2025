export CONDA_ENVS_DIRS="/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/conda/envs/"
export CONDA_PKGS_DIRS="/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Active/conda/pkgs/"
export LSF_DOCKER_VOLUMES="/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/:/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/"
export PATH="/opt/conda/bin:$PATH"
#export PATH="/home/jichang/st1/anaconda3//bin:$PATH"

bsub -Is -G compute-gjrandolph -n 2 -q general-interactive -R 'rusage[mem=20GB]' -a 'docker(continuumio/anaconda3:2021.11)' /bin/bash

conda activate hdWGCNA






export LSF_DOCKER_VOLUMES='/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/:/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/'
bsub -Is -G compute-gjrandolph -n 2 -q general-interactive -R "select[hname!='compute1-exec-130.ris.wustl.edu']" -R 'rusage[mem=20GB]' -a 'docker(jichanghan/scrna_seurat5_monocle2:01292024)' /bin/bash

library(Seurat)
library(dplyr)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)   
library(sctransform)
library(harmony)
library(metap)
library(multtest)
library(biomaRt) 
library(monocle)
library(GSVA)
library(clustree)
library(rlang)
library(nichenetr)
library(tidyr)
library(tibble)

rm(list=ls()) 




#[1]Expression merge
quantile_normalization <- function(df){
  df_rank <- apply(df,2,rank,ties.method="min")
  df_sorted <- data.frame(apply(df, 2, sort))
  df_mean <- apply(df_sorted, 1, mean)
  
  index_to_mean <- function(my_index, my_mean){
    return(my_mean[my_index])
  }
  
  df_final <- apply(df_rank, 2, index_to_mean, my_mean=df_mean)
  rownames(df_final) <- rownames(df)
  return(df_final)
}

#source("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Software/Probe_data_affy.R")


myinf1 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/SR003137/download/all.gene_TPM.tsv"
data1<-read.table(myinf1, sep="\t",quote = "", header=T)
tag1<-grep('sample',colnames(data1))
data1<-unique(data1[,c(3,tag1)])
tag1<-duplicated(data1$external_gene_name)
data1<-data1[!tag1,]
row.names(data1)<-data1$external_gene_name
tag1<-grep('sample',colnames(data1))
data1<-data1[,tag1]



myinf2 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/SR004337/all.gene_TPM.tsv"
data2<-read.table(myinf2, sep="\t",quote = "", header=T)
tag2<-grep('sample',colnames(data2))
data2<-unique(data2[,c(3,tag2)])
tag2<-duplicated(data2$external_gene_name)
data2<-data2[!tag2,]
row.names(data2)<-data2$external_gene_name
tag2<-grep('sample',colnames(data2))
data2<-data2[,tag2]


tag<-colnames(data1)%in%colnames(data2)
data1<-data1[,!tag]


xx1 <- row.names(data1)
xx2 <- row.names(data2)

com <- Reduce(intersect, list(xx1,xx2))

data1<-data1[com,]
data2<-data2[com,]

data1<-log2(data1+1)
data2<-log2(data2+1)

data_log2 <- cbind(data1,data2)
tag<-apply(data_log2,1,function(x){sum(x)==0})
data_log2<-data_log2[!tag,]




tag1<-c('sample.L9','sample.L13','sample.L24','sample.L29','sample.L1','sample.L34')
tag2<-c('sample.L14','sample.L21','sample.L6','sample.L18','sample.L2','sample.L31','sample.L10','sample.L26','sample.L36')
tag3<-c('sample.L3','sample.L22','sample.L37','sample.L11','sample.L27','sample.L7','sample.L15','sample.L19','sample.L32')
tag4<-c('sample.A9','sample.A13','sample.A29')
tag5<-c('sample.A2','sample.A21','sample.A31','sample.A36')
tag6<-c('sample.A3','sample.A22','sample.A7','sample.A19')
tag7<-c('sample.D9','sample.D24','sample.D29')
tag8<-c('sample.D31','sample.D6','sample.D18','sample.D36')
tag9<-c('sample.D3','sample.D37','sample.D7','sample.D32')
colnames(data_log2)[colnames(data_log2)%in%tag1]<-paste0('Liver.Sham_',seq(1,length(tag1),1))
colnames(data_log2)[colnames(data_log2)%in%tag2]<-paste0('Liver.SBR.Vehicle_',seq(1,length(tag2),1))
colnames(data_log2)[colnames(data_log2)%in%tag3]<-paste0('Liver.SBR.WUSTL0717_',seq(1,length(tag3),1))
colnames(data_log2)[colnames(data_log2)%in%tag4]<-paste0('Ileum.Postana.Sham_',seq(1,length(tag4),1))
colnames(data_log2)[colnames(data_log2)%in%tag5]<-paste0('Ileum.Postana.SBR.Vehicle_',seq(1,length(tag5),1))
colnames(data_log2)[colnames(data_log2)%in%tag6]<-paste0('Ileum.Postana.SBR.WUSTL0717_',seq(1,length(tag6),1))
colnames(data_log2)[colnames(data_log2)%in%tag7]<-paste0('Duodenum.Sham_',seq(1,length(tag7),1))
colnames(data_log2)[colnames(data_log2)%in%tag8]<-paste0('Duodenum.SBR.Vehicle_',seq(1,length(tag8),1))
colnames(data_log2)[colnames(data_log2)%in%tag9]<-paste0('Duodenum.SBR.WUSTL0717_',seq(1,length(tag9),1))

myoutf1 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/Analysis_final_Exp1/With_all_liver_sample/Log2_metadata_before_combat_quantile.txt"
write.table(data_log2,myoutf1,sep="\t",quote=F)


myoutf <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/Analysis_final_Exp1/With_all_liver_sample/Normalization_check_before_combat_quantile.pdf"
pdf(myoutf,width=25,height=15)
boxplot(data_log2, las=2, main="Normalization check before combat quantile")
dev.off()



gene_means <- rowMeans(data_log2)
gene_variances <- apply(data_log2, 1, var)

# Step 2: Set thresholds to filter low-quality genes
# For example, keep genes with mean expression > 1 and variance > 0.5
mean_threshold <- 1
variance_threshold <- 0.5
filter_criteria <- (gene_means > mean_threshold) & (gene_variances > variance_threshold)

# Step 3: Apply the filter
data_log2 <- data_log2[filter_criteria, ]

myoutf <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/Analysis_final_Exp1/With_all_liver_sample/Normalization_check_before_combat_quantile_good_quality_genes.pdf"
pdf(myoutf,width=25,height=15)
boxplot(data_log2, las=2, main="Normalization check before combat quantile with only good quality genes")
dev.off()


sample_names <- colnames(data_log2)
batch <- ifelse(grepl("Liver", sample_names), "Liver", "Other")
batch <- as.factor(batch)
batch_df <- data.frame(batch = batch)
mod <- model.matrix(~1, data = batch_df)
data_log2 <- ComBat(dat = as.matrix(data_log2), batch = batch, mod = mod)
myoutf <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/Analysis_final_Exp1/With_all_liver_sample/Normalization_check_after_combat_good_quality_genes.pdf"
pdf(myoutf,width=25,height=15)
boxplot(data_log2, las=2, main="Batch Corrected Data with ComBat")
dev.off()


myoutf1 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/Analysis_final_Exp1/With_all_liver_sample/Log2_metadata_after_combat_no_quantile.txt"
write.table(data_log2,myoutf1,sep="\t",quote=F)





library(ggplot2)

tag<-apply(data_log2,1,function(x){sum(x)==0})
data_log2<-data_log2[!tag,]

res = prcomp(t(data_log2),scale. = TRUE)

df = res$x
df = df[,c("PC1","PC2")]
df = as.data.frame(df)
df$Group = unlist(strsplit(row.names(df), "_",fixed=T))[seq(1,92,2)]
#tag<-grep('Liver.SBR.WUSTL0717',df$Group)
#for (i in 1:length(tag))
#{
##  df$Group[tag[i]]<-paste0('Liver.SBR.WUSTL0717.',i)
#}

df$Group = factor(df$Group,levels = unique(df$Group))

myoutf = "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/Analysis_final_Exp1/With_all_liver_sample/PCA_plot_ggplot2_after_combat_all_samples.pdf"
pdf(myoutf,width=8,height=5)
p = ggplot(df,aes(x=PC1,y=PC2,color=Group))+ geom_point()
p
dev.off()










tag<-grep('Liver',colnames(data_log2))
data_log2_liver<-data_log2[,tag]
tag<-apply(data_log2_liver,1,function(x){sum(x)==0})
data_log2_liver<-data_log2_liver[!tag,]

res = prcomp(t(data_log2_liver),center=TRUE,scale. = TRUE)
df = res$x
df = df[,c("PC1","PC2")]
df = as.data.frame(df)
df$Group = unlist(strsplit(row.names(df), "_",fixed=T))[seq(1,48,2)]
df$Group = factor(df$Group,levels = unique(df$Group))

myoutf = "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/Analysis_final_Exp1/With_all_liver_sample/PCA_plot_liver_ggplot2_after_combatall_samples.pdf"
pdf(myoutf,width=8,height=5)
p = ggplot(df,aes(x=PC1,y=PC2,color=Group))+ geom_point()
p
dev.off()









tag<-grep('Liver',colnames(data_log2))
data_log2_liver<-data_log2[,tag]
tag<-apply(data_log2_liver,1,function(x){sum(x)==0})
data_log2_liver<-data_log2_liver[!tag,]

res = prcomp(t(data_log2_liver),center=TRUE,scale. = TRUE)
df = res$x
df = df[,c("PC1","PC2")]
df = as.data.frame(df)
df$Group = factor(row.names(df),levels = row.names(df))

myoutf = "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/Analysis_final_Exp1/With_all_liver_sample/PCA_plot_liver_ggplot2_after_combatall_samples_individual_sample.pdf"
pdf(myoutf,width=8,height=5)
p = ggplot(df,aes(x=PC1,y=PC2,color=Group))+ geom_point()
p
dev.off()












grep('Liver.Sham_6',colnames(data_log2))
data_log2<-data_log2[,1:45]

library(GSEABase)
myinf<-'/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/revision/Myeloid_PW_Adult_only_Human_nCountRNA_modified/PCx/res.0.6/GSVA/TF_activity/MSigDB/c3.tft.v2023.1.Hs.symbols.gmt'
database<-getGmt(myinf)

TF_use<-database[c(grep("^LXR", names(database)))]

LXR_Q3 <- geneIds(TF_use[[1]])
LXR_Q3_DR4 <- geneIds(TF_use[[2]])
LXR_overlap<-LXR_Q3[LXR_Q3%in%LXR_Q3_DR4]
LXRalpha <- c("ABCA1", "ABCD2", "ABCG1", "AEBP1", "BARD1", "BRCA1", "CETP", "CRP", "DHCR24", "ENG", "FABP4", "FOXO1", "HTT", 
"IFNG", "IGF1", "IL10", "JAK1", "LIPG", "LPCAT3", "LPL", "MIR1-1", "MIR206", "MIR613", "NFKB1", "NR0B2", "NR1H2", "NR1H3", "PLIN2",
 "PPARA", "PPARG", "PPARGC1B", "PRKACA", "PRMT3", "PTX3", "RORA", "RORC", "RXRA", "SCARB1", "SCD", "SREBF1", "STAT1", "THRB", "UGT1A1")


LXR_taget<-list(LXR_Q3,LXR_Q3_DR4,LXR_overlap,LXRalpha)
names(LXR_taget)<-c('LXR_Q3','LXR_Q3_DR4','LXR_overlap','LXRalpha')
#screen -r 238076

row.names(data_log2)<-toupper(row.names(data_log2))
gsva.es <- gsva(as.matrix(data_log2), LXR_taget, verbose=FALSE)

#gsva.es<-gsva.es[,c(1:8,10:34)]
#gsva.es <- apply(gsva.es,1, function(arg) (arg-mean(arg))/sd(arg))

myoutf = "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/Analysis_final_Exp1/With_all_liver_sample/GSVA_LXR_scale_with_liver_sample_No-scale.xls"

 write.table(gsva.es,myoutf,sep="\t",quote=F)





tag<-c(grep('Duodenum',colnames(data_log2)),grep('Ileum',colnames(data_log2)))
data_log2_intestine<-data_log2[,tag]

gsva.es_intestine <- gsva(as.matrix(data_log2_intestine), LXR_taget, verbose=FALSE)
#gsva.es_intestine <- apply(gsva.es_intestine,1, function(arg) (arg-mean(arg))/sd(arg))


myoutf = "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/Analysis_final_Exp1/With_all_liver_sample/GSVA_LXR_scale_only_intestine_no_scale.xls"
 write.table(gsva.es_intestine,myoutf,sep="\t",quote=F)






################Exclude the outliers and perform downstream analysis#######################
###########################################################################################
######3@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@2##################################
tar<-c('Liver.SBR.WUSTL0717_6','Liver.SBR.WUSTL0717_7','Liver.SBR.WUSTL0717_8','Liver.SBR.WUSTL0717_9','Liver.SBR.Vehicle_9','Liver.Sham_6')
data_log2<-data_log2[,!colnames(data_log2)%in%tar]

myoutf1 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/Analysis_final_Exp1/Exclude_outliers/Log2_metadata_after_combat_no_quantile.txt"
write.table(data_log2,myoutf1,sep="\t",quote=F)

myinf1 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/Analysis_final_Exp1/Exclude_outliers/Log2_metadata_after_combat_no_quantile.txt"
data_log2 <- read.table(myinf1,sep="\t",quote=NULL)

######GSVA for LXR analysis no liver outlier########

library(GSEABase)
myinf<-'/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/revision/Myeloid_PW_Adult_only_Human_nCountRNA_modified/PCx/res.0.6/GSVA/TF_activity/MSigDB/c3.tft.v2023.1.Hs.symbols.gmt'
database<-getGmt(myinf)

TF_use<-database[c(grep("^LXR", names(database)))]

LXR_Q3 <- geneIds(TF_use[[1]])
LXR_Q3_DR4 <- geneIds(TF_use[[2]])
LXR_overlap<-LXR_Q3[LXR_Q3%in%LXR_Q3_DR4]
LXRalpha <- c("ABCA1", "ABCD2", "ABCG1", "AEBP1", "BARD1", "BRCA1", "CETP", "CRP", "DHCR24", "ENG", "FABP4", "FOXO1", "HTT", 
"IFNG", "IGF1", "IL10", "JAK1", "LIPG", "LPCAT3", "LPL", "MIR1-1", "MIR206", "MIR613", "NFKB1", "NR0B2", "NR1H2", "NR1H3", "PLIN2",
 "PPARA", "PPARG", "PPARGC1B", "PRKACA", "PRMT3", "PTX3", "RORA", "RORC", "RXRA", "SCARB1", "SCD", "SREBF1", "STAT1", "THRB", "UGT1A1")


LXR_taget<-list(LXR_Q3,LXR_Q3_DR4,LXR_overlap,LXRalpha)
names(LXR_taget)<-c('LXR_Q3','LXR_Q3_DR4','LXR_overlap','LXRalpha')
#screen -r 238076

row.names(data_log2)<-toupper(row.names(data_log2))
gsva.es <- gsva(as.matrix(data_log2), LXR_taget, verbose=FALSE)

#gsva.es <- apply(gsva.es,1, function(arg) (arg-mean(arg))/sd(arg))


myoutf = "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/Analysis_final_Exp1/Exclude_outliers/GSVA_LXR_scale_with_liver_sample_no_scale.xls"

 write.table(gsva.es,myoutf,sep="\t",quote=F)












#1.  Liver  analysis

tag<-grep('Liver',colnames(data_log2))
data_log2_liver<-data_log2[,tag]


myoutf1 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/Analysis_final_Exp1/Exclude_outliers/Liver_TPMnormalized_log2.txt"
write.table(data_log2_liver,myoutf1,sep="\t",quote=F)

data_log2_liver_scaled <- apply(data_log2_liver, 1, function(x) {(x - mean(x)) / sd(x)})
data_log2_liver_scaled <- t(data_log2_liver_scaled)

myoutf1 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/Analysis_final_Exp1/Exclude_outliers/Liver_TPMnormalized_log2_scaled.txt"
write.table(data_log2_liver_scaled,myoutf1,sep="\t",quote=F)



zero_var_cols <- apply(data_log2_liver, 1, function(x) var(x) == 0)

data_log2_liver_filtered <- data_log2_liver[!zero_var_cols, ]
res = prcomp(t(data_log2_liver_filtered), center = TRUE, scale. = TRUE)

df = res$x
df = df[,c("PC1","PC2")]
df = as.data.frame(df)
df$Group = unlist(strsplit(row.names(df), "_",fixed=T))[seq(1,36,2)]
df$Group = factor(df$Group,levels = unique(df$Group))

myoutf = "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/Analysis_final_Exp1/Exclude_outliers/liver/PCA_plot_liver_ggplot2_after_combat_no_outlier_L34.pdf"
pdf(myoutf,width=8,height=5)
p = ggplot(df,aes(x=PC1,y=PC2,color=Group))+ geom_point()
p
dev.off()

#PCA plot with circles

myoutf = "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/Analysis_final_Exp1/Exclude_outliers/liver/PCA_plot_liver_ggplot2_after_combat_no_outlier_L34_with_circles.pdf"
pdf(myoutf,width=8,height=5)
p <- ggplot(df, aes(x = PC1, y = PC2, color = Group)) +
     geom_point(size = 3) +  # Add points for each sample
     stat_ellipse(aes(fill = Group), type = "t", linetype = 1, alpha = 0.5, 
                  level = 0.8) +  # Add ellipses with less overlap
     scale_color_manual(values = c("#1b9e77", "#d95f02", "#7570b3")) +  # Professional colors
     scale_fill_manual(values = c("#1b9e77", "#d95f02", "#7570b3")) +   # Matching fill colors
     labs(title = "PCA Plot", x = "PC1", y = "PC2") +  # Add labels
     theme_minimal() +  # Use a minimal theme for better appearance
     theme(legend.position = "right")
# Show the plot
print(p)
dev.off()



#2. GSVA for ECM and Fibrosis pathways


library(GSEABase)
myinf<-'/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/PM_KA_manuscript/Human_bulkRNA/analysis/GSEA/c2.all.v6.2.symbols_modified.gmt.txt'
database<-getGmt(myinf)

TF_use_collagen<-database[c(grep("COLLAGEN", names(database)))]
TF_use_ECM<-database[c(grep("ECM", names(database)))]
TF_use_Fibro<-database[c(grep("FIBROBLAST", names(database)))]


TF_use_combined <- c(TF_use_collagen, TF_use_ECM, TF_use_Fibro)
TF_use_combined_gsc <- GeneSetCollection(TF_use_combined)


row.names(data_log2_liver_filtered)<-toupper(row.names(data_log2_liver_filtered))
gsva.es <- gsva(as.matrix(data_log2_liver_filtered), TF_use_combined_gsc, verbose=FALSE)

#gsva.es<-gsva.es[,c(1:8,10:34)]
gsva.es <- apply(gsva.es,1, function(arg) (arg-mean(arg))/sd(arg))
myoutf = "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/Analysis_final_Exp1/Exclude_outliers/liver/GSVA_ECM_Collagen_Fibrosis.xls"
 write.table(gsva.es,myoutf,sep="\t",quote=F)




# 3. Calculate DEGs and plot volcano plot for certain genes


myinf1 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/Analysis_final_Exp1/Exclude_outliers/Log2_metadata_after_combat_no_quantile.txt"
data_log2 <- read.table(myinf1,sep="\t",quote=NULL)

tag<-grep('Liver',colnames(data_log2))
data_log2_liver<-data_log2[,tag]

zero_var_cols <- apply(data_log2_liver, 1, function(x) var(x) == 0)
data_log2_liver_filtered <- data_log2_liver[!zero_var_cols, ]
library(EnhancedVolcano)

res = prcomp(t(data_log2_liver_filtered), center = TRUE, scale. = TRUE)
df = res$x
df = df[,c("PC1","PC2")]
df = as.data.frame(df)
df$Group = unlist(strsplit(row.names(df), "_", fixed=T))[seq(1,36,2)]
df$Group = factor(df$Group, levels = unique(df$Group))

group <- unique(as.character(df$Group))
all_comb <- combn(group, 2)

for (j in 1:ncol(all_comb)) {
  
  res = matrix(0, nrow(data_log2_liver_filtered), 6)
  colnames(res) = c("Avg_xx1","Avg_xx2","Log2FC","T.pvalue","T.score","W.pval")
  row.names(res) = row.names(data_log2_liver_filtered)
  res = as.data.frame(res)
  
  for(i in 1:nrow(data_log2_liver_filtered)) {
    cat("\r", i)
    tag1 <- grep(all_comb[1, j], colnames(data_log2_liver_filtered))
    tag2 <- grep(all_comb[2, j], colnames(data_log2_liver_filtered))
    xx1 = data_log2_liver_filtered[i, tag1]
    xx2 = data_log2_liver_filtered[i, tag2]
    
    xx1 = as.numeric(xx1)
    xx2 = as.numeric(xx2)
    
    Avg_xx1 = mean(xx1)
    Avg_xx2 = mean(xx2)
    Log2FC = Avg_xx1 - Avg_xx2
    
    # Perform t-test and Wilcoxon test with error handling
    Fit1 <- try(t.test(xx1, xx2), silent = TRUE)
    Fit2 <- try(wilcox.test(xx1, xx2), silent = TRUE)
    
    if (inherits(Fit1, "try-error")) {
      T.pvalue <- NA
      T.score <- NA
    } else {
      T.pvalue <- Fit1$p.value
      T.score <- Fit1$statistic
    }
    
    if (inherits(Fit2, "try-error")) {
      W.pvalue <- NA
    } else {
      W.pvalue <- Fit2$p.value
    }
    
    res[i, "Avg_xx1"] = Avg_xx1
    res[i, "Avg_xx2"] = Avg_xx2
    res[i, "Log2FC"] = Log2FC
    res[i, "T.pvalue"] = T.pvalue
    res[i, "T.score"] = T.score
    res[i, "W.pval"] = W.pvalue
  }
  
  res <- res[order(-res$Log2FC),]
  res[,"T.pval.adj"] = p.adjust(res[,"T.pvalue"], method = "BH")
  
  myoutf = paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/Analysis_final_Exp1/Exclude_outliers/liver/DEG/", all_comb[1,j], "_VS_", all_comb[2,j], "_DEG_All.xls")
  write.table(res, myoutf, sep="\t", quote=F)
  
  tag <- res$T.pval.adj < 0.05
  res_sig <- res[tag,]
  myoutf = paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/Analysis_final_Exp1/Exclude_outliers/liver/DEG/", all_comb[1,j], "_VS_", all_comb[2,j], "_DEG_Sig.xls")
  write.table(res_sig, myoutf, sep="\t", quote=F)
  
  myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/Analysis_final_Exp1/Exclude_outliers/liver/DEG/", all_comb[1,j], "_VS_", all_comb[2,j], "_Volcano.pdf")
  pdf(myoutf, width=20, height=20)
  EnhancedVolcano(res,
    lab = row.names(res),
    x = 'Log2FC',
    y = 'T.pval.adj',
    xlim = c(min(res[['Log2FC']], na.rm = TRUE) - 0.5, max(res[['Log2FC']], na.rm = TRUE) + 0.5),
    ylim = c(0, max(-log10(res[['T.pval.adj']]), na.rm = TRUE) + 1),
    xlab = bquote(~Log[2]~ 'fold change'),
    title = paste0(all_comb[1,j], "_VS_", all_comb[2,j]),
    pCutoff = 5e-2,
    FCcutoff = 0.5,
    pointSize = 4.0,
    boxedLabels = F,
    labSize = 6.0,
    labCol = 'black',
    labFace = 'bold',
    col = c('black', 'black', 'black', 'red3'),
    colAlpha = 1,
    legendPosition = 'right',
    legendLabSize = 14,
    legendIconSize = 4.0,
    drawConnectors = TRUE,
    widthConnectors = 1.0,
    colConnectors = 'black')
  dev.off()
  
}





































library(EnhancedVolcano)

# Assume data_log2_liver_filtered is your preprocessed data matrix

# Identify the columns for the two groups
tag2 <- grep("Liver.SBR.Vehicle", colnames(data_log2_liver_filtered))
tag1 <- grep("Liver.SBR.WUSTL0717", colnames(data_log2_liver_filtered))

# Initialize the results matrix
res <- matrix(0, nrow(data_log2_liver_filtered), 6)
colnames(res) <- c("Avg_xx1", "Avg_xx2", "Log2FC", "T.pvalue", "T.score", "W.pval")
row.names(res) <- row.names(data_log2_liver_filtered)
res <- as.data.frame(res)

# Loop through each gene to calculate differential expression
for (i in 1:nrow(data_log2_liver_filtered)) {
  cat("\r", i)
  xx1 <- data_log2_liver_filtered[i, tag1]
  xx2 <- data_log2_liver_filtered[i, tag2]
  
  xx1 <- as.numeric(xx1)
  xx2 <- as.numeric(xx2)
  
  Avg_xx1 <- mean(xx1)
  Avg_xx2 <- mean(xx2)
  Log2FC <- Avg_xx1 - Avg_xx2
  
  Fit1 <- try(t.test(xx1, xx2), silent = TRUE)
  Fit2 <- try(wilcox.test(xx1, xx2), silent = TRUE)
  
  if (inherits(Fit1, "try-error")) {
    T.pvalue <- NA
    T.score <- NA
  } else {
    T.pvalue <- Fit1$p.value
    T.score <- Fit1$statistic
  }
  
  if (inherits(Fit2, "try-error")) {
    W.pvalue <- NA
  } else {
    W.pvalue <- Fit2$p.value
  }
  
  res[i, "Avg_xx1"] <- Avg_xx1
  res[i, "Avg_xx2"] <- Avg_xx2
  res[i, "Log2FC"] <- Log2FC
  res[i, "T.pvalue"] <- T.pvalue
  res[i, "T.score"] <- T.score
  res[i, "W.pval"] <- W.pvalue
}

# Adjust p-values for multiple testing
res <- res[order(-res$Log2FC),]
res[,"T.pval.adj"] <- p.adjust(res[,"T.pvalue"], method = "BH")

# Define the genes to label
genes_to_label <- c("Col1a1", "Col3a1", "Tnf", "Tgfb", "Pdgfa", "Pdgfb", "Cyp2c54", "Cyp4a12a", "Cd163", "Cyp2b10", "Traj14","Krt19","Krt7", "Cyp7a1", "Abca11",'Abca1','Abcb11', "Abcb4", "Slc51a", "Slc51b", "Fgf15")

# Create a column for labels, only labeling the specified genes
res$label <- ifelse(row.names(res) %in% genes_to_label, row.names(res), NA)

# Save the results
myoutf_all <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/Analysis_final_Exp1/Exclude_outliers/liver/DEG/SBR_vehicle_VS_SBR_WUSTL0717_DEG_All.xls"
write.table(res, myoutf_all, sep = "\t", quote = F)

# Filter for significant genes
res_sig <- res[res$T.pvalue < 0.05,]
myoutf_sig <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/Analysis_final_Exp1/Exclude_outliers/liver/DEG/SBR_vehicle_VS_SBR_WUSTL0717_DEG_Sig.xls"
write.table(res_sig, myoutf_sig, sep = "\t", quote = F)

# Generate the volcano plot

myoutf_volcano <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/Analysis_final_Exp1/Exclude_outliers/liver/DEG/SBR_vehicle_VS_SBR_WUSTL0717_Volcano_comprehensive_gene_labeled_1120.pdf"
pdf(myoutf_volcano, width = 15, height = 10)

EnhancedVolcano(res,
  lab = ifelse(res$T.pvalue < 0.05 & abs(res$Log2FC) > 0.1, res$label, NA),  # 只label通过cutoff的基因
  x = 'Log2FC',
  y = 'T.pvalue',
  xlim = c(min(res[['Log2FC']], na.rm = TRUE) - 0.5, max(res[['Log2FC']], na.rm = TRUE) + 0.5),
  ylim = c(0, max(-log10(res[['T.pvalue']]), na.rm = TRUE) + 1),
  xlab = bquote(~Log[2]~ 'fold change'),
  title = 'SBR WUSTL071 VS SBR vehicle',
  pCutoff = 0.05,
  FCcutoff = 0.1,
  pointSize = 2.0,
  boxedLabels = TRUE,  # 将标签框起来以提高清晰度
  labSize = 6.0,
  labCol = 'black',
  labFace = 'bold',
  col = c('lightgray', 'lightgray', 'lightgray', '#D55E00'),  # 修改颜色，第4种颜色为专业的红色
  colAlpha = 1,
  legendPosition = 'right',
  legendLabSize = 14,
  legendIconSize = 4.0,
  drawConnectors = TRUE,  # 连接标签和点
  widthConnectors = 1.0,
  colConnectors = 'black'
)

dev.off()

















library(EnhancedVolcano)

# Assume data_log2_liver_filtered is your preprocessed data matrix

# Identify the columns for the two groups
tag1 <- grep("Liver.SBR.Vehicle", colnames(data_log2_liver_filtered))
tag2 <- grep("Liver.SBR.WUSTL0717", colnames(data_log2_liver_filtered))

# Initialize the results matrix
res <- matrix(0, nrow(data_log2_liver_filtered), 6)
colnames(res) <- c("Avg_xx1", "Avg_xx2", "Log2FC", "T.pvalue", "T.score", "W.pval")
row.names(res) <- row.names(data_log2_liver_filtered)
res <- as.data.frame(res)

# Loop through each gene to calculate differential expression
for (i in 1:nrow(data_log2_liver_filtered)) {
  cat("\r", i)
  xx1 <- data_log2_liver_filtered[i, tag1]
  xx2 <- data_log2_liver_filtered[i, tag2]
  
  xx1 <- as.numeric(xx1)
  xx2 <- as.numeric(xx2)
  
  Avg_xx1 <- mean(xx1)
  Avg_xx2 <- mean(xx2)
  Log2FC <- Avg_xx1 - Avg_xx2
  
  Fit1 <- try(t.test(xx1, xx2), silent = TRUE)
  Fit2 <- try(wilcox.test(xx1, xx2), silent = TRUE)
  
  if (inherits(Fit1, "try-error")) {
    T.pvalue <- NA
    T.score <- NA
  } else {
    T.pvalue <- Fit1$p.value
    T.score <- Fit1$statistic
  }
  
  if (inherits(Fit2, "try-error")) {
    W.pvalue <- NA
  } else {
    W.pvalue <- Fit2$p.value
  }
  
  res[i, "Avg_xx1"] <- Avg_xx1
  res[i, "Avg_xx2"] <- Avg_xx2
  res[i, "Log2FC"] <- Log2FC
  res[i, "T.pvalue"] <- T.pvalue
  res[i, "T.score"] <- T.score
  res[i, "W.pval"] <- W.pvalue
}

# Adjust p-values for multiple testing
res <- res[order(-res$Log2FC),]
res[,"T.pval.adj"] <- p.adjust(res[,"T.pvalue"], method = "BH")

# Define the genes to label
genes_to_label <- c("Krt19", "Cyp7a1", "Abcb11", "Abcb4", "Slc51a", "Slc51b", "Fgf15")

# Create a column for labels, only labeling the specified genes
res$label <- ifelse(row.names(res) %in% genes_to_label, row.names(res), NA)

# Filter for significant genes
res_sig <- res[res$T.pvalue < 0.05,]
# Generate the volcano plot

myoutf_volcano <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/Analysis_final_Exp1/Exclude_outliers/liver/DEG/SBR_vehicle_VS_SBR_WUSTL0717_Volcano_specific_gene_labeled_krt19.pdf"
pdf(myoutf_volcano, width = 15, height = 10)

EnhancedVolcano(res,
  lab = ifelse(res$T.pvalue < 0.05 & abs(res$Log2FC) > 0.1, res$label, NA),  # 只label通过cutoff的基因
  x = 'Log2FC',
  y = 'T.pvalue',
  xlim = c(min(res[['Log2FC']], na.rm = TRUE) - 0.5, max(res[['Log2FC']], na.rm = TRUE) + 0.5),
  ylim = c(0, max(-log10(res[['T.pvalue']]), na.rm = TRUE) + 1),
  xlab = bquote(~Log[2]~ 'fold change'),
  title = 'SBR vehicle vs. SBR WUSTL0717',
  pCutoff = 0.05,
  FCcutoff = 0.1,
  pointSize = 2.0,
  boxedLabels = TRUE,  # 将标签框起来以提高清晰度
  labSize = 6.0,
  labCol = 'black',
  labFace = 'bold',
  col = c('lightgray', 'lightgray', 'lightgray', '#D55E00'),  # 修改颜色，第4种颜色为专业的红色
  colAlpha = 1,
  legendPosition = 'right',
  legendLabSize = 14,
  legendIconSize = 4.0,
  drawConnectors = TRUE,  # 连接标签和点
  widthConnectors = 1.0,
  colConnectors = 'black'
)

dev.off()


















library(EnhancedVolcano)

# Assume data_log2_liver_filtered is your preprocessed data matrix

# Identify the columns for the two groups
tag1 <- grep("Liver.SBR.Vehicle", colnames(data_log2_liver_filtered))
tag2 <- grep("Liver.SBR.WUSTL0717", colnames(data_log2_liver_filtered))

# Initialize the results matrix
res <- matrix(0, nrow(data_log2_liver_filtered), 6)
colnames(res) <- c("Avg_xx1", "Avg_xx2", "Log2FC", "T.pvalue", "T.score", "W.pval")
row.names(res) <- row.names(data_log2_liver_filtered)
res <- as.data.frame(res)

# Loop through each gene to calculate differential expression
for (i in 1:nrow(data_log2_liver_filtered)) {
  cat("\r", i)
  xx1 <- data_log2_liver_filtered[i, tag1]
  xx2 <- data_log2_liver_filtered[i, tag2]
  
  xx1 <- as.numeric(xx1)
  xx2 <- as.numeric(xx2)
  
  Avg_xx1 <- mean(xx1)
  Avg_xx2 <- mean(xx2)
  Log2FC <- Avg_xx1 - Avg_xx2
  
  Fit1 <- try(t.test(xx1, xx2), silent = TRUE)
  Fit2 <- try(wilcox.test(xx1, xx2), silent = TRUE)
  
  if (inherits(Fit1, "try-error")) {
    T.pvalue <- NA
    T.score <- NA
  } else {
    T.pvalue <- Fit1$p.value
    T.score <- Fit1$statistic
  }
  
  if (inherits(Fit2, "try-error")) {
    W.pvalue <- NA
  } else {
    W.pvalue <- Fit2$p.value
  }
  
  res[i, "Avg_xx1"] <- Avg_xx1
  res[i, "Avg_xx2"] <- Avg_xx2
  res[i, "Log2FC"] <- Log2FC
  res[i, "T.pvalue"] <- T.pvalue
  res[i, "T.score"] <- T.score
  res[i, "W.pval"] <- W.pvalue
}

# Adjust p-values for multiple testing
res <- res[order(-res$Log2FC),]
res[,"T.pval.adj"] <- p.adjust(res[,"T.pvalue"], method = "BH")

# Define the genes to label
genes_to_label <- c("Col1a1", "Col3a1", "Tnf", "Tgfb", "Pdgfa", "Pdgfb", "Cyp2c54", "Cyp4a12a", "Cd163", "Cyp2b10", "Traj14","Krt19", "Cyp7a1", "Abcb11", "Abcb4", "Slc51a", "Slc51b", "Fgf15")

# Create a column for labels, only labeling the specified genes
res$label <- ifelse(row.names(res) %in% genes_to_label, row.names(res), NA)

# Filter for significant genes
res_sig <- res[res$T.pvalue < 0.05,]
# Generate the volcano plot

myoutf_volcano <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/Analysis_final_Exp1/Exclude_outliers/liver/DEG/SBR_vehicle_VS_SBR_WUSTL0717_Volcano_specific_gene_labeled_krt19_Cola1.pdf"
pdf(myoutf_volcano, width = 15, height = 10)

EnhancedVolcano(res,
  lab = ifelse(res$T.pvalue < 0.05 & abs(res$Log2FC) > 0.1, res$label, NA),  # 只label通过cutoff的基因
  x = 'Log2FC',
  y = 'T.pvalue',
  xlim = c(min(res[['Log2FC']], na.rm = TRUE) - 0.5, max(res[['Log2FC']], na.rm = TRUE) + 0.5),
  ylim = c(0, max(-log10(res[['T.pvalue']]), na.rm = TRUE) + 1),
  xlab = bquote(~Log[2]~ 'fold change'),
  title = 'SBR vehicle vs. SBR WUSTL0717',
  pCutoff = 0.05,
  FCcutoff = 0.1,
  pointSize = 2.0,
  boxedLabels = TRUE,  # 将标签框起来以提高清晰度
  labSize = 6.0,
  labCol = 'black',
  labFace = 'bold',
  col = c('lightgray', 'lightgray', 'lightgray', '#D55E00'),  # 修改颜色，第4种颜色为专业的红色
  colAlpha = 1,
  legendPosition = 'right',
  legendLabSize = 14,
  legendIconSize = 4.0,
  drawConnectors = TRUE,  # 连接标签和点
  widthConnectors = 1.0,
  colConnectors = 'black'
)

dev.off()










# 定义基因组
genes_group1 <- c("Col1a1", "Col3a1", "Tgfb", "Pdgfa", "Pdgfb", "Tnf")
genes_group2 <- c("Krt19", "Cyp7a1", "Abcb11", "Abcb4", "Slc51a", "Slc51b", "Fgf15")

# 从数据中提取所需基因的表达值
data_group1 <- data_log2_liver_filtered[rownames(data_log2_liver_filtered) %in% genes_group1, ]
data_group2 <- data_log2_liver_filtered[rownames(data_log2_liver_filtered) %in% genes_group2, ]

# Z-score标准化
data_group1_z <- t(scale(t(data_group1)))
data_group2_z <- t(scale(t(data_group2)))

# 定义输出文件路径和文件名
output_folder <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/Analysis_final_Exp1/Exclude_outliers/liver/DEG/"
output_file1 <- paste0(output_folder, "Group1_Genes_Expression_Z_Score.xls")
output_file2 <- paste0(output_folder, "Group2_Genes_Expression_Z_Score.xls")

# 将数据写入文件
write.table(data_group1_z, file = output_file1, sep = "\t", quote = FALSE, col.names = NA)
write.table(data_group2_z, file = output_file2, sep = "\t", quote = FALSE, col.names = NA)







# 定义基因组
genes_group1 <- c("Col1a1", "Col3a1", "Tgfb", "Pdgfa", "Pdgfb", "Tnf")
genes_group2 <- c("Krt19", "Cyp7a1", "Abcb11", "Abcb4", "Slc51a", "Slc51b", "Fgf15")

# 提取所需基因的原始表达值
data_group1_raw <- data_log2_liver_filtered[rownames(data_log2_liver_filtered) %in% genes_group1, ]
data_group2_raw <- data_log2_liver_filtered[rownames(data_log2_liver_filtered) %in% genes_group2, ]

# 定义输出文件路径和文件名
output_folder <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/Analysis_final_Exp1/Exclude_outliers/liver/DEG/"
output_file1_raw <- paste0(output_folder, "Group1_Genes_Expression_Raw.xls")
output_file2_raw <- paste0(output_folder, "Group2_Genes_Expression_Raw.xls")

# 保存原始数据
write.table(data_group1_raw, file = output_file1_raw, sep = "\t", quote = FALSE, col.names = NA)
write.table(data_group2_raw, file = output_file2_raw, sep = "\t", quote = FALSE, col.names = NA)


























































library(ggplot2)

tag<-apply(data_log2,1,function(x){sum(x)==0})
data_log2<-data_log2[!tag,]

res = prcomp(t(data_log2),scale. = TRUE)

df = res$x
df = df[,c("PC1","PC2")]
df = as.data.frame(df)
df$Group = unlist(strsplit(row.names(df), "_",fixed=T))[seq(1,92,2)]
#tag<-grep('Liver.SBR.WUSTL0717',df$Group)
#for (i in 1:length(tag))
#{
##  df$Group[tag[i]]<-paste0('Liver.SBR.WUSTL0717.',i)
#}

df$Group = factor(df$Group,levels = unique(df$Group))

myoutf = "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/Analysis_final_Exp1/Exclude_outliers/PCA_plot_ggplot2_no_outlier.pdf"
pdf(myoutf,width=8,height=5)
p = ggplot(df,aes(x=PC1,y=PC2,color=Group))+ geom_point()
p
dev.off()






myoutf2 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/Analysis_final_Exp1/Exclude_outliers/Log2_metadata_after_combat_quantile.txt"

data_log2 <- quantile_normalization(data_log2)
batch <- c(rep(1,ncol(data1)),rep(2,ncol(data2)))
library(sva)
metadata <- ComBat(dat=as.matrix(data_log2),batch=batch,mod=NULL,par.prior=T)
write.table(metadata,myoutf2,sep="\t",quote=F)


library(ggplot2)
res = prcomp(t(metadata),scale. = TRUE)

df = res$x
df = df[,c("PC1","PC2")]
df = as.data.frame(df)
df$Group = unlist(strsplit(row.names(df), "_",fixed=T))[seq(1,92,2)]


myoutf = "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/Analysis_final_Exp1/Exclude_outliers/PCA_plot_ggplot2_after_combat.pdf"
pdf(myoutf,width=8,height=5)
p = ggplot(df,aes(x=PC1,y=PC2,color=Group))+ geom_point()
p
dev.off()










########Generate PCA plot for all liver samples only##########





##################Liver.SBR.WUSTL0717.3 is the outlier sample###########################

################Generate GSEA for each two groups comparision with or without liver SBR WUSTL 3#####################
#1. with all samples

all_sample_group<-unique(unlist(strsplit(colnames(data_log2), "_",fixed=T))[seq(1,68,2)])

all_comb<-combn(all_sample_group,2)

for (j in 1:ncol(all_comb))
{

res = matrix(0, nrow(data_log2),6)
colnames(res) = c("Avg_xx1","Avg_xx2","Log2FC","T.pvalue","T.score","W.pval")
row.names(res) = row.names(data_log2)
res = as.data.frame(res)

  for(i in 1 : nrow(data_log2))
   {
	cat("\r",i)
   tag1<-grep(all_comb[1,j],colnames(data_log2))
   tag2<-grep(all_comb[2,j],colnames(data_log2))
	xx1 = data_log2[i,tag1]
	xx2 = data_log2[i,tag2]
	
	xx1 = as.numeric(xx1)
	xx2 = as.numeric(xx2)
	
	Avg_xx1 = mean(xx1)
	Avg_xx2 = mean(xx2)
	Log2FC = Avg_xx1 - Avg_xx2
	
	Fit1 = t.test(xx1,xx2)
	Fit2 = wilcox.test(xx1,xx2)
	
	T.pvalue = Fit1$p.value
	T.score = Fit1$statistic
	W.pvalue = Fit2$p.value
	
	res[i,"Avg_xx1"] = Avg_xx1
	res[i,"Avg_xx2"] = Avg_xx2
	res[i,"Log2FC"] = Log2FC
	res[i,"T.pvalue"] = T.pvalue
	res[i,"T.score"] = T.score
	res[i,"W.pval"] = W.pvalue

  }
  
  
  res<-res[order(-res$Log2FC),]
  res[,"T.pval.adj"] = p.adjust(res[,"T.pvalue"] ,method="BH")

 myoutf = paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/Analysis_final_Exp1/Exclude_outliers/GSEA/marker_all_sample/",all_comb[1,j],"_VS_",all_comb[2,j],"_DEG_All.xls")
 write.table(res,myoutf,sep="\t",quote=F)

}







#2. without outlier
data_log2<-data_log2[,c(1:8,10:34)]

all_sample_group<-unique(unlist(strsplit(colnames(data_log2), "_",fixed=T))[seq(1,66,2)])

all_comb<-combn(all_sample_group,2)

for (j in 1:ncol(all_comb))
{

res = matrix(0, nrow(data_log2),6)
colnames(res) = c("Avg_xx1","Avg_xx2","Log2FC","T.pvalue","T.score","W.pval")
row.names(res) = row.names(data_log2)
res = as.data.frame(res)

  for(i in 1 : nrow(data_log2))
   {
	cat("\r",i)
   tag1<-grep(all_comb[1,j],colnames(data_log2))
   tag2<-grep(all_comb[2,j],colnames(data_log2))
	xx1 = data_log2[i,tag1]
	xx2 = data_log2[i,tag2]
	
	xx1 = as.numeric(xx1)
	xx2 = as.numeric(xx2)
	
	Avg_xx1 = mean(xx1)
	Avg_xx2 = mean(xx2)
	Log2FC = Avg_xx1 - Avg_xx2
	
	Fit1 = t.test(xx1,xx2)
	Fit2 = wilcox.test(xx1,xx2)
	
	T.pvalue = Fit1$p.value
	T.score = Fit1$statistic
	W.pvalue = Fit2$p.value
	
	res[i,"Avg_xx1"] = Avg_xx1
	res[i,"Avg_xx2"] = Avg_xx2
	res[i,"Log2FC"] = Log2FC
	res[i,"T.pvalue"] = T.pvalue
	res[i,"T.score"] = T.score
	res[i,"W.pval"] = W.pvalue

  }
  
  
  res<-res[order(-res$Log2FC),]
  res[,"T.pval.adj"] = p.adjust(res[,"T.pvalue"] ,method="BH")

 myoutf = paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/Analysis_final_Exp1/Exclude_outliers/GSEA/marker_without_outlier/",all_comb[1,j],"_VS_",all_comb[2,j],"_DEG_All.xls")
 write.table(res,myoutf,sep="\t",quote=F)

}





dir <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/Analysis_final_Exp1/Exclude_outliers/GSEA/marker_without_outlier/"
files <- list.files(dir)
file_nam <- gsub(".xls","",files)
for(i in 1:length(files))
{
  cat("\r",i)
  tmpinf <- paste0(dir,files[i])
  res <- read.table(tmpinf,sep="\t",quote=NULL)
  tag<-res$Log2FC==0
  res<-res[!tag,]
  res[,'ranking_matric']<-res$Log2FC*(-log10(res$T.pvalue))
  res <- res[order(res$ranking_matric,decreasing=T),]
  
  xx <- res$ranking_matric
  names(xx) <- row.names(res)
  names(xx) <- toupper(names(xx))
  
  tag = duplicated(names(xx))
  res = res[tag==0,]
  row.names(res) <- toupper(row.names(res))
  res[,"Name"] <- row.names(res)
  res[,"metric"] <- res$ranking_matric
  res <- res[,c("Name","metric")]
  
  myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/Analysis_final_Exp1/Exclude_outliers/GSEA/preranked_GSEA_no_outlier/",file_nam[i],".rnk")
  write.table(res,file=myoutf,quote=F,sep="\t",row.names=F)
}






dir <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/Analysis_final_Exp1/Exclude_outliers/GSEA/marker_without_outlier/"
files <- list.files(dir)
file_nam <- gsub(".xls","",files)
for(i in 1:length(files))
{
  cat("\r",i)
  tmpinf <- paste0(dir,files[i])
  res <- read.table(tmpinf,sep="\t",quote=NULL)
  tag<-res$Log2FC==0
  res<-res[!tag,]
  res[,'ranking_matric']<-res$Log2FC*(-log10(res$T.pvalue))*(-1)
  res <- res[order(res$ranking_matric,decreasing=T),]
  
  xx <- res$ranking_matric
  names(xx) <- row.names(res)
  names(xx) <- toupper(names(xx))
  
  tag = duplicated(names(xx))
  res = res[tag==0,]
  row.names(res) <- toupper(row.names(res))
  res[,"Name"] <- row.names(res)
  res[,"metric"] <- res$ranking_matric
  res <- res[,c("Name","metric")]
  
  myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/Analysis_final_Exp1/Exclude_outliers/GSEA/preranked_GSEA_no_outlier_reversed/",unlist(strsplit(files[i], "_",fixed=T))[3],
                   '_VS_',unlist(strsplit(files[i], "_",fixed=T))[1],".rnk")
  write.table(res,file=myoutf,quote=F,sep="\t",row.names=F)
}







#####Try quantile nomrlaization and combat?##############
#######Do not use#################3
data_log2_quan<-quantile_normalization(data_log2)

res = prcomp(t(data_log2_quan),scale. = TRUE)

df = res$x
df = df[,c("PC1","PC2")]
df = as.data.frame(df)
df$Group = unlist(strsplit(row.names(df), "_",fixed=T))[seq(1,68,2)]

df$Group = factor(df$Group,levels = unique(df$Group))

myoutf = "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/Analysis_final_Exp1/Exclude_outliers/PCA_plot_ggplot2_quantile.pdf"
pdf(myoutf,width=8,height=5)
p = ggplot(df,aes(x=PC1,y=PC2,color=Group))+ geom_point()
p
dev.off()



batch<-c(1,1,1,1,2,2,3,3,3,4,4,4,5,5,5,6,6,7,7,7,8,9,9,2,2,3,5,6,6,8,8,8,9,9)

library(sva)
metadata3 <- ComBat(dat=as.matrix(data_log2_quan),batch=batch,mod=NULL,par.prior=T)

res = prcomp(t(metadata3),scale. = TRUE)

df = res$x
df = df[,c("PC1","PC2")]
df = as.data.frame(df)
df$Group = unlist(strsplit(row.names(df), "_",fixed=T))[seq(1,68,2)]

df$Group = factor(df$Group,levels = unique(df$Group))

myoutf = "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/Analysis_final_Exp1/Exclude_outliers/PCA_plot_ggplot2_quantile_combat.pdf"
pdf(myoutf,width=8,height=5)
p = ggplot(df,aes(x=PC1,y=PC2,color=Group))+ geom_point()
p
dev.off()




#########@@@@@@@@@@@@@@@LIVER Analysis########@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
######PCA only on liver#########
tag<-grep('Liver',colnames(data_log2))
data_log2_liver<-data_log2[,tag]

tag<-apply(data_log2_liver,1,function(x){sum(x)==0})
data_log2_liver<-data_log2_liver[!tag,]


gene_means <- rowMeans(data_log2_liver)
gene_variances <- apply(data_log2_liver, 1, var)

# Step 2: Set thresholds to filter low-quality genes
# For example, keep genes with mean expression > 1 and variance > 0.5
mean_threshold <- 0.2
variance_threshold <- 0.2
filter_criteria <- (gene_means > mean_threshold) & (gene_variances > variance_threshold)

# Step 3: Apply the filter
data_log2_liver <- data_log2_liver[filter_criteria, ]




res = prcomp(t(data_log2_liver),center=TRUE,scale. = TRUE)
df = res$x
df = df[,c("PC1","PC2")]
df = as.data.frame(df)
df$Group = unlist(strsplit(row.names(df), "_",fixed=T))[seq(1,36,2)]
df$Group = factor(df$Group,levels = unique(df$Group))

myoutf = "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/Analysis_final_Exp1/Exclude_outliers/liver/PCA_plot_liver_ggplot2_before_combat_no_outlier_L34.pdf"
pdf(myoutf,width=8,height=5)
p = ggplot(df,aes(x=PC1,y=PC2,color=Group))+ geom_point()
p
dev.off()







sample_names <- colnames(data_log2_liver)
batch <- ifelse(grepl("Sham", sample_names), "Sham", 
                ifelse(grepl("Vehicle", sample_names), "Vehicle", "WUSTL0717"))

min_samples <- ncol(data_log2_liver) / 2
average_threshold <- 1
filter_criteria <- rowSums(data_log2_liver > 1) >= min_samples & rowMeans(data_log2_liver) > average_threshold
filtered_log2_TPM <- data_log2_liver[filter_criteria, ]
mod <- model.matrix(~1, data=as.data.frame(colnames(filtered_log2_TPM)))
combat_data <- ComBat(dat = filtered_log2_TPM, batch = batch, mod = mod)

res = prcomp(t(combat_data), center = TRUE, scale. = TRUE)

df = res$x
df = df[,c("PC1","PC2")]
df = as.data.frame(df)
df$Group = unlist(strsplit(row.names(df), "_",fixed=T))[seq(1,48,2)]


myoutf = "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/Analysis_final_Exp1/Exclude_outliers/liver/PCA_plot_liver_ggplot2_after_combat.pdf"
pdf(myoutf,width=8,height=5)
p = ggplot(df,aes(x=PC1,y=PC2,color=Group))+ geom_point()
p
dev.off()












res = prcomp(t(data_log2_liver),center=TRUE,scale. = TRUE)
df = res$x
df = df[,c("PC1","PC2")]
df = as.data.frame(df)
df$Group = unlist(strsplit(row.names(df), "_",fixed=T))[seq(1,48,2)]
df$Group = factor(df$Group,levels = unique(df$Group))

myoutf = "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/Analysis_final_Exp1/Exclude_outliers/liver/PCA_plot_liver_ggplot2_after_combat_after_quantile.pdf"
pdf(myoutf,width=8,height=5)
p = ggplot(df,aes(x=PC1,y=PC2,color=Group))+ geom_point()
p
dev.off()







##########DEG and Volcano plot###############
#######Generate gene list for GSEA analysis########### 


library(EnhancedVolcano)


group<-unique(as.character(df$Group))
all_comb<-combn(group,2)
for (j in 1:ncol(all_comb))
{

res = matrix(0, nrow(data_log2_liver),6)
colnames(res) = c("Avg_xx1","Avg_xx2","Log2FC","T.pvalue","T.score","W.pval")
row.names(res) = row.names(data_log2_liver)
res = as.data.frame(res)

  for(i in 1 : nrow(data_log2_liver))
   {
	cat("\r",i)
   tag1<-grep(all_comb[1,j],colnames(data_log2_liver))
   tag2<-grep(all_comb[2,j],colnames(data_log2_liver))
	xx1 = data_log2_liver[i,tag1]
	xx2 = data_log2_liver[i,tag2]
	
	xx1 = as.numeric(xx1)
	xx2 = as.numeric(xx2)
	
	Avg_xx1 = mean(xx1)
	Avg_xx2 = mean(xx2)
	Log2FC = Avg_xx1 - Avg_xx2
	
	Fit1 = t.test(xx1,xx2)
	Fit2 = wilcox.test(xx1,xx2)
	
	T.pvalue = Fit1$p.value
	T.score = Fit1$statistic
	W.pvalue = Fit2$p.value
	
	res[i,"Avg_xx1"] = Avg_xx1
	res[i,"Avg_xx2"] = Avg_xx2
	res[i,"Log2FC"] = Log2FC
	res[i,"T.pvalue"] = T.pvalue
	res[i,"T.score"] = T.score
	res[i,"W.pval"] = W.pvalue

  }
  
  
  res<-res[order(-res$Log2FC),]
  res[,"T.pval.adj"] = p.adjust(res[,"T.pvalue"] ,method="BH")

 myoutf = paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/SR003137/download/JH_processed/liver/DEG/",all_comb[1,j],"_VS_",all_comb[2,j],"_DEG_All.xls")
 write.table(res,myoutf,sep="\t",quote=F)
 tag<-res$T.pval.adj<0.05
  res_sig<-res[tag,]
 myoutf = paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/SR003137/download/JH_processed/liver/DEG/",all_comb[1,j],"_VS_",all_comb[2,j],"_DEG_Sig.xls")
 write.table(res_sig,myoutf,sep="\t",quote=F)


 myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/SR003137/download/JH_processed/liver/DEG/",all_comb[1,j],"_VS_",all_comb[2,j],"_Volcano.pdf")
 pdf(myoutf,width=20,height=20)
  EnhancedVolcano(res,
    lab = row.names(res),
    x = 'Log2FC',
    y = 'T.pval.adj',
	 xlim = c(min(res[['Log2FC']], na.rm = TRUE) - 0.5, max(res[['Log2FC']], na.rm = TRUE) +
    0.5),
  ylim = c(0, max(-log10(res[['T.pval.adj']]), na.rm = TRUE) + 1),
    xlab = bquote(~Log[2]~ 'fold change'),
    title = paste0(all_comb[1,j],"_VS_",all_comb[2,j]),
    pCutoff = 5e-2,
    FCcutoff = 0.5,
    pointSize = 4.0,
	 boxedLabels = F,
	  labSize = 6.0,
    labCol = 'black',
    labFace = 'bold',
    col=c('black', 'black', 'black', 'red3'),
    colAlpha = 1,
	 legendPosition = 'right',
    legendLabSize = 14,
    legendIconSize = 4.0,
    drawConnectors = TRUE,
    widthConnectors = 1.0,
    colConnectors = 'black')
 dev.off()
 
}


























#####2. Ileum Analysis Oct 2024      Ileum.Postana region

#(1)##########Generate vocalno plot for SBR vs SBR WUSTL for Ileum.Postana region###############

myinf1 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/Analysis_final_Exp1/Exclude_outliers/Log2_metadata_after_combat_no_quantile.txt"
data_log2 <- read.table(myinf1,sep="\t",quote=NULL)


tag<-grep('Ileum.Postana',colnames(data_log2))
data_log2_ileum<-data_log2[,tag]


myoutf1 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/Analysis_final_Exp1/Exclude_outliers/Ileum_TPMnormalized_log2.txt"
write.table(data_log2_ileum,myoutf1,sep="\t",quote=F)

data_log2_ileum_scaled <- apply(data_log2_ileum, 1, function(x) {(x - mean(x)) / sd(x)})
data_log2_ileum_scaled <- t(data_log2_ileum_scaled)

myoutf1 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/Analysis_final_Exp1/Exclude_outliers/Ileum_TPMnormalized_log2_scaled.txt"
write.table(data_log2_ileum_scaled,myoutf1,sep="\t",quote=F)




library(EnhancedVolcano)

tag<-apply(data_log2_ileum,1,function(x){sum(x)==0})
data_log2_ileum<-data_log2_ileum[!tag,]

zero_var_cols <- apply(data_log2_ileum, 1, function(x) var(x) == 0)
data_log2_ileum_filtered <- data_log2_ileum[!zero_var_cols, ]


res = prcomp(t(data_log2_ileum_filtered), center = TRUE, scale. = TRUE)
df = res$x
df = df[,c("PC1","PC2")]
df = as.data.frame(df)
df$Group = unlist(strsplit(row.names(df), "_",fixed=T))[seq(1,22,2)]



# Define the genes to label
genes_to_label <- c('Srebf1', 'Abca1', 'Scd1','Ppara','Lpcat3','Abcg8','Abcg5','Abcg1')


group<-unique(as.character(df$Group))
all_comb<-combn(group,2)
for (j in 1:ncol(all_comb))
{

res = matrix(0, nrow(data_log2_ileum_filtered),6)
colnames(res) = c("Avg_xx1","Avg_xx2","Log2FC","T.pvalue","T.score","W.pval")
row.names(res) = row.names(data_log2_ileum_filtered)
res = as.data.frame(res)
for(i in 1 : nrow(data_log2_ileum_filtered)) {
    cat("\r", i)
    
    tag1 <- grep(all_comb[2,j], colnames(data_log2_ileum_filtered))
    tag2 <- grep(all_comb[1,j], colnames(data_log2_ileum_filtered))
    
    xx1 = data_log2_ileum_filtered[i, tag1]
    xx2 = data_log2_ileum_filtered[i, tag2]
    
    xx1 = as.numeric(xx1)
    xx2 = as.numeric(xx2)
    
    # Skip this gene if the variance is too small (near constant values)
    if (var(xx1) < 1e-6 || var(xx2) < 1e-6) {
        next
    }
    
    Avg_xx1 = mean(xx1)
    Avg_xx2 = mean(xx2)
    Log2FC = Avg_xx1 - Avg_xx2
    
    Fit1 = t.test(xx1, xx2)
    Fit2 = wilcox.test(xx1, xx2)
    
    T.pvalue = Fit1$p.value
    T.score = Fit1$statistic
    W.pvalue = Fit2$p.value
    
    res[i,"Avg_xx1"] = Avg_xx1
    res[i,"Avg_xx2"] = Avg_xx2
    res[i,"Log2FC"] = Log2FC
    res[i,"T.pvalue"] = T.pvalue
    res[i,"T.score"] = T.score
    res[i,"W.pval"] = W.pvalue
}
  
  res<-res[order(-res$Log2FC),]
  res[,"T.pval.adj"] = p.adjust(res[,"T.pvalue"] ,method="BH")

 myoutf = paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/Analysis_final_Exp1/Exclude_outliers/ileum/DEG/",all_comb[2,j],"_VS_",all_comb[1,j],"_DEG_All.xls")
 write.table(res,myoutf,sep="\t",quote=F)
 res$label <- ifelse(row.names(res) %in% genes_to_label, row.names(res), NA)
 res<-res[res$T.pvalue != 0,]
 res_sig <- res[res$T.pvalue < 0.05,]
 res_sig <- res_sig[abs(res_sig$Log2FC) > 0,]
 

 myoutf = paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/Analysis_final_Exp1/Exclude_outliers/ileum/DEG/",all_comb[2,j],"_VS_",all_comb[1,j],"_DEG_Sig.xls")
 write.table(res_sig,myoutf,sep="\t",quote=F)


 myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/Analysis_final_Exp1/Exclude_outliers/ileum/Volcano_plots/",all_comb[2,j],"_VS_",all_comb[1,j],"_Volcano_p005.pdf")
 pdf(myoutf,width=15,height=13)
 
 EnhancedVolcano(res,
  lab = ifelse(res$T.pvalue < 0.05 & abs(res$Log2FC) > 0.1, res$label, NA),  # 只label通过cutoff的基因
  x = 'Log2FC',
  y = 'T.pvalue',
  xlim = c(min(res[['Log2FC']], na.rm = TRUE) - 0.5, max(res[['Log2FC']], na.rm = TRUE) + 0.5),
  ylim = c(0, max(-log10(res[['T.pvalue']]), na.rm = TRUE) + 1),
  xlab = bquote(~Log[2]~ 'fold change'),
  title = 'SBR WUSTL0717 VS SBR vehicle',
  pCutoff = 0.05,
  FCcutoff = 0.1,
  pointSize = 2.0,
  boxedLabels = TRUE, 
  labSize = 6.0,
  labCol = 'black',
  labFace = 'bold',
  col = c('lightgray', 'lightgray', 'lightgray', '#D55E00'),  # 修改颜色，第4种颜色为专业的红色
  colAlpha = 1,
  legendPosition = 'right',
  legendLabSize = 14,
  legendIconSize = 4.0,
  drawConnectors = TRUE,  # 连接标签和点
  widthConnectors = 1.0,
  colConnectors = 'black'
)
    dev.off()
 }



#(2) GSEA analyisis

myinf1 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/Analysis_final_Exp1/Exclude_outliers/Log2_metadata_after_combat_no_quantile.txt"
data_log2 <- read.table(myinf1,sep="\t",quote=NULL)
tag<-grep('Ileum.Postana',colnames(data_log2))
data_log2_ileum<-data_log2[,tag]

tag<-apply(data_log2_ileum,1,function(x){sum(x)==0})
data_log2_ileum<-data_log2_ileum[!tag,]
zero_var_cols <- apply(data_log2_ileum, 1, function(x) var(x) == 0)
data_log2_ileum_filtered <- data_log2_ileum[!zero_var_cols, ]

res = prcomp(t(data_log2_ileum_filtered), center = TRUE, scale. = TRUE)
df = res$x
df = df[,c("PC1","PC2")]
df = as.data.frame(df)
df$Group = unlist(strsplit(row.names(df), "_",fixed=T))[seq(1,22,2)]


group<-unique(as.character(df$Group))
all_comb<-combn(group,2)
for (j in 1:ncol(all_comb))
{

res = matrix(0, nrow(data_log2_ileum_filtered),6)
colnames(res) = c("Avg_xx1","Avg_xx2","Log2FC","T.pvalue","T.score","W.pval")
row.names(res) = row.names(data_log2_ileum_filtered)
res = as.data.frame(res)
for(i in 1 : nrow(data_log2_ileum_filtered)) {
    cat("\r", i)
    
    tag1 <- grep(all_comb[1,j], colnames(data_log2_ileum_filtered))
    tag2 <- grep(all_comb[2,j], colnames(data_log2_ileum_filtered))
    
    xx1 = data_log2_ileum_filtered[i, tag1]
    xx2 = data_log2_ileum_filtered[i, tag2]
    
    xx1 = as.numeric(xx1)
    xx2 = as.numeric(xx2)
    
    # Skip this gene if the variance is too small (near constant values)
    if (var(xx1) < 1e-6 || var(xx2) < 1e-6) {
        next
    }
    
    Avg_xx1 = mean(xx1)
    Avg_xx2 = mean(xx2)
    Log2FC = Avg_xx1 - Avg_xx2
    
    Fit1 = t.test(xx1, xx2)
    Fit2 = wilcox.test(xx1, xx2)
    
    T.pvalue = Fit1$p.value
    T.score = Fit1$statistic
    W.pvalue = Fit2$p.value
    
    res[i,"Avg_xx1"] = Avg_xx1
    res[i,"Avg_xx2"] = Avg_xx2
    res[i,"Log2FC"] = Log2FC
    res[i,"T.pvalue"] = T.pvalue
    res[i,"T.score"] = T.score
    res[i,"W.pval"] = W.pvalue
}
  
  res<-res[order(-res$Log2FC),]
  res[,"T.pval.adj"] = p.adjust(res[,"T.pvalue"] ,method="BH")
 
 res<-res[abs(res$Log2FC) > 0,]
 res<-res[res$T.pvalue != 0,]

 myoutf = paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/Analysis_final_Exp1/Exclude_outliers/ileum/GSEA/marker/",all_comb[1,j],"_VS_",all_comb[2,j],"_DEG_All.xls")
 write.table(res,myoutf,sep="\t",quote=F)
}



dir <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/Analysis_final_Exp1/Exclude_outliers/ileum/GSEA/marker/"
files <- list.files(dir)
file_nam <- gsub(".xls","",files)
for(i in 1:length(files))
{
  cat("\r",i)
  tmpinf <- paste0(dir,files[i])
  res <- read.table(tmpinf,sep="\t",quote=NULL)
  tag<-res$Log2FC==0
  res<-res[!tag,]
  res[,'ranking_matric']<-res$Log2FC*(-log10(res$T.pvalue))
  res <- res[order(res$ranking_matric,decreasing=T),]
  
  xx <- res$ranking_matric
  names(xx) <- row.names(res)
  names(xx) <- toupper(names(xx))
  
  tag = duplicated(names(xx))
  res = res[tag==0,]
  row.names(res) <- toupper(row.names(res))
  res[,"Name"] <- row.names(res)
  res[,"metric"] <- res$ranking_matric
  res <- res[,c("Name","metric")]
  
  myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/Analysis_final_Exp1/Exclude_outliers/ileum/GSEA/prerank/weighted/",file_nam[i],".rnk")
  write.table(res,file=myoutf,quote=F,sep="\t",row.names=F)
}






dir <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/Analysis_final_Exp1/Exclude_outliers/ileum/GSEA/marker/"
files <- list.files(dir)
file_nam <- gsub(".xls","",files)
for(i in 1:length(files))
{
  cat("\r",i)
  tmpinf <- paste0(dir,files[i])
  res <- read.table(tmpinf,sep="\t",quote=NULL)
  tag<-res$Log2FC==0
  res<-res[!tag,]
  res[,'ranking_matric']<-res$Log2FC
  res <- res[order(res$ranking_matric,decreasing=T),]
  
  xx <- res$ranking_matric
  names(xx) <- row.names(res)
  names(xx) <- toupper(names(xx))
  
  tag = duplicated(names(xx))
  res = res[tag==0,]
  row.names(res) <- toupper(row.names(res))
  res[,"Name"] <- row.names(res)
  res[,"metric"] <- res$ranking_matric
  res <- res[,c("Name","metric")]
  
  myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/Analysis_final_Exp1/Exclude_outliers/ileum/GSEA/prerank/non-weighted/",file_nam[i],".rnk")
  write.table(res,file=myoutf,quote=F,sep="\t",row.names=F)
}
































#######3. Sham Duodenum vs ileum, show different pathways including that nutrition is higher in duodenum and ileum has enriched pathways for anti-inflammation


tag<-grep('Duodenum',colnames(data_log2))
data_log2_Duodenum<-data_log2[,tag]


myoutf1 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/Analysis_final_Exp1/Exclude_outliers/Duodenum_TPMnormalized_log2.txt"
write.table(data_log2_Duodenum,myoutf1,sep="\t",quote=F)

data_log2_Duodenum_scaled <- apply(data_log2_Duodenum, 1, function(x) {(x - mean(x)) / sd(x)})
data_log2_Duodenum_scaled <- t(data_log2_Duodenum_scaled)

myoutf1 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/Analysis_final_Exp1/Exclude_outliers/Duodenum_TPMnormalized_log2_scaled.txt"
write.table(data_log2_Duodenum_scaled,myoutf1,sep="\t",quote=F)




myinf1 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/Analysis_final_Exp1/Exclude_outliers/Log2_metadata_after_combat_no_quantile.txt"
data_log2 <- read.table(myinf1,sep="\t",quote=NULL)
tag1<-grep('Duodenum.Sham',colnames(data_log2))
tag2<-grep('Ileum.Postana.Sham',colnames(data_log2))

data_log2_ileum<-data_log2[,tag2]
data_log2_Duodenum<-data_log2[,tag1]
data_log2_ileum_vs_duo_sham<-cbind(data_log2_ileum,data_log2_Duodenum)

tag<-apply(data_log2_ileum_vs_duo_sham,1,function(x){sum(x)==0})
data_log2_ileum_vs_duo_sham<-data_log2_ileum_vs_duo_sham[!tag,]
zero_var_cols <- apply(data_log2_ileum_vs_duo_sham, 1, function(x) var(x) == 0)
data_log2_ileum_vs_duo_sham_filtered <- data_log2_ileum_vs_duo_sham[!zero_var_cols, ]

res = prcomp(t(data_log2_ileum_vs_duo_sham_filtered), center = TRUE, scale. = TRUE)
df = res$x
df = df[,c("PC1","PC2")]
df = as.data.frame(df)
df$Group = unlist(strsplit(row.names(df), "_",fixed=T))[seq(1,12,2)]


group<-unique(as.character(df$Group))
all_comb<-combn(group,2)
for (j in 1:ncol(all_comb))
{

res = matrix(0, nrow(data_log2_ileum_vs_duo_sham_filtered),6)
colnames(res) = c("Avg_xx1","Avg_xx2","Log2FC","T.pvalue","T.score","W.pval")
row.names(res) = row.names(data_log2_ileum_vs_duo_sham_filtered)
res = as.data.frame(res)
for(i in 1 : nrow(data_log2_ileum_vs_duo_sham_filtered)) {
    cat("\r", i)
    
    tag1 <- grep(all_comb[1,j], colnames(data_log2_ileum_vs_duo_sham_filtered))
    tag2 <- grep(all_comb[2,j], colnames(data_log2_ileum_vs_duo_sham_filtered))
    
    xx1 = data_log2_ileum_vs_duo_sham_filtered[i, tag1]
    xx2 = data_log2_ileum_vs_duo_sham_filtered[i, tag2]
    
    xx1 = as.numeric(xx1)
    xx2 = as.numeric(xx2)
    
    # Skip this gene if the variance is too small (near constant values)
    if (var(xx1) < 1e-6 || var(xx2) < 1e-6) {
        next
    }
    
    Avg_xx1 = mean(xx1)
    Avg_xx2 = mean(xx2)
    Log2FC = Avg_xx1 - Avg_xx2
    
    Fit1 = t.test(xx1, xx2)
    Fit2 = wilcox.test(xx1, xx2)
    
    T.pvalue = Fit1$p.value
    T.score = Fit1$statistic
    W.pvalue = Fit2$p.value
    
    res[i,"Avg_xx1"] = Avg_xx1
    res[i,"Avg_xx2"] = Avg_xx2
    res[i,"Log2FC"] = Log2FC
    res[i,"T.pvalue"] = T.pvalue
    res[i,"T.score"] = T.score
    res[i,"W.pval"] = W.pvalue
}
  
  res<-res[order(-res$Log2FC),]
  res[,"T.pval.adj"] = p.adjust(res[,"T.pvalue"] ,method="BH")
 
 res<-res[abs(res$Log2FC) > 0,]
 res<-res[res$T.pvalue != 0,]

 myoutf = paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/Analysis_final_Exp1/Exclude_outliers/ileum_vs_duodenum/GSEA/marker/",all_comb[1,j],"_VS_",all_comb[2,j],"_DEG_All.xls")
 write.table(res,myoutf,sep="\t",quote=F)
}



dir <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/Analysis_final_Exp1/Exclude_outliers/ileum_vs_duodenum/GSEA/marker/"
files <- list.files(dir)
file_nam <- gsub(".xls","",files)
for(i in 1:length(files))
{
  cat("\r",i)
  tmpinf <- paste0(dir,files[i])
  res <- read.table(tmpinf,sep="\t",quote=NULL)
  tag<-res$Log2FC==0
  res<-res[!tag,]
  res[,'ranking_matric']<-res$Log2FC*(-log10(res$T.pvalue))
  res <- res[order(res$ranking_matric,decreasing=T),]
  
  xx <- res$ranking_matric
  names(xx) <- row.names(res)
  names(xx) <- toupper(names(xx))
  
  tag = duplicated(names(xx))
  res = res[tag==0,]
  row.names(res) <- toupper(row.names(res))
  res[,"Name"] <- row.names(res)
  res[,"metric"] <- res$ranking_matric
  res <- res[,c("Name","metric")]
  
  myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/Analysis_final_Exp1/Exclude_outliers/ileum_vs_duodenum/GSEA/prerank/weighted/",file_nam[i],".rnk")
  write.table(res,file=myoutf,quote=F,sep="\t",row.names=F)
}






dir <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/Analysis_final_Exp1/Exclude_outliers/ileum_vs_duodenum/GSEA/marker/"
files <- list.files(dir)
file_nam <- gsub(".xls","",files)
for(i in 1:length(files))
{
  cat("\r",i)
  tmpinf <- paste0(dir,files[i])
  res <- read.table(tmpinf,sep="\t",quote=NULL)
  tag<-res$Log2FC==0
  res<-res[!tag,]
  res[,'ranking_matric']<-res$Log2FC
  res <- res[order(res$ranking_matric,decreasing=T),]
  
  xx <- res$ranking_matric
  names(xx) <- row.names(res)
  names(xx) <- toupper(names(xx))
  
  tag = duplicated(names(xx))
  res = res[tag==0,]
  row.names(res) <- toupper(row.names(res))
  res[,"Name"] <- row.names(res)
  res[,"metric"] <- res$ranking_matric
  res <- res[,c("Name","metric")]
  
  myoutf <- paste0("/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/Analysis_final_Exp1/Exclude_outliers/ileum_vs_duodenum/GSEA/prerank/non-weighted/",file_nam[i],".rnk")
  write.table(res,file=myoutf,quote=F,sep="\t",row.names=F)
}



























































######################Heatmap#############################
#####1. LXR targets#################Use this on10/25/2024 for Ayoung data

data_log2 <- cbind(data1,data2)
tag<-apply(data_log2,1,function(x){sum(x)==0})
data_log2<-data_log2[!tag,]

tag1<-c('sample.L9','sample.L13','sample.L24','sample.L29','sample.L1','sample.L34')
tag2<-c('sample.L14','sample.L21','sample.L6','sample.L18','sample.L2','sample.L31','sample.L10','sample.L26','sample.L36')
tag3<-c('sample.L3','sample.L22','sample.L37','sample.L11','sample.L27','sample.L7','sample.L15','sample.L19','sample.L32')
tag4<-c('sample.A9','sample.A13','sample.A29')
tag5<-c('sample.A2','sample.A21','sample.A31','sample.A36')
tag6<-c('sample.A3','sample.A22','sample.A7','sample.A19')
tag7<-c('sample.D9','sample.D24','sample.D29')
tag8<-c('sample.D31','sample.D6','sample.D18','sample.D36')
tag9<-c('sample.D3','sample.D37','sample.D7','sample.D32')
colnames(data_log2)[colnames(data_log2)%in%tag1]<-paste0('Liver.Sham_',seq(1,length(tag1),1))
colnames(data_log2)[colnames(data_log2)%in%tag2]<-paste0('Liver.SBR.Vehicle_',seq(1,length(tag2),1))
colnames(data_log2)[colnames(data_log2)%in%tag3]<-paste0('Liver.SBR.WUSTL0717_',seq(1,length(tag3),1))
colnames(data_log2)[colnames(data_log2)%in%tag4]<-paste0('Ileum.Postana.Sham_',seq(1,length(tag4),1))
colnames(data_log2)[colnames(data_log2)%in%tag5]<-paste0('Ileum.Postana.SBR.Vehicle_',seq(1,length(tag5),1))
colnames(data_log2)[colnames(data_log2)%in%tag6]<-paste0('Ileum.Postana.SBR.WUSTL0717_',seq(1,length(tag6),1))
colnames(data_log2)[colnames(data_log2)%in%tag7]<-paste0('Duodenum.Sham_',seq(1,length(tag7),1))
colnames(data_log2)[colnames(data_log2)%in%tag8]<-paste0('Duodenum.SBR.Vehicle_',seq(1,length(tag8),1))
colnames(data_log2)[colnames(data_log2)%in%tag9]<-paste0('Duodenum.SBR.WUSTL0717_',seq(1,length(tag9),1))



tar<-c('Liver.SBR.WUSTL0717_6','Liver.SBR.WUSTL0717_7','Liver.SBR.WUSTL0717_8','Liver.SBR.WUSTL0717_9','Liver.SBR.Vehicle_9','Liver.Sham_6')
data_log2<-data_log2[,!colnames(data_log2)%in%tar]





gene_means <- rowMeans(data_log2)
gene_variances <- apply(data_log2, 1, var)

# Step 2: Set thresholds to filter low-quality genes
# For example, keep genes with mean expression > 1 and variance > 0.5
mean_threshold <- 1
variance_threshold <- 0.5
filter_criteria <- (gene_means > mean_threshold) & (gene_variances > variance_threshold)

# Step 3: Apply the filter
data_log2 <- data_log2[filter_criteria, ]

sample_names <- colnames(data_log2)
batch <- ifelse(grepl("Liver", sample_names), "Liver", "Other")
batch <- as.factor(batch)
batch_df <- data.frame(batch = batch)
mod <- model.matrix(~1, data = batch_df)
data_log2 <- ComBat(dat = as.matrix(data_log2), batch = batch, mod = mod)
myoutf <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/Analysis_final_Exp1/Exclude_outliers/Combat_batch_correct_check_after_removing_outliers.pdf"
pdf(myoutf,width=25,height=15)
boxplot(data_log2, las=2, main="Batch Corrected Data with ComBat")
dev.off()





#data_log2 <- normalizeBetweenArrays(data_log2, method="quantile")
#batch <- as.factor(batch)
#batch_df <- data.frame(batch = batch)
#mod <- model.matrix(~1, data = batch_df)
#data_log2 <- ComBat(dat = as.matrix(data_log2), batch = batch, mod = mod)
#myoutf <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/Analysis_final_Exp1/Exclude_outliers/Combat_batch_correct_check_after_quantile.pdf"
#pdf(myoutf,width=25,height=15)
#boxplot(data_log2, las=2, main="Batch Corrected Data with ComBat")
#dev.off()





library(GSEABase)
myinf<-'/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/revision/Myeloid_PW_Adult_only_Human_nCountRNA_modified/PCx/res.0.6/GSVA/TF_activity/MSigDB/c3.tft.v2023.1.Hs.symbols.gmt'
database<-getGmt(myinf)

TF_use<-database[c(grep("^LXR", names(database)))]

LXR_Q3 <- geneIds(TF_use[[1]])
LXR_Q3_DR4 <- geneIds(TF_use[[2]])
LXR_overlap<-LXR_Q3[LXR_Q3%in%LXR_Q3_DR4]
LXRalpha <- c("ABCA1", "ABCD2", "ABCG1", "AEBP1", "BARD1", "BRCA1", "CETP", "CRP", "DHCR24", "ENG", "FABP4", "FOXO1", "HTT", 
"IFNG", "IGF1", "IL10", "JAK1", "LIPG", "LPCAT3", "LPL", "MIR1-1", "MIR206", "MIR613", "NFKB1", "NR0B2", "NR1H2", "NR1H3", "PLIN2",
 "PPARA", "PPARG", "PPARGC1B", "PRKACA", "PRMT3", "PTX3", "RORA", "RORC", "RXRA", "SCARB1", "SCD", "SREBF1", "STAT1", "THRB", "UGT1A1")


LXR_taget<-list(LXR_Q3,LXR_Q3_DR4,LXR_overlap,LXRalpha)
names(LXR_taget)<-c('LXR_Q3','LXR_Q3_DR4','LXR_overlap','LXRalpha')
#screen -r 238076

row.names(data_log2)<-toupper(row.names(data_log2))
gsva.es <- gsva(as.matrix(data_log2), LXR_taget, verbose=FALSE)

#gsva.es<-gsva.es[,c(1:8,10:34)]
gsva.es <- apply(gsva.es,1, function(arg) (arg-mean(arg))/sd(arg))
myoutf = "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/Analysis_final_Exp1/Exclude_outliers/GSVA_LXR_1.xls"
 write.table(gsva.es,myoutf,sep="\t",quote=F)







data_use<-data_log2[LXR_overlap,]
samples <- c(
  "Liver.Sham_1",
  "Liver.Sham_2",
  "Liver.Sham_3",
  "Liver.Sham_4",
  "Liver.SBR.Vehicle_1",
  "Liver.SBR.Vehicle_2",
  "Liver.SBR.Vehicle_3",
  "Liver.SBR.Vehicle_4",
  "Liver.SBR.WUSTL0717_1",
  "Liver.SBR.WUSTL0717_2",
  "Liver.SBR.WUSTL0717_4",
  "Duodenum.Sham_1",
  "Duodenum.Sham_2",
  "Duodenum.Sham_3",
  "Duodenum.SBR.Vehicle_1",
  "Duodenum.SBR.Vehicle_2",
  "Duodenum.SBR.Vehicle_3",
  "Duodenum.SBR.Vehicle_4",
  "Duodenum.SBR.WUSTL0717_1",
  "Duodenum.SBR.WUSTL0717_2",
  "Duodenum.SBR.WUSTL0717_3",
  "Duodenum.SBR.WUSTL0717_4",
  "Ileum.Postana.Sham_1",
  "Ileum.Postana.Sham_2",
  "Ileum.Postana.Sham_3",
  "Ileum.Postana.SBR.Vehicle_1",
  "Ileum.Postana.SBR.Vehicle_2",
  "Ileum.Postana.SBR.Vehicle_3",
  "Ileum.Postana.SBR.Vehicle_4",
  "Ileum.Postana.SBR.WUSTL0717_1",
  "Ileum.Postana.SBR.WUSTL0717_2",
  "Ileum.Postana.SBR.WUSTL0717_3",
  "Ileum.Postana.SBR.WUSTL0717_4"
)
data_use<-data_use[,samples]
data_use_liver<-data_use[,grep('Liver',colnames(data_use))]
data_use_duo<-data_use[,grep('Duodenum',colnames(data_use))]
data_use_ileum<-data_use[,grep('Ileum',colnames(data_use))]

data_use_intestine<-cbind(data_use_duo,data_use_ileum)


data_use_liver <- apply(data_use_liver,1, function(arg) (arg-mean(arg))/sd(arg))
data_use_duo <- apply(data_use_duo,1, function(arg) (arg-mean(arg))/sd(arg))
data_use_ileum <- apply(data_use_ileum,1, function(arg) (arg-mean(arg))/sd(arg))
data_use_intestine<-apply(data_use_intestine,1, function(arg) (arg-mean(arg))/sd(arg))

res_liver <- t(data_use_liver)
res_duo <- t(data_use_duo)
res_ileum <- t(data_use_ileum)
data_use_intestine<-t(data_use_intestine)
#res<-res[,c("2",'5','4','8','9','0','1','7','3','6')]

#ts<-res[sort(row.names(res)),]
breaklist_1<-seq(-2,0,0.1)
breaklist_2<-seq(0,2,0.1)
breaklist<-unique(c(breaklist_1,breaklist_2))

myoutf <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/Analysis_final_Exp1/Exclude_outliers/heatmap/LXR_overlap_GSEAdataset_liver.pdf"
pdf(myoutf,width=25,height=15)
pheatmap(mat=res_liver,
         color = colorRampPalette(rev(brewer.pal(n = 9, name ="RdYlBu")))(length(breaklist)),
         breaks = breaklist,
         cluster_rows=FALSE, 
         cluster_cols=FALSE,
         border_color      = NA,
         show_colnames     = T,
         show_rownames     = T,
         drop_levels       = TRUE,
         scale='none',
         
         fontsize=20)
dev.off()



myoutf <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/Analysis_final_Exp1/Exclude_outliers/heatmap/LXR_overlap_GSEAdataset_duor.pdf"
pdf(myoutf,width=25,height=15)
pheatmap(mat=res_duo,
         color = colorRampPalette(rev(brewer.pal(n = 9, name ="RdYlBu")))(length(breaklist)),
         breaks = breaklist,
         cluster_rows=FALSE, 
         cluster_cols=FALSE,
         border_color      = NA,
         show_colnames     = T,
         show_rownames     = T,
         drop_levels       = TRUE,
         scale='none',
         
         fontsize=20)
dev.off()


myoutf <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/Analysis_final_Exp1/Exclude_outliers/heatmap/LXR_overlap_GSEAdataset_ileum.pdf"
pdf(myoutf,width=10,height=15)
pheatmap(mat=res_ileum,
         color = colorRampPalette(rev(brewer.pal(n = 9, name ="RdYlBu")))(length(breaklist)),
         breaks = breaklist,
         cluster_rows=FALSE, 
         cluster_cols=FALSE,
         border_color      = NA,
         show_colnames     = T,
         show_rownames     = T,
         drop_levels       = TRUE,
         scale='none',
         
         fontsize=20)
dev.off()




myoutf <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/Analysis_final_Exp1/Exclude_outliers/heatmap/LXR_overlap_GSEAdataset_intestine.pdf"
pdf(myoutf,width=10,height=15)
pheatmap(mat=data_use_intestine,
         color = colorRampPalette(rev(brewer.pal(n = 9, name ="RdYlBu")))(length(breaklist)),
         breaks = breaklist,
         cluster_rows=FALSE, 
         cluster_cols=TRUE,
         border_color      = NA,
         show_colnames     = T,
         show_rownames     = T,
         drop_levels       = TRUE,
         scale='none',
         
         fontsize=20)
dev.off()








enriched <- unique( c("Chmp2b", "Tfeb", "Slc2a4", "Acaca", "Lrp3", "Syt17", "Rara", "Mid1ip1", "Lpcat3", "Sgk1",
                 "Adamts4", "Acsl3", "Abcg1", "Apoc1", "Abca1","Sod1", "Aldh2", "Dbi", "Pck1", "Hp", "Gng11", "Rbp4",
                  "Adrb3", "Cfb", "Apoe", "Nnat", "Cfd", "Fabp4", "Cyp2e1", "Thrsp", "Srebf1", "Scd1", "Apoc1", "Abca1","Cyp2c55", "Srebf1", "Scd3", "Abcg1", "Scd2", "Adipoq", "Ces2a",  "Cyp3a11", "Acaa1b", 
          "Cyp2e1", "G6pc", "Pnpla3", 
            "Fabp4", "Fabp1", "Cyp3a44",
           "Slc6a2",
           "Cyp2c66","Abcg5", "Vwa1", "Apoc3","Ap2a2", "Apoa2", "Pcsk6", "Alb", "Scarb1", "Pltp", "Clta", "Cubn", "Creb3l3", "Nr1h3", "Bmp1", "Amn", "Npc1", "Nceh1", "Apoe", "Apoa1", 
           "Lpl", "Apoa4", "Cidec", "Apoc3", "Mylip", "Abcg1", "Apoc1", "Abca1"))





data_use<-data_log2[enriched,]
samples <- c(
  "Liver.Sham_1",
  "Liver.Sham_2",
  "Liver.Sham_3",
  "Liver.Sham_4",
  "Liver.SBR.Vehicle_1",
  "Liver.SBR.Vehicle_2",
  "Liver.SBR.Vehicle_3",
  "Liver.SBR.Vehicle_4",
  "Liver.SBR.WUSTL0717_1",
  "Liver.SBR.WUSTL0717_2",
  "Liver.SBR.WUSTL0717_4",
  "Duodenum.Sham_1",
  "Duodenum.Sham_2",
  "Duodenum.Sham_3",
  "Duodenum.SBR.Vehicle_1",
  "Duodenum.SBR.Vehicle_2",
  "Duodenum.SBR.Vehicle_3",
  "Duodenum.SBR.Vehicle_4",
  "Duodenum.SBR.WUSTL0717_1",
  "Duodenum.SBR.WUSTL0717_2",
  "Duodenum.SBR.WUSTL0717_3",
  "Duodenum.SBR.WUSTL0717_4",
  "Ileum.Postana.Sham_1",
  "Ileum.Postana.Sham_2",
  "Ileum.Postana.Sham_3",
  "Ileum.Postana.SBR.Vehicle_1",
  "Ileum.Postana.SBR.Vehicle_2",
  "Ileum.Postana.SBR.Vehicle_3",
  "Ileum.Postana.SBR.Vehicle_4",
  "Ileum.Postana.SBR.WUSTL0717_1",
  "Ileum.Postana.SBR.WUSTL0717_2",
  "Ileum.Postana.SBR.WUSTL0717_3",
  "Ileum.Postana.SBR.WUSTL0717_4"
)
data_use<-data_use[,samples]
data_use_liver<-data_use[,grep('Liver',colnames(data_use))]
data_use_duo<-data_use[,grep('Duodenum',colnames(data_use))]
data_use_ileum<-data_use[,grep('Ileum',colnames(data_use))]

data_use_intestine<-cbind(data_use_duo,data_use_ileum)


data_use_liver <- apply(data_use_liver,1, function(arg) (arg-mean(arg))/sd(arg))
data_use_duo <- apply(data_use_duo,1, function(arg) (arg-mean(arg))/sd(arg))
data_use_ileum <- apply(data_use_ileum,1, function(arg) (arg-mean(arg))/sd(arg))
data_use_intestine<-apply(data_use_intestine,1, function(arg) (arg-mean(arg))/sd(arg))

res_liver <- t(data_use_liver)
res_duo <- t(data_use_duo)
res_ileum <- t(data_use_ileum)
data_use_intestine<-t(data_use_intestine)
#res<-res[,c("2",'5','4','8','9','0','1','7','3','6')]

#ts<-res[sort(row.names(res)),]
breaklist_1<-seq(-2,0,0.1)
breaklist_2<-seq(0,2,0.1)
breaklist<-unique(c(breaklist_1,breaklist_2))


myoutf <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/Analysis_final_Exp1/Exclude_outliers/heatmap/LXR_GSEA_enriched_GSEAdataset_ileum.pdf"
pdf(myoutf,width=10,height=20)
pheatmap(mat=res_ileum,
         color = colorRampPalette(rev(brewer.pal(n = 9, name ="RdYlBu")))(length(breaklist)),
         breaks = breaklist,
         cluster_rows=FALSE, 
         cluster_cols=FALSE,
         border_color      = NA,
         show_colnames     = T,
         show_rownames     = T,
         drop_levels       = TRUE,
         scale='none',
         
         fontsize=20)
dev.off()








############Dot PLot to plot GSEA pathways#################

#(1) Show that ileum versus duodenum GSEA
myinf<-'/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/Analysis_final_Exp1/Exclude_outliers/ileum_vs_duodenum/GSEA/TSV/Ileum_vs_Duo_Sham_UP.tsv'
res <- read.table(myinf,sep="\t",quote=NULL,header=TRUE)
res[,'enriched_percentage']<-rep(0,nrow(res))

for (i in 1:nrow(res))
{
  xx<-unlist(strsplit(res$LEADING.EDGE[i], "=",fixed=T))[2]
  xx<-as.numeric(unlist(strsplit(xx, "%,",fixed=T))[1])
  res[,'enriched_percentage'][i]<-0.01*xx
}

res$number_genes_enriched<-round(res$SIZE*res$enriched_percentage,digits=0)




myinf<-'/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/Analysis_final_Exp1/Exclude_outliers/ileum_vs_duodenum/GSEA/TSV/Ileum_vs_Duo_Sham_DN.tsv'
res_1 <- read.table(myinf,sep="\t",quote=NULL,header=TRUE)
res_1[,'enriched_percentage']<-rep(0,nrow(res_1))

for (i in 1:nrow(res_1))
{
  xx<-unlist(strsplit(res_1$LEADING.EDGE[i], "=",fixed=T))[2]
  xx<-as.numeric(unlist(strsplit(xx, "%,",fixed=T))[1])
  res_1[,'enriched_percentage'][i]<-0.01*xx
}

res_1$number_genes_enriched<-round(res_1$SIZE*res_1$enriched_percentage,digits=0)



res_UP_DN<-rbind(res,res_1)






tag1<-grep('TGFBR',res_UP_DN$NAME)
tag2<-grep('TNF',res_UP_DN$NAME)
tag3<-grep('TOLL_RECEPTOR',res_UP_DN$NAME)
tag4<-grep('TLR4',res_UP_DN$NAME)
tag5<-grep('COSTIMULATION',res_UP_DN$NAME)
tag6<-grep('DENDRITIC',res_UP_DN$NAME)
tag7<-grep('CXCR4',res_UP_DN$NAME)
tag8<-grep('FC_GAMMA_R',res_UP_DN$NAME)
tag9<-grep('B_LYMPHOCYTE',res_UP_DN$NAME)
tag10<-grep('NFKAPPAB',res_UP_DN$NAME)
tag11<-grep('T_CELL',res_UP_DN$NAME)
tag12<-grep('INNATE_IMMUNE_SYSTEM',res_UP_DN$NAME)
tag13<-grep('TCR',res_UP_DN$NAME)
tag14<-grep('FOXP3',res_UP_DN$NAME)

tag15<-grep('TCA',res_UP_DN$NAME)
tag16<-grep('GLYCO',res_UP_DN$NAME)
tag17<-grep('FATTY_ACID',res_UP_DN$NAME)
tag18<-grep('ADIPO',res_UP_DN$NAME)
tag19<-grep('CHOLESTEROL',res_UP_DN$NAME)



tag<-c(tag1,tag2,tag3,tag4,tag5,tag6,tag7,tag8,tag9,tag10,tag11,tag12,tag13,tag14,tag15,tag16,tag17,tag18,tag19)
res_use<-res_UP_DN[tag,]
res_use<-res_use[res_use$FDR.q.val<0.1,]

res_use<-res_use[order(-res_use$NES),]
res_use$NAME<-factor(res_use$NAME,levels=unique(res_use$NAME))

myoutf<-'/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/Analysis_final_Exp1/Exclude_outliers/ileum_vs_duodenum/GSEA/Dotplot/Ileum_VS_Duo_Shan_UP_DN.pdf'
pdf(myoutf,width=15,height=15)
ggplot(res_use, aes(x = NES, y = NAME)) + 
               geom_point(aes(size = number_genes_enriched, color = FDR.q.val)) +
               theme_minimal(base_size = 16) +  # 使用 theme_minimal 去掉背景格子，并调大基础字号
               scale_size(range = c(4, 10)) +  # 增大 dot size 的范围
               scale_colour_gradient(limits=c(0, 0.10), low="red") +
               ylab(NULL) +
               ggtitle("GSEA") +
               geom_vline(xintercept = 0, linetype="dashed", color = "black") +  # 在 NES=0 位置添加一条竖直虚线
               theme(axis.text.y = element_text(size = 14),  # 调整y轴文字的大小
                     axis.text.x = element_text(size = 14),  # 调整x轴文字的大小
                     plot.title = element_text(size = 18, hjust = 0.5))

dev.off()
























#(2) GSEA Sham versus SBR 
myinf<-'/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/Analysis_final_Exp1/Exclude_outliers/ileum/GSEA/TSV/SBR_vehicle_VS_WUSTL0717/Ileum.Postana.Sham_VS_Ileum.Postana.SBR.Vehicle_UP.tsv'
res <- read.table(myinf,sep="\t",quote=NULL,header=TRUE)
res[,'enriched_percentage']<-rep(0,nrow(res))

for (i in 1:nrow(res))
{
  xx<-unlist(strsplit(res$LEADING.EDGE[i], "=",fixed=T))[2]
  xx<-as.numeric(unlist(strsplit(xx, "%,",fixed=T))[1])
  res[,'enriched_percentage'][i]<-0.01*xx
}

res$number_genes_enriched<-round(res$SIZE*res$enriched_percentage,digits=0)




myinf<-'/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/Analysis_final_Exp1/Exclude_outliers/ileum/GSEA/TSV/SBR_vehicle_VS_WUSTL0717/Ileum.Postana.Sham_VS_Ileum.Postana.SBR.Vehicle_DN.tsv'
res_1 <- read.table(myinf,sep="\t",quote=NULL,header=TRUE)
res_1[,'enriched_percentage']<-rep(0,nrow(res_1))

for (i in 1:nrow(res_1))
{
  xx<-unlist(strsplit(res_1$LEADING.EDGE[i], "=",fixed=T))[2]
  xx<-as.numeric(unlist(strsplit(xx, "%,",fixed=T))[1])
  res_1[,'enriched_percentage'][i]<-0.01*xx
}

res_1$number_genes_enriched<-round(res_1$SIZE*res_1$enriched_percentage,digits=0)



res_UP_DN<-rbind(res,res_1)






tag1<-grep('ANTIGEN',res_UP_DN$NAME)
tag2<-grep('LYMPHOCYTE',res_UP_DN$NAME)
tag3<-grep('B_CELL',res_UP_DN$NAME)
tag4<-grep('T_CELL',res_UP_DN$NAME)
tag5<-grep('TCR',res_UP_DN$NAME)
tag6<-grep('FC_GAMMA_R',res_UP_DN$NAME)
tag7<-grep('NATURAL_KILLER',res_UP_DN$NAME)
tag8<-grep('CD40',res_UP_DN$NAME)
tag9<-grep('B_LYMPHOCYTE',res_UP_DN$NAME)
tag10<-grep('TNF',res_UP_DN$NAME)
tag11<-grep('INFLAMMATORY',res_UP_DN$NAME)
tag12<-grep('IL12',res_UP_DN$NAME)
tag13<-grep('BCR',res_UP_DN$NAME)
tag14<-grep('SREBF',res_UP_DN$NAME)

tag15<-grep('TCA',res_UP_DN$NAME)
tag16<-grep('GLYCO',res_UP_DN$NAME)
tag17<-grep('FATTY_ACID',res_UP_DN$NAME)
tag18<-grep('ADIPO',res_UP_DN$NAME)
tag19<-grep('CHOLESTEROL',res_UP_DN$NAME)



tag<-c(tag1,tag2,tag3,tag4,tag5,tag6,tag7,tag8,tag9,tag10,tag11,tag12,tag13,tag14,tag15,tag16,tag17,tag18,tag19)
res_use<-res_UP_DN[tag,]
res_use<-res_use[res_use$FDR.q.val<0.1,]

res_use<-res_use[order(-res_use$NES),]
res_use$NAME<-factor(res_use$NAME,levels=unique(res_use$NAME))

myoutf<-'/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/Analysis_final_Exp1/Exclude_outliers/ileum/GSEA/Dotplot/Ileum_sham_vs_SBR_vehicle.pdf'
pdf(myoutf,width=20,height=15)
ggplot(res_use, aes(x = NES, y = NAME)) + 
               geom_point(aes(size = number_genes_enriched, color = FDR.q.val)) +
               theme_minimal(base_size = 16) +  # 使用 theme_minimal 去掉背景格子，并调大基础字号
               scale_size(range = c(4, 10)) +  # 增大 dot size 的范围
               scale_colour_gradient(limits=c(0, 0.10), low="red") +
               ylab(NULL) +
               ggtitle("GSEA") +
               geom_vline(xintercept = 0, linetype="dashed", color = "black") +  # 在 NES=0 位置添加一条竖直虚线
               theme(axis.text.y = element_text(size = 14),  # 调整y轴文字的大小
                     axis.text.x = element_text(size = 14),  # 调整x轴文字的大小
                     plot.title = element_text(size = 18, hjust = 0.5))

dev.off()












#(3) GSEA SBR vehicle versus WUSTL0717
myinf<-'/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/Analysis_final_Exp1/Exclude_outliers/ileum/GSEA/TSV/SBR_vehicle_VS_WUSTL0717/Ileum.Postana.SBR.Vehicle_VS_Ileum.Postana.SBR.WUSTL0717_UP.tsv'
res <- read.table(myinf,sep="\t",quote=NULL,header=TRUE)
res[,'enriched_percentage']<-rep(0,nrow(res))

for (i in 1:nrow(res))
{
  xx<-unlist(strsplit(res$LEADING.EDGE[i], "=",fixed=T))[2]
  xx<-as.numeric(unlist(strsplit(xx, "%,",fixed=T))[1])
  res[,'enriched_percentage'][i]<-0.01*xx
}

res$number_genes_enriched<-round(res$SIZE*res$enriched_percentage,digits=0)




myinf<-'/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/Analysis_final_Exp1/Exclude_outliers/ileum/GSEA/TSV/SBR_vehicle_VS_WUSTL0717/Ileum.Postana.SBR.Vehicle_VS_Ileum.Postana.SBR.WUSTL0717_DN.tsv'
res_1 <- read.table(myinf,sep="\t",quote=NULL,header=TRUE)
res_1[,'enriched_percentage']<-rep(0,nrow(res_1))

for (i in 1:nrow(res_1))
{
  xx<-unlist(strsplit(res_1$LEADING.EDGE[i], "=",fixed=T))[2]
  xx<-as.numeric(unlist(strsplit(xx, "%,",fixed=T))[1])
  res_1[,'enriched_percentage'][i]<-0.01*xx
}

res_1$number_genes_enriched<-round(res_1$SIZE*res_1$enriched_percentage,digits=0)



res_UP_DN<-rbind(res,res_1)






tag1<-grep('ANTIGEN',res_UP_DN$NAME)
tag2<-grep('INTERFERON',res_UP_DN$NAME)
tag3<-grep('B_CELL',res_UP_DN$NAME)
tag4<-grep('IFN_GAMMA',res_UP_DN$NAME)
tag5<-grep('TCR',res_UP_DN$NAME)
tag6<-grep('NEUTROPHIL',res_UP_DN$NAME)
tag7<-grep('TNF',res_UP_DN$NAME)
tag8<-grep('CHEMOKINE',res_UP_DN$NAME)
tag9<-grep('B_LYMPHOCYTE',res_UP_DN$NAME)
tag10<-grep('IL4',res_UP_DN$NAME)
tag11<-grep('INFLAMMATORY',res_UP_DN$NAME)


tag14<-grep('PPAR',res_UP_DN$NAME)

tag15<-grep('LIPOPROTEIN',res_UP_DN$NAME)
tag16<-grep('ABC',res_UP_DN$NAME)
tag17<-grep('FATTY_ACID',res_UP_DN$NAME)
tag18<-grep('ADIPO',res_UP_DN$NAME)
tag19<-grep('CHOLESTEROL',res_UP_DN$NAME)
tag20<-grep('PYRUVATE',res_UP_DN$NAME)
tag21<-grep('GLUCO',res_UP_DN$NAME)


tag<-c(tag1,tag2,tag3,tag4,tag5,tag6,tag7,tag8,tag9,tag10,tag11,tag12,tag13,tag14,tag15,tag16,tag17,tag18,tag19,tag20,tag21)
res_use<-res_UP_DN[tag,]
res_use<-res_use[res_use$FDR.q.val<0.1,]

res_use<-res_use[order(-res_use$NES),]
res_use$NAME<-factor(res_use$NAME,levels=unique(res_use$NAME))

myoutf<-'/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/Analysis_final_Exp1/Exclude_outliers/ileum/GSEA/Dotplot/Ileum_SBR_vehicle_vs_WUSTL0717.pdf'
pdf(myoutf,width=20,height=17)
ggplot(res_use, aes(x = NES, y = NAME)) + 
               geom_point(aes(size = number_genes_enriched, color = FDR.q.val)) +
               theme_minimal(base_size = 16) +  # 使用 theme_minimal 去掉背景格子，并调大基础字号
               scale_size(range = c(4, 10)) +  # 增大 dot size 的范围
               scale_colour_gradient(limits=c(0, 0.10), low="red") +
               ylab(NULL) +
               ggtitle("GSEA") +
               geom_vline(xintercept = 0, linetype="dashed", color = "black") +  # 在 NES=0 位置添加一条竖直虚线
               theme(axis.text.y = element_text(size = 14),  # 调整y轴文字的大小
                     axis.text.x = element_text(size = 14),  # 调整x轴文字的大小
                     plot.title = element_text(size = 18, hjust = 0.5))

dev.off()

































myinf<-'/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/Analysis_final_Exp1/Exclude_outliers/GSEA/liver_results/SBR_vehicle_vs_WUSTL/Vehicle_VS_WUSTL_UP.tsv'
res <- read.table(myinf,sep="\t",quote=NULL,header=TRUE)
res[,'enriched_percentage']<-rep(0,nrow(res))

for (i in 1:nrow(res))
{
  xx<-unlist(strsplit(res$LEADING.EDGE[i], "=",fixed=T))[2]
  xx<-as.numeric(unlist(strsplit(xx, "%,",fixed=T))[1])
  res[,'enriched_percentage'][i]<-0.01*xx
}

res$number_genes_enriched<-round(res$SIZE*res$enriched_percentage,digits=0)

tag1<-grep('FIBRO',res$NAME)
tag2<-grep('ECM',res$NAME)
tag3<-grep('FOCAL',res$NAME)
tag4<-grep('COLLAGEN',res$NAME)
tag5<-grep('PI3K',res$NAME)

tag<-c(tag1,tag2,tag3,tag4,tag5)
res_use<-res[tag,]
res_use<-res_use[res_use$FDR.q.val<0.1,]

res_use<-res_use[order(-res_use$NES),]
res_use$NAME<-factor(res_use$NAME,levels=unique(res_use$NAME))

myoutf<-'/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/Analysis_final_Exp1/Exclude_outliers/GSEA/liver_results/SBR_vehicle_vs_WUSTL/Vehicle_VS_WUSTL_ECM_fibrosis_DotPlot.pdf'
pdf(myoutf,width=20,height=10)
ggplot(res_use, aes(x = NES, y = NAME)) + 
               geom_point(aes(size = number_genes_enriched, color = FDR.q.val)) +
               theme_bw(base_size = 14) +
        scale_colour_gradient(limits=c(0, 0.10), low="red") +
        ylab(NULL) +
        ggtitle("GSEA")

dev.off()












myinf<-'/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/Analysis_final_Exp1/Exclude_outliers/GSEA/liver_results/SBR_vehicle_vs_WUSTL/Vehicle_VS_WUSTL_UP.tsv'
res <- read.table(myinf,sep="\t",quote=NULL,header=TRUE)
res[,'enriched_percentage']<-rep(0,nrow(res))

for (i in 1:nrow(res))
{
  xx<-unlist(strsplit(res$LEADING.EDGE[i], "=",fixed=T))[2]
  xx<-as.numeric(unlist(strsplit(xx, "%,",fixed=T))[1])
  res[,'enriched_percentage'][i]<-0.01*xx
}

res$number_genes_enriched<-round(res$SIZE*res$enriched_percentage,digits=0)

tag1<-grep('FIBRO',res$NAME)
tag2<-grep('ECM',res$NAME)
tag3<-grep('FOCAL',res$NAME)
tag4<-grep('COLLAGEN',res$NAME)
tag5<-grep('PI3K',res$NAME)

tag<-c(tag1,tag2,tag3,tag4,tag5)
res_use<-res[tag,]
res_use<-res_use[res_use$FDR.q.val<0.1,]

res_use<-res_use[order(-res_use$NES),]
res_use$NAME<-factor(res_use$NAME,levels=unique(res_use$NAME))

myoutf<-'/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/Analysis_final_Exp1/Exclude_outliers/GSEA/liver_results/SBR_vehicle_vs_WUSTL/Vehicle_VS_WUSTL_ECM_fibrosis_DotPlot.pdf'
pdf(myoutf,width=20,height=10)
ggplot(res_use, aes(x = NES, y = NAME)) + 
               geom_point(aes(size = number_genes_enriched, color = FDR.q.val)) +
               theme_bw(base_size = 14) +
        scale_colour_gradient(limits=c(0, 0.10), low="red") +
        ylab(NULL) +
        ggtitle("GSEA")

dev.off()
















myinf<-'/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/Analysis_final_Exp1/Exclude_outliers/GSEA/liver_results/SBR_vehicle_vs_WUSTL/Vehicle_VS_WUSTL_UP.tsv'
res <- read.table(myinf,sep="\t",quote=NULL,header=TRUE)
res[,'enriched_percentage']<-rep(0,nrow(res))

for (i in 1:nrow(res))
{
  xx<-unlist(strsplit(res$LEADING.EDGE[i], "=",fixed=T))[2]
  xx<-as.numeric(unlist(strsplit(xx, "%,",fixed=T))[1])
  res[,'enriched_percentage'][i]<-0.01*xx
}

res$number_genes_enriched<-round(res$SIZE*res$enriched_percentage,digits=0)

tag1<-grep('FIBRO',res$NAME)
tag2<-grep('ECM',res$NAME)
tag3<-grep('FOCAL',res$NAME)
tag4<-grep('COLLAGEN',res$NAME)
tag5<-grep('PI3K',res$NAME)

tag<-c(tag1,tag2,tag3,tag4,tag5)
res_use<-res[tag,]
res_use<-res_use[res_use$FDR.q.val<0.1,]

res_use<-res_use[order(-res_use$NES),]
res_use$NAME<-factor(res_use$NAME,levels=unique(res_use$NAME))

myoutf<-'/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/Analysis_final_Exp1/Exclude_outliers/GSEA/liver_results/SBR_vehicle_vs_WUSTL/Vehicle_VS_WUSTL_ECM_fibrosis_DotPlot.pdf'
pdf(myoutf,width=20,height=10)
ggplot(res_use, aes(x = NES, y = NAME)) + 
               geom_point(aes(size = number_genes_enriched, color = FDR.q.val)) +
               theme_bw(base_size = 14) +
        scale_colour_gradient(limits=c(0, 0.10), low="red") +
        ylab(NULL) +
        ggtitle("GSEA")

dev.off()





























#######3GSVA 10252024 Use this, for small intestine data, do not need to combat. For liver data, first exclude samples then combat
#######SI only

myinf1 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/SR003137/download/all.gene_TPM.tsv"
data1<-read.table(myinf1, sep="\t",quote = "", header=T)
tag1<-grep('sample',colnames(data1))
data1<-unique(data1[,c(3,tag1)])
tag1<-duplicated(data1$external_gene_name)
data1<-data1[!tag1,]
row.names(data1)<-data1$external_gene_name
tag1<-grep('sample',colnames(data1))
data1<-data1[,tag1]



myinf2 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/SR004337/all.gene_TPM.tsv"
data2<-read.table(myinf2, sep="\t",quote = "", header=T)
tag2<-grep('sample',colnames(data2))
data2<-unique(data2[,c(3,tag2)])
tag2<-duplicated(data2$external_gene_name)
data2<-data2[!tag2,]
row.names(data2)<-data2$external_gene_name
tag2<-grep('sample',colnames(data2))
data2<-data2[,tag2]


tag<-colnames(data1)%in%colnames(data2)
data1<-data1[,!tag]


xx1 <- row.names(data1)
xx2 <- row.names(data2)

com <- Reduce(intersect, list(xx1,xx2))

data1<-data1[com,]
data2<-data2[com,]

data1<-log2(data1+1)
data2<-log2(data2+1)

data_log2 <- cbind(data1,data2)
tag<-apply(data_log2,1,function(x){sum(x)==0})
data_log2<-data_log2[!tag,]




tag1<-c('sample.L9','sample.L13','sample.L24','sample.L29','sample.L1','sample.L34')
tag2<-c('sample.L14','sample.L21','sample.L6','sample.L18','sample.L2','sample.L31','sample.L10','sample.L26','sample.L36')
tag3<-c('sample.L3','sample.L22','sample.L37','sample.L11','sample.L27','sample.L7','sample.L15','sample.L19','sample.L32')
tag4<-c('sample.A9','sample.A13','sample.A29')
tag5<-c('sample.A2','sample.A21','sample.A31','sample.A36')
tag6<-c('sample.A3','sample.A22','sample.A7','sample.A19')
tag7<-c('sample.D9','sample.D24','sample.D29')
tag8<-c('sample.D31','sample.D6','sample.D18','sample.D36')
tag9<-c('sample.D3','sample.D37','sample.D7','sample.D32')
colnames(data_log2)[colnames(data_log2)%in%tag1]<-paste0('Liver.Sham_',seq(1,length(tag1),1))
colnames(data_log2)[colnames(data_log2)%in%tag2]<-paste0('Liver.SBR.Vehicle_',seq(1,length(tag2),1))
colnames(data_log2)[colnames(data_log2)%in%tag3]<-paste0('Liver.SBR.WUSTL0717_',seq(1,length(tag3),1))
colnames(data_log2)[colnames(data_log2)%in%tag4]<-paste0('Ileum.Postana.Sham_',seq(1,length(tag4),1))
colnames(data_log2)[colnames(data_log2)%in%tag5]<-paste0('Ileum.Postana.SBR.Vehicle_',seq(1,length(tag5),1))
colnames(data_log2)[colnames(data_log2)%in%tag6]<-paste0('Ileum.Postana.SBR.WUSTL0717_',seq(1,length(tag6),1))
colnames(data_log2)[colnames(data_log2)%in%tag7]<-paste0('Duodenum.Sham_',seq(1,length(tag7),1))
colnames(data_log2)[colnames(data_log2)%in%tag8]<-paste0('Duodenum.SBR.Vehicle_',seq(1,length(tag8),1))
colnames(data_log2)[colnames(data_log2)%in%tag9]<-paste0('Duodenum.SBR.WUSTL0717_',seq(1,length(tag9),1))

tag<-apply(data_log2,1,function(x){sum(x)==0})
data_log2<-data_log2[!tag,]


tag<-c(grep('Duodenum',colnames(data_log2)),grep('Ileum',colnames(data_log2)))
data_log2_intestine<-data_log2[,tag]


gene_means <- rowMeans(data_log2_intestine)
gene_variances <- apply(data_log2_intestine, 1, var)

# Step 2: Set thresholds to filter low-quality genes
# For example, keep genes with mean expression > 1 and variance > 0.5
mean_threshold <- 1
variance_threshold <- 0.5
filter_criteria <- (gene_means > mean_threshold) & (gene_variances > variance_threshold)

# Step 3: Apply the filter
data_log2_intestine <- data_log2_intestine[filter_criteria, ]



library(GSEABase)
myinf<-'/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/revision/Myeloid_PW_Adult_only_Human_nCountRNA_modified/PCx/res.0.6/GSVA/TF_activity/MSigDB/c3.tft.v2023.1.Hs.symbols.gmt'
database<-getGmt(myinf)

TF_use<-database[c(grep("^LXR", names(database)))]

LXR_Q3 <- geneIds(TF_use[[1]])
LXR_Q3_DR4 <- geneIds(TF_use[[2]])
LXR_overlap<-LXR_Q3[LXR_Q3%in%LXR_Q3_DR4]
LXRalpha <- c("ABCA1", "ABCD2", "ABCG1", "AEBP1", "BARD1", "BRCA1", "CETP", "CRP", "DHCR24", "ENG", "FABP4", "FOXO1", "HTT", 
"IFNG", "IGF1", "IL10", "JAK1", "LIPG", "LPCAT3", "LPL", "MIR1-1", "MIR206", "MIR613", "NFKB1", "NR0B2", "NR1H2", "NR1H3", "PLIN2",
 "PPARA", "PPARG", "PPARGC1B", "PRKACA", "PRMT3", "PTX3", "RORA", "RORC", "RXRA", "SCARB1", "SCD", "SREBF1", "STAT1", "THRB", "UGT1A1")


LXR_taget<-list(LXR_Q3,LXR_Q3_DR4,LXR_overlap,LXRalpha)
names(LXR_taget)<-c('LXR_Q3','LXR_Q3_DR4','LXR_overlap','LXRalpha')
#screen -r 238076

row.names(data_log2_intestine)<-toupper(row.names(data_log2_intestine))
gsva.es <- gsva(as.matrix(data_log2_intestine), LXR_taget, verbose=FALSE)

#gsva.es<-gsva.es[,c(1:8,10:34)]
gsva.es <- apply(gsva.es,1, function(arg) (arg-mean(arg))/sd(arg))
myoutf = "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/Analysis_final_Exp1/Exclude_outliers/GSVA_LXR_scale_only_intestine.xls"
 write.table(gsva.es,myoutf,sep="\t",quote=F)





#######only exclude sham liver 6


myinf1 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/SR003137/download/all.gene_TPM.tsv"
data1<-read.table(myinf1, sep="\t",quote = "", header=T)
tag1<-grep('sample',colnames(data1))
data1<-unique(data1[,c(3,tag1)])
tag1<-duplicated(data1$external_gene_name)
data1<-data1[!tag1,]
row.names(data1)<-data1$external_gene_name
tag1<-grep('sample',colnames(data1))
data1<-data1[,tag1]



myinf2 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/SR004337/all.gene_TPM.tsv"
data2<-read.table(myinf2, sep="\t",quote = "", header=T)
tag2<-grep('sample',colnames(data2))
data2<-unique(data2[,c(3,tag2)])
tag2<-duplicated(data2$external_gene_name)
data2<-data2[!tag2,]
row.names(data2)<-data2$external_gene_name
tag2<-grep('sample',colnames(data2))
data2<-data2[,tag2]


tag<-colnames(data1)%in%colnames(data2)
data1<-data1[,!tag]


xx1 <- row.names(data1)
xx2 <- row.names(data2)

com <- Reduce(intersect, list(xx1,xx2))

data1<-data1[com,]
data2<-data2[com,]

data1<-log2(data1+1)
data2<-log2(data2+1)

data_log2 <- cbind(data1,data2)
tag<-apply(data_log2,1,function(x){sum(x)==0})
data_log2<-data_log2[!tag,]




tag1<-c('sample.L9','sample.L13','sample.L24','sample.L29','sample.L1','sample.L34')
tag2<-c('sample.L14','sample.L21','sample.L6','sample.L18','sample.L2','sample.L31','sample.L10','sample.L26','sample.L36')
tag3<-c('sample.L3','sample.L22','sample.L37','sample.L11','sample.L27','sample.L7','sample.L15','sample.L19','sample.L32')
tag4<-c('sample.A9','sample.A13','sample.A29')
tag5<-c('sample.A2','sample.A21','sample.A31','sample.A36')
tag6<-c('sample.A3','sample.A22','sample.A7','sample.A19')
tag7<-c('sample.D9','sample.D24','sample.D29')
tag8<-c('sample.D31','sample.D6','sample.D18','sample.D36')
tag9<-c('sample.D3','sample.D37','sample.D7','sample.D32')
colnames(data_log2)[colnames(data_log2)%in%tag1]<-paste0('Liver.Sham_',seq(1,length(tag1),1))
colnames(data_log2)[colnames(data_log2)%in%tag2]<-paste0('Liver.SBR.Vehicle_',seq(1,length(tag2),1))
colnames(data_log2)[colnames(data_log2)%in%tag3]<-paste0('Liver.SBR.WUSTL0717_',seq(1,length(tag3),1))
colnames(data_log2)[colnames(data_log2)%in%tag4]<-paste0('Ileum.Postana.Sham_',seq(1,length(tag4),1))
colnames(data_log2)[colnames(data_log2)%in%tag5]<-paste0('Ileum.Postana.SBR.Vehicle_',seq(1,length(tag5),1))
colnames(data_log2)[colnames(data_log2)%in%tag6]<-paste0('Ileum.Postana.SBR.WUSTL0717_',seq(1,length(tag6),1))
colnames(data_log2)[colnames(data_log2)%in%tag7]<-paste0('Duodenum.Sham_',seq(1,length(tag7),1))
colnames(data_log2)[colnames(data_log2)%in%tag8]<-paste0('Duodenum.SBR.Vehicle_',seq(1,length(tag8),1))
colnames(data_log2)[colnames(data_log2)%in%tag9]<-paste0('Duodenum.SBR.WUSTL0717_',seq(1,length(tag9),1))

tag<-apply(data_log2,1,function(x){sum(x)==0})
data_log2<-data_log2[!tag,]


grep('Liver.Sham_6',colnames(data_log2))
data_log2<-data_log2[,1:45]


gene_means <- rowMeans(data_log2)
gene_variances <- apply(data_log2, 1, var)

# Step 2: Set thresholds to filter low-quality genes
# For example, keep genes with mean expression > 1 and variance > 0.5
mean_threshold <- 1
variance_threshold <- 0.5
filter_criteria <- (gene_means > mean_threshold) & (gene_variances > variance_threshold)

# Step 3: Apply the filter
data_log2 <- data_log2[filter_criteria, ]

sample_names <- colnames(data_log2)
batch <- ifelse(grepl("Liver", sample_names), "Liver", "Other")
batch <- as.factor(batch)
batch_df <- data.frame(batch = batch)
mod <- model.matrix(~1, data = batch_df)
data_log2 <- ComBat(dat = as.matrix(data_log2), batch = batch, mod = mod)




library(GSEABase)
myinf<-'/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Analysis/revision/Myeloid_PW_Adult_only_Human_nCountRNA_modified/PCx/res.0.6/GSVA/TF_activity/MSigDB/c3.tft.v2023.1.Hs.symbols.gmt'
database<-getGmt(myinf)

TF_use<-database[c(grep("^LXR", names(database)))]

LXR_Q3 <- geneIds(TF_use[[1]])
LXR_Q3_DR4 <- geneIds(TF_use[[2]])
LXR_overlap<-LXR_Q3[LXR_Q3%in%LXR_Q3_DR4]
LXRalpha <- c("ABCA1", "ABCD2", "ABCG1", "AEBP1", "BARD1", "BRCA1", "CETP", "CRP", "DHCR24", "ENG", "FABP4", "FOXO1", "HTT", 
"IFNG", "IGF1", "IL10", "JAK1", "LIPG", "LPCAT3", "LPL", "MIR1-1", "MIR206", "MIR613", "NFKB1", "NR0B2", "NR1H2", "NR1H3", "PLIN2",
 "PPARA", "PPARG", "PPARGC1B", "PRKACA", "PRMT3", "PTX3", "RORA", "RORC", "RXRA", "SCARB1", "SCD", "SREBF1", "STAT1", "THRB", "UGT1A1")


LXR_taget<-list(LXR_Q3,LXR_Q3_DR4,LXR_overlap,LXRalpha)
names(LXR_taget)<-c('LXR_Q3','LXR_Q3_DR4','LXR_overlap','LXRalpha')
#screen -r 238076

row.names(data_log2)<-toupper(row.names(data_log2))
gsva.es <- gsva(as.matrix(data_log2), LXR_taget, verbose=FALSE)

#gsva.es<-gsva.es[,c(1:8,10:34)]
gsva.es <- apply(gsva.es,1, function(arg) (arg-mean(arg))/sd(arg))

myoutf = "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/Analysis_final_Exp1/With_all_liver_sample/GSVA_LXR_scale_without_sham_6_scale_use_this.xls"

 write.table(gsva.es,myoutf,sep="\t",quote=F)




























######Plot genes Gwen wants:

myinf1 <- "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/Analysis_final_Exp1/Exclude_outliers/Log2_metadata_after_combat_no_quantile.txt"
data_log2 <- read.table(myinf1,sep="\t",quote=NULL)


data_use<-data_log2
data_use_liver<-data_use[,grep('Liver',colnames(data_use))]
data_use_duo<-data_use[,grep('Duodenum',colnames(data_use))]
data_use_ileum<-data_use[,grep('Ileum',colnames(data_use))]

data_use_intestine<-cbind(data_use_duo,data_use_ileum)

data_log2<-apply(data_log2,1, function(arg) (arg-mean(arg))/sd(arg))
data_use_liver <- apply(data_use_liver,1, function(arg) (arg-mean(arg))/sd(arg))
data_use_duo <- apply(data_use_duo,1, function(arg) (arg-mean(arg))/sd(arg))
data_use_ileum <- apply(data_use_ileum,1, function(arg) (arg-mean(arg))/sd(arg))
data_use_intestine<-apply(data_use_intestine,1, function(arg) (arg-mean(arg))/sd(arg))

data_log2<-t(data_log2)
res_liver <- t(data_use_liver)
res_duo <- t(data_use_duo)
res_ileum <- t(data_use_ileum)
data_use_intestine<-t(data_use_intestine)
#res<-res[,c("2",'5','4','8','9','0','1','7','3','6')]


genes_plot<- c('Slc10a2','Reg3g','Reg3b','Fabp1','Mttp','Apoa1','Abca1')
data_plot<-data_log2[genes_plot,]
liver_plot<-res_liver[genes_plot,]
due_plot<-res_duo[genes_plot,]
ileum_plot<-res_ileum[genes_plot,]




myoutf = "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/Analysis_final_Exp1/Exclude_outliers/genes_plot/scaled_selected_gene_exp_liver.xls"
write.table(t(liver_plot),myoutf,sep="\t",quote=F)



myoutf = "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/Analysis_final_Exp1/Exclude_outliers/genes_plot/scaled_selected_gene_exp_duo.xls"
write.table(t(due_plot),myoutf,sep="\t",quote=F)


myoutf = "/storage1/fs1/gjrandolph/Active/RNAseq_backup/jhan/Ayoung/Analysis_final_Exp1/Exclude_outliers/genes_plot/scaled_selected_gene_exp_ileum-1.xls"
write.table(t(ileum_plot),myoutf,sep="\t",quote=F)

