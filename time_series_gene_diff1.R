DIR=getwd()
setwd("E:\\Rcode\\data")

clinical_choosed <- read.csv("clinical_merged.csv",row.names = 1)
summary(clinical_choosed$tumor_stage)

sample_stage_ia_c <- clinical_choosed$submitter_id[which(clinical_choosed$tumor_stage=="stage ia")]
sample_stage_ib_c <- clinical_choosed$submitter_id[which(clinical_choosed$tumor_stage=="stage ib")]
sample_stage_iia_c <- clinical_choosed$submitter_id[which(clinical_choosed$tumor_stage=="stage iia")]
sample_stage_iib_c <- clinical_choosed$submitter_id[which(clinical_choosed$tumor_stage=="stage iib")]
sample_stage_iiia_c <- clinical_choosed$submitter_id[which(clinical_choosed$tumor_stage=="stage iiia")]
sample_stage_iiib_c <- clinical_choosed$submitter_id[which(clinical_choosed$tumor_stage=="stage iiib")]
sample_stage_iiic_c <- clinical_choosed$submitter_id[which(clinical_choosed$tumor_stage=="stage iiic")]
sample_stage_iv_c <- clinical_choosed$submitter_id[which(clinical_choosed$tumor_stage=="stage iv")]

path_expr_data <- "sum.csv"
data <- read.csv(path_expr_data)
rownames(data) <- data$ENSGgeneID
names(data) <- gsub(".","-",names(data),fixed = TRUE)

all_sample <- read.csv("SampleID_and_Type.csv")

sample_stage_ia_s <- all_sample[which(all_sample$Case.ID%in%sample_stage_ia_c),]
sample_stage_ib_s <- all_sample[which(all_sample$Case.ID%in%sample_stage_ib_c),]
sample_stage_iia_s <- all_sample[which(all_sample$Case.ID%in%sample_stage_iia_c),]
sample_stage_iib_s <- all_sample[which(all_sample$Case.ID%in%sample_stage_iib_c),]
sample_stage_iiia_s <- all_sample[which(all_sample$Case.ID%in%sample_stage_iiia_c),]
sample_stage_iiib_s <- all_sample[which(all_sample$Case.ID%in%sample_stage_iiib_c),]
sample_stage_iiic_s <- all_sample[which(all_sample$Case.ID%in%sample_stage_iiic_c),]
sample_stage_iv_s <- all_sample[which(all_sample$Case.ID%in%sample_stage_iv_c),]

expr_stage_ia <- data[,which(names(data)%in%sample_stage_ia_s$Sample.ID)]
expr_stage_ib <- data[,which(names(data)%in%sample_stage_ib_s$Sample.ID)]
expr_stage_iia <- data[,which(names(data)%in%sample_stage_iia_s$Sample.ID)]
expr_stage_iib <- data[,which(names(data)%in%sample_stage_iib_s$Sample.ID)]
expr_stage_iiia <- data[,which(names(data)%in%sample_stage_iiia_s$Sample.ID)]
expr_stage_iiib <- data[,which(names(data)%in%sample_stage_iiib_s$Sample.ID)]
expr_stage_iiic <- data[,which(names(data)%in%sample_stage_iiic_s$Sample.ID)]
expr_stage_iv <- data[,which(names(data)%in%sample_stage_iv_s$Sample.ID)]

library(limma)
library(edgeR)

expr_data_annotation <- read.csv("ensembel_gene_all.csv")
expr_data_choosed.protein <- expr_stage_ia[which(expr_data_annotation$gene_biotype=="protein_coding"),]
sig_gene <- read.csv("sig_gene_v3.csv")

gene_id2name <- function(result=result.up,annotation = expr_data_annotation){
  gene_id_v <- row.names(result)
  result$gene_name <- gene_id_v
  row.names(annotation) <- annotation$gene_id_v
  result.annotation <- annotation[gene_id_v,]
  return(result.annotation)
}

deg_analysis<-function(expr_stage_data=expr_stage_ia){
  counts <- as.matrix(expr_stage_data[-c(1:5),])
  dge <- DGEList(counts = counts)
  dge <- calcNormFactors(dge)
  
  # 归一化
  counts<- cpm(dge, log=TRUE, prior.count=3)
  
  #hist(counts)
  #boxplot(counts, col=rainbow(4),main="expression value",las=2)
  
  group_list <- gsub("..*11A$",replacement = "Normal",colnames(counts))
  group_list <- gsub("..*01A$",replacement = "Tumor",group_list)
  group_list <- gsub("..*11B$",replacement = "Normal",group_list)
  group_list <- gsub("..*01B$",replacement = "Tumor",group_list)
  group_list <- factor(group_list)
  design <- model.matrix(~0+group_list)
  
  cont.matrix <- makeContrasts(group_listTumor-group_listNormal,levels = design)
  row.names(cont.matrix) <- levels(group_list)
  colnames(design) <- levels(group_list)
  rownames(design) <- colnames(counts)
  fit <- lmFit(counts, design)
  fit2 <- contrasts.fit(fit,contrasts = cont.matrix)
  fit2 <- eBayes(fit2)  
  result <- topTable(fit2,adjust="BH",n=Inf)
  
  sum(abs(result$logFC)>=1&result$P.Value<0.01)
  result.sig <- result[abs(result$logFC)>=1&result$P.Value<0.01,]
  
  result.annotation <- gene_id2name(result = result)
  result <- cbind(result,result.annotation)
  result.protein <- result[which(abs(result$logFC)>=1&result$P.Value<0.05&result$gene_biotype=="protein_coding"),]
  
  diff_protein <- result.protein[,c(1,9,10)]
  
  res<- diff_protein[which(diff_protein$gene_name%in%sig_gene$geneName),]
  return(res)
}

res_ia<-deg_analysis()
#res_ib<-deg_analysis(expr_stage_ib)
res_iia<-deg_analysis(expr_stage_iia)
res_iib<-deg_analysis(expr_stage_iib)
res_iiia<-deg_analysis(expr_stage_iiia)
res_iiib<-deg_analysis(expr_stage_iiib)
#res_iiic<-deg_analysis(expr_stage_iiic)
res_iv<-deg_analysis(expr_stage_iv)

write.csv(res_ia,"res_ia.csv")
write.csv(res_iia,"res_iia.csv")
write.csv(res_iib,"res_iib.csv")
write.csv(res_iiia,"res_iiia.csv")
write.csv(res_iiib,"res_iiib.csv")
write.csv(res_iv,"res_iv.csv")
save.image("time_sig_gene_diff1.RData")
#load("time_sig_gene_diff.RData")
