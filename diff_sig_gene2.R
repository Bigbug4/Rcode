# using limma
library(limma)
library(edgeR)
DIR=getwd()
setwd("E:\\Rcode\\data")

path_expr_data <- "sum.csv"
sum_expr_data <- read.csv(path_expr_data)
rownames(sum_expr_data) <- sum_expr_data$ENSGgeneID
names(sum_expr_data) <- gsub(".","-",names(sum_expr_data),fixed = TRUE)

sample_choosed <- read.csv("sample_list.csv")
names(sample_choosed) <- "Sample.ID"

expr_data_choosed <- sum_expr_data[-c(1:5),which(names(sum_expr_data)%in%sample_choosed$Sample.ID)]
counts <- as.matrix(expr_data_choosed)
# counts_t <- t(counts)
# write.table(counts_t,"counts_t.txt",sep = "\t")
write.csv(counts,"counts.csv")

dge <- DGEList(counts = counts)
dge <- calcNormFactors(dge)
CPM <- cpm(dge, log=FALSE, prior.count=2)
# 归一化
logCPM <- cpm(dge, log=TRUE, prior.count=3)
group_list <- gsub("..*11A$",replacement = "Normal",names(expr_data_choosed))
group_list <- gsub("..*01A$",replacement = "Tumor",group_list)
group_list <- factor(group_list)
design <- model.matrix(~0+group_list)
write.csv(design,"design.csv")

cont.matrix <- makeContrasts(group_listTumor-group_listNormal,levels = design)
row.names(cont.matrix) <- levels(group_list)
colnames(design) <- levels(group_list)
rownames(design) <- colnames(counts)
fit <- lmFit(logCPM, design)

##
# fit <- eBayes(fit, trend=TRUE)
# output <- topTable(fit, coef=2,n=Inf)
# sum(output$adj.P.Val<0.01)
# sum(output$logFC>=2&output$adj.P.Val<0.01)
# sum(output$logFC<=-1&output$adj.P.Val<0.01)
##

fit2 <- contrasts.fit(fit,contrasts = cont.matrix)
fit2 <- eBayes(fit2)  
result <- topTable(fit2,adjust="BH",n=Inf)

sum(result$logFC>=1&result$adj.P.Val<0.01&result$AveExpr>=5&result$gene_biotype=="protein_coding")
sum(result$logFC<=-1&result$adj.P.Val<0.01&result$AveExpr>=5)
#result.up <- result[which(result$logFC>=1&result$adj.P.Val<0.01&result$AveExpr>=5),]
#result.down <- result[which(result$logFC<=-1&result$adj.P.Val<0.01&result$AveExpr>=5),]

# ע??
expr_data_annotation <- read.csv("ensembel_gene_all.csv")
expr_data_choosed.annotation <- expr_data_annotation[expr_data_annotation$gene_id_v==row.names(expr_data_choosed),]

expr_data_choosed.protein <- expr_data_choosed[which(expr_data_annotation$gene_biotype=="protein_coding"),]
#expr_data_choosed.mirna <- expr_data_choosed[which(expr_data_annotation$gene_biotype=="miRNA"),]
gene_id2name <- function(result=result.up,annotation = expr_data_annotation){
  gene_id_v <- row.names(result)
  result$gene_name <- gene_id_v
  row.names(annotation) <- annotation$gene_id_v
  result.annotation <- annotation[gene_id_v,]
  return(result.annotation)
}

result.annotation <- gene_id2name(result = result)
result <- cbind(result,result.annotation)
result.up_protein <- result[which(result$logFC>=1&result$adj.P.Val<0.01&result$gene_biotype=="protein_coding"&result$AveExpr>=5),]
result.down_protein <- result[which(result$logFC<=-1&result$adj.P.Val<0.01&result$gene_biotype=="protein_coding"&result$AveExpr>=5),]

result.protein <- rbind(result.up_protein,result.down_protein)
diff_protein <- result.protein[,c(1,9,10)]

# result.mirna <- result[which((result$logFC<=-0.5|result$logFC>=0.5)&result$adj.P.Val<0.05&result$gene_biotype=="miRNA"),]
# result.lincrna <- result[which((result$logFC<=-1|result$logFC>=1)&result$adj.P.Val<0.05&result$gene_biotype=="lincRNA"),]
# diff_lincrna <- result.lincrna[,c(1,9,10)]

result.lncrna <- result[which((result$logFC<=-1|result$logFC>=1)&result$adj.P.Val<0.05&result$gene_biotype%in%c("lincRNA","non coding","3prime_overlapping_ncRNA","antisense","retained_intron","sense_intronic","sense_overlapping","macro_lncRNA")),]
diff_lncrna <- result.lncrna[,c(1,9,10)]

# result.snrna <- result[which((result$logFC<=-0.5|result$logFC>=0.5)&result$adj.P.Val<0.05&result$gene_biotype=="snRNA"),]
# result.snorna <- result[which((result$logFC<=-0.5|result$logFC>=0.5)&result$adj.P.Val<0.05&result$gene_biotype=="snoRNA"),]
# 
# house_keeping_gene <- result[which(result$logFC<0.1&result$logFC>-0.1&result$P.Value<0.05),]
# house_keeping_mirna <- result.mirna[which(result.mirna$logFC<0.01&result.mirna>-0.01&result.mirna$adj.P.Val<0.05),]

# write.csv(result.up,"diff_expr_up.csv")
# write.csv(result.down,"diff_expr_down.csv")
# write.csv(result.down_protein,"diff_expr_down_protein.csv")
# write.csv(result.up_protein,"diff_expr_up_protein.csv")
# write.csv(result.lincrna,"diff_expr_lincrna.csv")
write.csv(result,"result_annotation_v2.csv")
write.csv(result.protein,"diff_expr_protein.csv")
write.csv(result.lncrna,"diff_expr_lncrna.csv")


# different expression of miRNA
expr_mirna <- read.csv("miRNA_read_count_merged.csv",header = TRUE)
# sample_choosed <- read.csv("sample_choosed.csv")

all_sample_mirna <- names(expr_mirna)
names(expr_mirna) <- gsub(".",replacement = "-",all_sample_mirna,fixed = TRUE)
expr_mirna_choosed <- expr_mirna[,which(names(expr_mirna)%in%sample_choosed$Sample.ID)]
row.names(expr_mirna_choosed) <- expr_mirna$miRNA_ID
expr_mirna <- expr_mirna_choosed
rm(expr_mirna_choosed)

mirna_counts <- as.matrix(expr_mirna)
dge_mirna <- DGEList(counts = mirna_counts)
dge_mirna <- calcNormFactors(dge_mirna)
CPM_mirna <- cpm(dge_mirna, log=FALSE, prior.count=2)
logCPM_mirna <- cpm(dge_mirna, log=TRUE, prior.count=3)

group_list <- gsub("..*11A$",replacement = "Normal",names(expr_mirna))
group_list <- gsub("..*01A$",replacement = "Tumor",group_list)
group_list <- factor(group_list)
design <- model.matrix(~0+group_list)

cont.matrix_mirna <- makeContrasts(group_listTumor-group_listNormal,levels = design)
row.names(cont.matrix_mirna) <- levels(group_list)
colnames(design) <- levels(group_list)
rownames(design) <- colnames(mirna_counts)
fit_mirna <- lmFit(logCPM_mirna, design)

##
# fit <- eBayes(fit, trend=TRUE)
# output <- topTable(fit, coef=2,n=Inf)
# sum(output$adj.P.Val<0.01)
# sum(output$logFC>=2&output$adj.P.Val<0.01)
# sum(output$logFC<=-1&output$adj.P.Val<0.01)
##


fit2_mirna <- contrasts.fit(fit = fit_mirna,contrasts = cont.matrix_mirna)
fit2_mirna <- eBayes(fit2_mirna)  
result.mirna <- topTable(fit2_mirna,adjust="BH",n=Inf)

sum(result.mirna$logFC>=0&result.mirna$adj.P.Val<0.01&result.mirna$P.Value<result.mirna$adj.P.Val)
sum(result.mirna$logFC<0&result.mirna$adj.P.Val<0.01&result.mirna$P.Value<result.mirna$adj.P.Val)

diff_mirna <- result.mirna[which((result.mirna$logFC>=1|result.mirna$logFC<=-1)&result.mirna$adj.P.Val<0.01&result.mirna$P.Value<result.mirna$adj.P.Val),]
expr_diff_mirna <- expr_mirna[which(row.names(expr_mirna)%in%row.names(diff_mirna)),]

write.csv(result.mirna,file = "result.mirna.csv")
write.csv(diff_mirna,file = "diff_expr_mirna.csv")
# write.table(row.names(diff_mirna),file = "diff_gene_mirna.txt",row.names = FALSE,col.names = FALSE,quote = FALSE)
write.csv(expr_diff_mirna,file = "expr_diff_mirna.csv")

# save.image("limma_result.Rdata")
# load("limma_result.Rdata")