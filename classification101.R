# 分类

DIR=getwd()
setwd("E:\\Rcode\\data")
sig_gene <- read.csv("sig_gene_101.csv")

expr_all <- read.csv("sum.csv",row.names = 1)
expr_all <- expr_all[as.character(sig_gene$gene_id_v),]
row.names(expr_all) <- sig_gene$gene_name
names(expr_all) <- gsub("\\.","-",names(expr_all))
write.csv(expr_all,"expr_all_stomach101.csv")

# scale
library(edgeR)
counts <- as.matrix(expr_all)
dge <- DGEList(counts = counts)
dge <- calcNormFactors(dge)
CPM <- cpm(dge, log=FALSE, prior.count=2)
logCPM <- cpm(dge, log=TRUE, prior.count=3)
write.csv(CPM,"CPM101.csv")
write.csv(logCPM,"logCPM101.csv")

sample_list <- read.csv("sample_list.csv")
sample_list <- as.vector(t(sample_list))
CPM <- as.data.frame(CPM)
logCPM <- as.data.frame(logCPM)
expr_data_choosed_cpm <- CPM[,which(names(CPM)%in%sample_list)]
write.csv(expr_data_choosed_cpm,file = "expr_gene_data101_cpm.csv")

expr_data_choosed <- read.csv("expr_data_choosed.csv",row.names = 1)
sig_gene <- sig_gene[order(sig_gene$logFC,decreasing = TRUE),]
expr_data_choosed <- expr_data_choosed[as.character(sig_gene$gene_id_v),]
row.names(expr_data_choosed) <- sig_gene$GeneName
names(expr_data_choosed) <- gsub("\\.","-",names(expr_data_choosed))
write.csv(expr_data_choosed,file = "expr_gene_data101.csv")


# case_choosed <- read.table("case_choosed.txt",col.names = "Sample")
sample_choosed <- read.csv("sample_choosed.csv")
case_choosed <- summary(sample_choosed$Case.ID)
case_choosed <- names(case_choosed)[which(case_choosed==2)]
s <- sample_choosed[which(sample_choosed$Case.ID%in%case_choosed),]
s <- s[order(s$Case.ID),]
sample_choosed_20 <- s[1:20,]
#sample_choosed_30 <- s[1:30,]

sig_gene_up <- sig_gene[which(sig_gene$logFC>0),]
sig_gene_down <- sig_gene[which(sig_gene$logFC<0),]
# heatmap(log10(as.matrix(expr_data_choosed_cpm[,sample_list])))

# expr_data_choosed_log <- log10(expr_data_choosed)
# expr_data_choosed_log <- expr_data_choosed_log/colSums(expr_data_choosed_log)
# expr_up <- expr_data_choosed_log[as.character(sig_gene_up$geneName),]
# expr_down <- expr_data_choosed_log[as.character(sig_gene_down$geneName),]

comp_in_ssample <- function(expr_up,expr_down){
  n <- dim(expr_down)[2]
  comp_vector <- rep(0,n)
  names(comp_vector) <- names(expr_down)
  for (i in 1:n) {
    comp_vector[i] <- mean(expr_up[,i])/mean(expr_down[,i])
    
  }
  return(comp_vector)
}


combine_sample <- function(expr_data_choosed,sample_choosed,n){
  combine_matrix <- expr_data_choosed[,1:n]
  names(combine_matrix) <- unique(sample_choosed$Case.ID)
  for (i in 1:n) {
    a <- as.numeric(expr_data_choosed[,as.character(sample_choosed$Sample.ID[i])])
    b <- as.numeric(expr_data_choosed[,as.character(sample_choosed$Sample.ID[2*i])])
    combine_matrix[,i] <- (a+b*sum(a)/sum(b))/2
  }
  return(combine_matrix)
}

# s1 <- comp_in_ssample(expr_up,expr_down)
s_combine1 <- combine_sample(expr_data_choosed,sample_choosed_20,10)
#s_combine2 <- combine_sample(expr_data_choosed,sample_choosed_30,15)
s1 <- comp_in_ssample(log10(s_combine1[1:41,]),log10(s_combine1[42:55,]))
c_factor1 <- exp(sum(log(s1))/length(s1))
# s2 <- comp_in_ssample(log10(s_combine2[1:41,]),log10(s_combine2[42:55,]))
# c_factor2 <- exp(sum(log(s2))/length(s2))

# 正负样本
expr_all_log <- log10(expr_all)
expr_all_log <- as.data.frame(t(expr_all_log))
expr_all_log$class <- rownames(expr_all_log)
expr_all_log$class[which(grepl("-01",expr_all_log$class))] <- "yes"
expr_all_log$class[which(grepl("-11",expr_all_log$class))] <- "no"
write.csv(expr_all_log,"expr_all_log101.csv",row.names = FALSE)

# expr_negative_log <- expr_all_log[which(grepl("-11",names(expr_all_log)))]
# expr_negative_log <-expr_all_log <-expr_negative_log))
# expr_negative_log$class <- "no"
# write.csv(expr_negative_log,"expr_negative_log.csv",row.names = FALSE)
# expr_positive_log <- expr_all_log[which(grepl("-01",names(expr_all_log)))]
# expr_positive_log <- as.data.frame(t(expr_positive_log))
# expr_positive_log $class <- "yes"
# write.csv(expr_positive_log,"expr_positive_log.csv",row.names = FALSE)

sample_all <- names(expr_all)
# sample_all <- sample_all[which(grepl("-IN-",sample_all)|grepl("-BR-",sample_all)|grepl("-HU-",sample_all))]
write.csv(sample_all,"Sample_all_stomach.csv",row.names = FALSE)
sample_all_test <- sample_all[which(!sample_all%in%sample_choosed_20$Sample.ID)]
write.csv(sample_all_test,"Sample_all_test.csv",row.names = FALSE)

s3 <- comp_in_ssample(log10(expr_all[1:41,]),log10(expr_all[42:55,]))
s3 <- s3[order(s3)]
s3 <- s3[which(names(s3)%in%sample_all_test)]
write.csv(s3,"test_result101.csv")
Positive <- s3[which(grepl("-01",names(s3)))]
Negative <- s3[which(grepl("-11",names(s3)))]

TP1 <- Positive[Positive>c_factor1]
TN1 <- Negative[Negative<c_factor1]
ACC1 <- (length(TP1)+length(TN1))/length(s3)
sensitive1 <- length(TP1)/length(Positive)
precision1 <- length(TN1)/length(Negative)

# TP2 <- Positive[Positive>c_factor2]
# TN2 <- Negative[Negative<c_factor2]
# ACC2 <- (length(TP2)+length(TN2))/length(s3)
# sensitive2 <- length(TP2)/length(Positive)
# precision2 <- length(TN2)/length(Negative)

# ROC
roc_plot1 <- names(s3)
roc_plot1[which(grepl("-01",roc_plot1))] <- 1
roc_plot1[which(grepl("-11",roc_plot1))] <- 0
roc_plot1 <- as.numeric(roc_plot1)

roc_plot2 <- s3
names(roc_plot2) <- NULL
roc_plot2[roc_plot2>c_factor] <- 1
roc_plot2[roc_plot2!=1] <- 0
roc_plot <- data.frame(roc_plot1, roc_plot2)
write.csv(roc_plot,file ="roc_test101.csv",row.names = FALSE)

#数据???
expr_gene <- read.csv("expr_gene_data101_cpm.csv", row.names = 1)
names(expr_gene) <- gsub("\\.","-",names(expr_gene))
expr_gene <- expr_gene[which(names(expr_gene)%in%sample_choosed_20$Sample.ID)]
expr_gene <- as.data.frame(t(expr_gene))
expr_gene$class <- rownames(expr_gene)
expr_gene$class[which(grepl("-01",expr_gene$class))] <- "yes"
expr_gene$class[which(grepl("-11",expr_gene$class))] <- "no"
write.csv(expr_gene,file = "expr_gene_choosed101.csv")
expr_gene <- read.csv("expr_gene_choosed101.csv",row.names = 1)
expr_gene_log <- log10(subset(expr_gene, select = -class))
expr_gene_log$class <- expr_gene$class
write.csv(expr_gene_log,file = "expr_gene_log_choosed101.csv",row.names = FALSE)

sample_all_test <- read.csv("Sample_all_test.csv")
expr_gene_all <- read.csv("CPM101.csv", row.names = 1)
names(expr_gene_all) <- gsub("\\.","-",names(expr_gene_all))
expr_gene_all <- expr_gene_all[which(names(expr_gene_all)%in%sample_all_test$x)]
expr_gene_all <- as.data.frame(t(expr_gene_all))
expr_gene_all$class <- rownames(expr_gene_all)
expr_gene_all$class[which(grepl(".01",expr_gene_all$class))] <- "yes"
expr_gene_all$class[which(grepl(".11",expr_gene_all$class))] <- "no"
write.csv(expr_gene_all,file = "expr_gene_all101.csv")
expr_gene_all_log <- log10(subset(expr_gene_all, select = -class))
expr_gene_all_log$class <- expr_gene_all$class
write.csv(expr_gene_all_log,file = "expr_gene_all_log101.csv",row.names = FALSE)

save.image("class101.RData")
load("class101.RData")
