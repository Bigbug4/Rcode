# 分类

DIR=getwd()
setwd("E:\\Rcode\\data")

expr_data_cpm <- read.csv("CPM.csv",row.names = 1)
names(expr_data_cpm) <- gsub("\\.","-",names(expr_data_cpm))
sample_list <- read.csv("sample_list.csv")
sample_list <- as.vector(t(sample_list))
expr_data_choosed_cpm <- expr_data_cpm[,which(names(expr_data_choosed_cpm)%in%sample_list)]
write.csv(expr_data_choosed_cpm,file = "expr_gene_data_cpm.csv")

sample_choosed <- read.csv("sample_choosed.csv")
case_choosed <- summary(sample_choosed$Case.ID)
case_choosed <- names(case_choosed)[which(case_choosed==2)]
s <- sample_choosed[which(sample_choosed$Case.ID%in%case_choosed),]
s <- s[order(s$Case.ID),]
sample_choosed_20 <- s[1:20,]


sig_gene <- read.csv("sig_gene_v3.csv")
sig_gene <- sig_gene[order(sig_gene$logFC,decreasing = TRUE),]
sig_gene_up <- sig_gene[which(sig_gene$logFC>0),]
sig_gene_down <- sig_gene[which(sig_gene$logFC<0),]

comp_in_ssample <- function(expr_up,expr_down){
  n <- dim(expr_down)[2]
  comp_vector <- rep(0,n)
  names(comp_vector) <- names(expr_down)
  for (i in 1:n) {
    comp_vector[i] <- sum(expr_up[,i])/sum(expr_down[,i])
    
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


s_combine1 <- combine_sample(expr_data_choosed_cpm,sample_choosed_20,10)
s1 <- comp_in_ssample(log10(s_combine1[1:41,]),log10(s_combine1[42:55,]))
c_factor1 <- exp(sum(log(s1))/length(s1))

# 正负样本
expr_cpm_log <- log10(expr_data_cpm)
expr_cpm_log <- as.data.frame(t(expr_cpm_log))
expr_cpm_log$class <- rownames(expr_cpm_log)
expr_cpm_log$class[which(grepl("-01",expr_cpm_log$class))] <- "yes"
expr_cpm_log$class[which(grepl("-11",expr_cpm_log$class))] <- "no"
write.csv(expr_cpm_log,"expr_cpm_log.csv",row.names = FALSE)

sample_all_test <- read.csv("Sample_all_test.csv",header = T)

s3 <- comp_in_ssample(log10(expr_data_cpm[1:41,]),log10(expr_data_cpm[42:55,]))
s3 <- s3[order(s3)]
s3 <- s3[which(names(s3)%in%sample_all_test)]
write.csv(s3,"test_cpm_result.csv")
Positive <- s3[which(grepl("-01",names(s3)))]
Negative <- s3[which(grepl("-11",names(s3)))]

TP1 <- Positive[Positive>c_factor1]
TN1 <- Negative[Negative<c_factor1]
ACC1 <- (length(TP1)+length(TN1))/length(s3)
sensitive1 <- length(TP1)/length(Positive)
precision1 <- length(TN1)/length(Negative)

save.image("classification_c.Rdata")
# load("classification_c.Rdata")

#数据???
expr_gene <- read.csv("expr_gene_data_cpm.csv", row.names = 1)
names(expr_gene) <- gsub("\\.","-",names(expr_gene))
expr_gene <- expr_gene[which(names(expr_gene)%in%sample_choosed_20$Sample.ID)]
expr_gene <- as.data.frame(t(expr_gene))
expr_gene$class <- rownames(expr_gene)
expr_gene$class[which(grepl("-01",expr_gene$class))] <- "yes"
expr_gene$class[which(grepl("-11",expr_gene$class))] <- "no"
write.csv(expr_gene,file = "expr_gene_choosed_cpm.csv")
expr_gene <- read.csv("expr_gene_choosed_cpm.csv",row.names = 1)
expr_gene_log <- log10(subset(expr_gene, select = -class))
expr_gene_log$class <- expr_gene$class
write.csv(expr_gene_log,file = "expr_gene_log_choosed_cpm.csv",row.names = FALSE)

sample_all_test <- read.csv("Sample_all_test.csv")
expr_gene_all <- read.csv("CPM.csv", row.names = 1)
names(expr_gene_all) <- gsub("\\.","-",names(expr_gene_all))
expr_gene_all <- expr_gene_all[which(names(expr_gene_all)%in%sample_all_test$x)]
expr_gene_all <- as.data.frame(t(expr_gene_all))
expr_gene_all$class <- rownames(expr_gene_all)
expr_gene_all$class[which(grepl(".01",expr_gene_all$class))] <- "yes"
expr_gene_all$class[which(grepl(".11",expr_gene_all$class))] <- "no"
write.csv(expr_gene_all,file = "expr_gene_all_cpm.csv")
expr_gene_all_log <- log10(subset(expr_gene_all, select = -class))
expr_gene_all_log$class <- expr_gene_all$class
write.csv(expr_gene_all_log,file = "expr_gene_all_log_cpm.csv",row.names = FALSE)
