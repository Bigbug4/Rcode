setwd("C:\\Users\\Administrator\\Desktop")

expr_data_annotation <- read.csv("ensembel_gene_all.csv",row.names = 1)
expr_data_choosed1 <- read.table("差异表达的mRNAs1.txt",sep = "\t",header = T)
expr_data_choosed2 <- read.table("差异表达的mRNAs2.txt",sep = "\t",header = T)

expr_data_annotation$gene_name<-as.character(expr_data_annotation$gene_name)
expr_data_choosed1$Name=as.character(expr_data_choosed1$Name)
expr_data_choosed2$Name=as.character(expr_data_choosed2$Name)

expr_data_choosed.annotation1 <- expr_data_annotation[which(expr_data_choosed1$Name%in%expr_data_annotation$gene_name),]
expr_data_choosed.annotation2 <- expr_data_annotation[which(expr_data_choosed2$Name%in%expr_data_annotation$gene_name),]
expr_data_choosed <- read.table("seq1.txt",sep = "\t",header = T,row.names = 2)
