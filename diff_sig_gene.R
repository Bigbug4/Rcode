# 选取样本

DIR=getwd()
setwd("E:\\Rcode\\data")

path_expr_data <- "sum.csv"
sum_expr_data <- read.csv(path_expr_data)
rownames(sum_expr_data) <- sum_expr_data$ENSGgeneID
names(sum_expr_data) <- gsub(".","-",names(sum_expr_data),fixed = TRUE)
          
# case <- read.csv("SampleID_and_Type.csv")
# case_normal <- case[case$Sample.Type=="Solid Tissue Normal",]
# sample_normal <- case_normal$Sample.ID
# case_choosed <- case_normal$Case.ID
# sample_choosed <- case[case$Case.ID%in%case_choosed,]
write.csv(sample_choosed,"sample_list.csv")

sample_choosed <- read.csv("sample_list.csv")
sample_normal <- sample_choosed[,1:20]
names(sample_choosed) <- "Sample.ID"
sample_tumor <- sample_choosed$Sample.ID[!sample_choosed$Sample.ID%in%sample_normal]
expr_data_choosed <- sum_expr_data[,which(names(sum_expr_data)%in%sample_choosed$Sample.ID)]

names(expr_data_choosed) <- gsub(".",replacement = "-",names(expr_data_choosed),fixed = TRUE)

expr_data_normal <- sum_expr_data[,which(names(sum_expr_data)%in%sample_normal)]
expr_data_tumor <- sum_expr_data[,which(names(sum_expr_data)%in%sample_tumor)]

expr_data_normal$mean <- rowMeans(expr_data_normal)
expr_data_tumor$mean <- rowMeans(expr_data_tumor)

expr_data_normal <- expr_data_normal[-c(1:5),]
expr_data_tumor <- expr_data_tumor[-c(1:5),]

# 清除表达量为零的部分
cut_empty <- function(expr_data){
  expr_cut <- expr_data
  for(i in dim(expr_cut)[1]:1)
  {
    if(expr_data$mean[i]==0){
      expr_cut <- expr_cut[-i,]
    }
  }
  return(expr_cut)
}

# 双样本t检
double_sample_t_test <- function(expr_normal,expr_tumor,pv = 0.001){
  n <- dim(expr_normal)[1]
  p_value <- rep(0,n)
  change <- rep(0,n)
  for(gene in 1:dim(expr_normal)[1]){
    gene_normal <- expr_normal[gene,1:22]
    gene_tumor <- expr_tumor[gene,1:20]
    r <- t.test(gene_normal,gene_tumor,var.equal = FALSE)
    p_value[gene] <- r$p.value
    if(p_value[gene]>=pv)
    {
      next
    }
    else{
      if(r$estimate[1]<r$estimate[2]){
        change[gene] <- 1
      }
      else{
        change[gene] <- -1
      }
    }
    
  }
  return(data.frame(p_value,change,row.names = row.names(expr_normal)))
}

expr_data_normal_cut <- cut_empty(expr_data = expr_data_normal)
expr_data_tumor_cut <- cut_empty(expr_data = expr_data_tumor)
diff_result <- double_sample_t_test(expr_data_normal_cut,expr_data_tumor_cut)
diff_result$gene_id_v <- row.names(diff_result)


# ensembel_gene <- read.csv("ncRNA_annotation\\ensembel_gene_all.csv",header = TRUE)
# diff_result <- merge(diff_result,ensembel_gene)
# diff_result_gene <- diff_result[,c(1,5,6,2,3,7)]
# 
# write.csv(diff_result_gene,file = "diff_expression_gene.csv",row.names = FALSE)
# write.csv(diff_result_gene[which(diff_result_gene$change>0),],file = "diff_expression_gene_up.csv",row.names = FALSE)
# write.csv(diff_result_gene[which(diff_result_gene$change<0),],file = "diff_expression_gene_down.csv",row.names = FALSE)
# 
# write.csv(diff_result_gene[which(diff_result_gene$change!=0),],file = "diff_expression_gene_change.csv",row.names = FALSE)
# write.csv(sample_choosed,file = "sample_choosed.csv",row.names = FALSE)


