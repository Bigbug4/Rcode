
#DIR=getwd()

setwd("E:\\Rcode\\data")

# sample choosing by clinical.

clinical_all  <- read.csv("clinical.tsv",sep = "\t",header = TRUE)

## make samples grouped as 
## not reported, stage_1, stage_1a, stage_1b, stage_2, stage_2a, stage_2b, 
## stage_3, stage_3a, stage_3b, stage_3c, stage_4

#sample_stage_i <- clinical_all$submitter_id[which(clinical_all$tumor_stage=="stage i")]
sample_stage_ia <- clinical_all$submitter_id[which(clinical_all$tumor_stage=="stage ia")]
sample_stage_ib <- clinical_all$submitter_id[which(clinical_all$tumor_stage=="stage ib")]
#sample_stage_ii <- clinical_all$submitter_id[which(clinical_all$tumor_stage=="stage ii")]
sample_stage_iia <- clinical_all$submitter_id[which(clinical_all$tumor_stage=="stage iia")]
sample_stage_iib <- clinical_all$submitter_id[which(clinical_all$tumor_stage=="stage iib")]
#sample_stage_iii <- clinical_all$submitter_id[which(clinical_all$tumor_stage=="stage iii")]
sample_stage_iiia <- clinical_all$submitter_id[which(clinical_all$tumor_stage=="stage iiia")]
sample_stage_iiib <- clinical_all$submitter_id[which(clinical_all$tumor_stage=="stage iiib")]
sample_stage_iiic <- clinical_all$submitter_id[which(clinical_all$tumor_stage=="stage iiic")]
sample_stage_iv <- clinical_all$submitter_id[which(clinical_all$tumor_stage=="stage iv")]

stage_series <- ls()[2:9] # get stage info from environment

path_expr_data <- "sum.csv"
sum_expr_data <- read.csv(path_expr_data)
rownames(sum_expr_data) <- sum_expr_data$ENSGgeneID
names(sum_expr_data) <- gsub(".","-",names(sum_expr_data),fixed = TRUE)

# scale
library(edgeR)
expr_all <- sum_expr_data[-c(1:5),-1]
counts <- as.matrix(expr_all)
dge <- DGEList(counts = counts)
dge <- calcNormFactors(dge)
CPM <- cpm(dge, log=FALSE, prior.count=2)
logCPM <- cpm(dge, log=TRUE, prior.count=3)
write.csv(CPM,"CPM_all.csv")
write.csv(logCPM,"logCPM_all.csv")

# find the key gene in time series

sig_gene_top <- read.csv("sig_gene_v3.csv",header = TRUE)
row.names(sig_gene_top) <- sig_gene_top[,1]
sig_gene_top <- sig_gene_top[,-1]

all_sample <- read.csv("SampleID_and_Type.csv")
all_sample_tumor <- all_sample[which(all_sample$Sample.Type=="Primary Tumor"),]

sum_expr_data <- as.data.frame(CPM)
expr_gene_top <- sum_expr_data[which(row.names(sum_expr_data)%in%sig_gene_top$gene_id_v),which(names(sum_expr_data)%in%all_sample_tumor$Sample.ID)]
#expr_keeping_gene <- sum_expr_data[which(row.names(sum_expr_data)%in%house_keeping_gene$gene_id_v),which(names(expr_keeping_gene)%in%all_sample_tumor$Sample.ID)]

get_stage_matrix <- function(expr_data=expr_gene_top,stage_series=stage_series){
  stage_matrix <- data.frame(row.names = row.names(expr_data))
  for(stage in 1:length(stage_series)){
    case <- get(stage_series[stage])
    sample <- all_sample_tumor$Sample.ID[which(all_sample_tumor$Case.ID%in%case)]
    stage_matrix[,stage] <- rowMeans(expr_data[,which(names(expr_data)%in%sample)])
  }
  names(stage_matrix) <- stage_series
  return(stage_matrix)
}


expr_stage_series <- get_stage_matrix(expr_gene_top,stage_series)
expr_stage_series_scale <- expr_stage_series
#col_sum <- colSums(expr_stage_series)
row_sum <- rowSums(expr_stage_series)

stage_info <- rep(0,8)
for(stage in 1:length(stage_series)){
  case <- get(stage_series[stage])
  s <- all_sample_tumor$Sample.ID[which(all_sample_tumor$Case.ID%in%case)]
  stage_info[stage] <- length(s)
}
rm(case)
names(stage_info) <- stage_series
stage_info <- data.frame(count = stage_info)
write.csv(stage_info,file = "stage_info.csv")

for(i in 1:55){
  for(j in 1:8){
    #expr_stage_series_scale[i,j] <- expr_stage_series[i,j]/col_sum[j]
    expr_stage_series_scale[i,j] <- expr_stage_series[i,j]/row_sum[i]
    
  }
}

# expr_stage_series_scale_diff <- expr_stage_series_scale[,1:7]
# for(j in 1:8){
#   for(i in 1:55){
#     expr_stage_series_scale_diff[i,j] <- expr_stage_series_scale[i,j]/expr_stage_series_scale[i,1]
#     expr_stage_series_scale_diff[i,j] <- log2(expr_stage_series_scale_diff[i,j])
#     
#     
#   }
#   names(expr_stage_series_scale_diff)[j] <- paste(names(expr_stage_series_scale)[1],names(expr_stage_series_scale)[j],sep = "-")
#   names(expr_stage_series_scale_diff) <- gsub("sample_","",names(expr_stage_series_scale_diff))
# }
# 
# expr_stage_series_scale_diff_1 <- expr_stage_series_scale[,1:7]
# for(j in 1:7){
#   for(i in 1:55){
#     expr_stage_series_scale_diff_1[i,j] <- expr_stage_series_scale[i,j+1]/expr_stage_series_scale[i,j]
#     expr_stage_series_scale_diff_1[i,j] <- log2(expr_stage_series_scale_diff_1[i,j])
#     
#   }
#   names(expr_stage_series_scale_diff_1)[j] <- paste(names(expr_stage_series_scale)[j],names(expr_stage_series_scale)[j+1],sep = "-")
#   names(expr_stage_series_scale_diff_1) <- gsub("sample_","",names(expr_stage_series_scale_diff_1))
# }
# 


write.csv(expr_stage_series,file = "expr_stage_series_v3.csv")
write.csv(expr_stage_series_scale,file = "expr_stage_series_scale_v3.csv")

# write.csv(expr_stage_series_scale_diff,file = "expr_stage_series_diff_v3.csv")
# write.csv(expr_stage_series_scale_diff_1,file = "expr_stage_series_scale_diff1_v3.csv")
# 

#
# find key gene in different stages.
#
expr_stage_series_scale <- read.csv("expr_stage_series_scale_v3.csv",header = TRUE)
row.names(expr_stage_series_scale) <- expr_stage_series_scale[,1]
expr_stage_series_scale[,1] <- NULL

key_gene_in_diff_stage <- function(expr_stage_series_scale = expr_stage_series_scale){
  v_sig_stage <- rep(0,nrow(expr_stage_series_scale))
  var_matrix <- expr_stage_series_scale
  for(gene in 1:dim(expr_stage_series_scale)[1]){
    for(stage in 1:dim(expr_stage_series_scale)[2]){
      v <- as.numeric(expr_stage_series_scale[gene,-stage])
      var_matrix[gene,stage] <- var(v/mean(v))
    }
    v_sig_stage[gene] <- names(var_matrix)[which(var_matrix[gene,]==min(var_matrix[gene,]))]
  }
  expr_stage_series_var <- var_matrix
  expr_stage_series_var$sig_stage <- v_sig_stage
  return(expr_stage_series_var)
}

expr_stage_series_var <- key_gene_in_diff_stage(expr_stage_series_scale)
write.csv(expr_stage_series_var,file = "expr_stage_series_var_v3.csv")


expr_stage_series_var <- read.csv("expr_stage_series_var_v3.csv",header = TRUE,row.names = 1)

expr_series_var_1 <- as.matrix(expr_stage_series_var[,1:8])
expr_stage_series_scale_1 <- as.matrix(expr_stage_series_scale)
expr_stage_series_1 <- as.matrix(expr_stage_series)

sig_gene <- sig_gene_top
row.names(sig_gene) <- sig_gene$gene_id_v
sig_gene <- sig_gene[as.character(row.names(expr_series_var_1)),]

colnames(expr_series_var_1) <- gsub("sample_","",colnames(expr_series_var_1))
row.names(expr_series_var_1) <- sig_gene$gene_name

colnames(expr_stage_series_scale_1) <- gsub("sample_","",colnames(expr_stage_series_scale_1))
row.names(expr_stage_series_scale_1) <- sig_gene$gene_name

colnames(expr_stage_series_1) <- gsub("sample_","",colnames(expr_stage_series_1))
row.names(expr_stage_series_1) <- sig_gene$gene_name

write.csv(expr_stage_series_scale_1,file = "expr_stage_series_scale_a3.csv")
write.csv(expr_series_var_1,file = "expr_stage_series_var_a3.csv")
write.csv(expr_stage_series_1,file = "expr_stage_series_a3.csv")

expr_series_var_2 <- expr_series_var_1

for(i in 1:dim(expr_series_var_2)[1]){
  for (j in 1:dim(expr_series_var_2)[2]) {
    expr_series_var_2[i,j] <- ifelse(expr_series_var_2[i,j]==min(expr_stage_series_var[i,1:8]),0,1)
    
  }
}

write.csv(expr_series_var_2,file = "expr_series_var_2.csv")


for(i in 1:dim(expr_series_var_1)[1]){
  for (j in 1:dim(expr_series_var_1)[2]) {
    expr_series_var_1[i,j] <- ifelse(expr_series_var_1[i,j]==min(expr_stage_series_var[i,1:8]),-1,1)
    if(expr_series_var_1[i,j]==-1){
      conf.int <- t.test(expr_stage_series_scale_1[i,-j])$conf.int
      if(expr_stage_series_scale_1[i,j]<conf.int[1]||expr_stage_series_scale_1>conf.int[2]){
        expr_series_var_1[i,j] <- 1
      }
      else if(expr_stage_series_scale_1[i,j]>=mean(expr_stage_series_scale_1[i,])){
        expr_series_var_1[i,j] <- 0
      }
    }
  }
}


signal_gene <- expr_series_var_1[which(rowMeans(expr_series_var_1)<1),]
write.csv(signal_gene,file = "signal_gene.csv")
time_sig_gene <- sig_gene[which(sig_gene$gene_name%in%row.names(signal_gene)),]

st1 <- 55 - colSums(expr_series_var_2)
pt1 <- st1/55
st2 <- 8 -colSums(signal_gene)
pt2 <- st2/8

df <- rbind(st1,pt1,st2,pt2)
write.csv(df,"time_gene_results.csv")

# save.image("time_sig_gene.RData")
# load("time_sig_gene.RData")
