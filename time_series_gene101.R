
DIR=getwd()

setwd("E:\\Rcode\\data")

# sample choosing by clinical.


stage_series <- read.csv("stage_info.csv")


path_expr_data <- "CPM_all.csv"
sum_expr_data <- read.csv(path_expr_data)
rownames(sum_expr_data) <- sum_expr_data$ENSGgeneID
names(sum_expr_data) <- gsub(".","-",names(sum_expr_data),fixed = TRUE)

# find the key gene in time series

sig_gene_top <- read.csv("sig_gene101.csv",header = TRUE,row.names = 1)

all_sample <- read.csv("SampleID_and_Type.csv")
all_sample_tumor <- all_sample[which(all_sample$Sample.Type=="Primary Tumor"),]

expr_gene_top <- sum_expr_data[which(row.names(sum_expr_data)%in%sig_gene_top$gene_id_v),which(names(sum_expr_data)%in%all_sample_tumor$Sample.ID)]
#expr_keeping_gene <- sum_expr_data[which(row.names(sum_expr_data)%in%house_keeping_gene$gene_id_v),which(names(expr_keeping_gene)%in%all_sample_tumor$Sample.ID)]

get_stage_matrix <- function(expr_data=expr_gene_top,stage_series=stage_series){
  stage_matrix <- data.frame(row.names = row.names(expr_data))
  for(stage in 1:length(stage_series)){
    case <- get(stage_series[stage])
    sample <- all_sample_tumor_stomach$Sample.ID[which(all_sample_tumor_stomach$Case.ID%in%case)]
    stage_matrix[,stage] <- rowMeans(expr_data[,which(names(expr_data)%in%sample)])
  }
  names(stage_matrix) <- stage_series
  return(stage_matrix)
}


expr_stage_series <- get_stage_matrix(expr_gene_top,stage_series)
expr_stage_series_scale <- expr_stage_series/colSums(expr_stage_series)
col_sum <- colSums(expr_stage_series)


for(i in 1:55){
  for(j in 1:8){
    expr_stage_series_scale[i,j] <- expr_stage_series[i,j]/col_sum[j]
    
  }
}



write.csv(expr_stage_series,file = "expr_stage_series101.csv")
write.csv(expr_stage_series_scale,file = "expr_stage_series_scale101.csv")

# find key gene in different stages.
#
expr_stage_series_scale <- read.csv("expr_stage_series_scale101.csv",header = TRUE, row.names = 1)

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
write.csv(expr_stage_series_var,file = "expr_stage_series_var101.csv")
