
DIR=getwd()
setwd("E:\\Rcode\\data")

# correlation matrix

sample_list <- names(expr_mirna)[order(group_list)]
expr_diff_protein <- expr_data_choosed[result.protein$gene_id_v,]
expr_diff_protein.annotation <- gene_id2name(result = expr_diff_protein)
expr_diff_protein <- cbind(expr_diff_protein,expr_diff_protein.annotation$gene_name)
row.names(expr_diff_protein) <- expr_diff_protein$`expr_diff_protein.annotation$gene_name`
expr_diff_protein <- expr_diff_protein[,sample_list]


expr_diff_lincrna <- expr_data_choosed[result.lincrna$gene_id_v,]
expr_diff_lincrna <- expr_diff_lincrna[,sample_list]

expr_diff_lncrna <- expr_data_choosed[result.lncrna$gene_id_v,]
expr_diff_lncrna <- expr_diff_lncrna[,sample_list]
write.csv(expr_diff_lncrna,"expr_diff_lncrna.csv")
write.csv(expr_data_choosed,"expr_data_choosed.csv")


# expr_diff_snorna <- expr_data_choosed[result.snorna$gene_id_v,]
# expr_diff_snrna <- expr_data_choosed[result.snrna$gene_id_v,]

expr_diff_mirna <- expr_mirna[row.names(diff_mirna),]
expr_diff_mirna <- expr_diff_mirna[,sample_list]

correlation_matrix <- function(expr_data1,expr_data2,pvalue = 0.05){
  gene_count1 <- dim(expr_data1)[1]
  gene_count2 <- dim(expr_data2)[1]
  cor_matrix <- matrix(data = 0,nrow = gene_count1,ncol = gene_count2)
  pvalue_matrix <- cor_matrix
  row.names(cor_matrix) <- row.names(expr_data1)
  colnames(cor_matrix) <- row.names(expr_data2)
  for (gene1 in 1:gene_count1) {
    for (gene2 in 1:gene_count2) {
      cor_result <- cor.test(as.numeric(expr_data1[gene1,]),as.numeric(expr_data2[gene2,]))
      #   if(cor_result$p.value<pvalue){
      #   cor_matrix[gene1,gene2] <- cor_result$estimate
      #   }
      # else{
      #   cor_matrix[gene1,gene2] <- -2
      # }
      cor_matrix[gene1,gene2] <- cor_result$estimate
      pvalue_matrix[gene1,gene2] <- cor_result$p.value
    }
  }
  return(list(cormatrix = cor_matrix,pvalue_matrix = pvalue_matrix))
}


correlation_matrixofgene <- correlation_matrix(expr_data = expr_diff_up_protein)
cor_mat_mirna <- correlation_matrix(expr_diff_mirna)

cor_mat_mirna_protein <- correlation_matrix(expr_diff_mirna,expr_diff_protein)
cor_mat_lincrna_protein <- correlation_matrix(expr_diff_lincrna,expr_diff_protein)
cor_mat_lncrna_protein <- correlation_matrix(expr_diff_lncrna,expr_diff_protein)


# find the key genes by give the signal entropy.

node_info <- diff_protein$logFC
names(node_info) <- diff_protein$gene_name

# m1 is cor_mat_mirna_protein, m2 is cor_mat_lincrna_protein, node_info denotes the degree of gene in the PPI network.
# col is gene,row is order rna.
# cor_mat$["cormatrix","pvalue_matrix"]

estimate_score <- function(cor_mat,score_list,pvalue){
  for(gene in 1:dim(cor_mat$cormatrix)[2]){
    gene_cor_vector <- abs(cor_mat$cormatrix[-which(cor_mat$pvalue_matrix>pvalue),gene])
    gene_cor_vector <- gene_cor_vector/sum(gene_cor_vector)
    #score_list[gene] <- -log10(1-sum(-gene_cor_vector*log(gene_cor_vector))/log(length(gene_cor_vector)))
    score_list[gene] <- sum(-gene_cor_vector*log(gene_cor_vector))
    
  }
  return(score_list)
}



keygenes_find <- function(node_info,m1,m2,pvalue=0.05){
  gene_counts <- length(node_info)
  score_list <- rep(0,gene_counts)
  names(score_list) <- names(node_info)
  score_list1 <- score_list
  score_list2 <- score_list
  
  score_list1 <- estimate_score(m1,score_list1,pvalue = pvalue)
  score_list2 <- estimate_score(m2,score_list2,pvalue = pvalue)
  score_list <- (score_list1+score_list2)
  
  return(score_list)
  
}

# names of vector node_info is different gene selected

estimate_score_with_PPI <- function(node_network,node_info = node_info){
  score_list_PPI <- rep(0,length(node_info))
  names(score_list_PPI) <- names(node_info)
  node_select <- unique(node_network$node1)
  for (edge in 1:length(node_select)){
    gene <- node_select[edge]
    score_list_PPI[gene] <- sum(node_network[which(node_network$node1==gene),]$combined_score)
  }
  score_list_PPI <- (score_list_PPI+1)
  return(score_list_PPI)
}


protein_interaction <- read.csv("string_protein_interactions_170.tsv",sep = "\t")
names(protein_interaction)[1] <- "node1"
node_network <- protein_interaction[,c("node1","combined_score")]

score_list_rna <- keygenes_find(node_info,cor_mat_mirna_protein,cor_mat_lncrna_protein)
score_list_PPI <- estimate_score_with_PPI(node_network,node_info)

#score_list_rna[order(score_list_rna,decreasing = TRUE)]
plot(score_list_rna[order(score_list_rna,decreasing = TRUE)])

score_result <- score_list_rna/max(abs(score_list_rna))+score_list_PPI/max(abs(score_list_PPI))

score_result_order <- score_result[order(score_result,decreasing = TRUE)]

x <- mutated_genes$Gene
y <- CNA_genes$Gene
var_gene <- score_result[union(x,y)]
var_gene <- var_gene[order(var_gene,decreasing = TRUE)]

score_result_order <- data.frame(score = score_result_order,row.names = names(score_result_order))
score_result_order$variation <- ifelse(row.names(score_result_order)%in%names(var_gene),1,0)
sig_gene <- score_result_order

write.csv(score_list_rna,file = "key_gene_score.csv")
write.csv(sig_gene,file = "sig_gene_mrna_v3.csv")


# in different case
sig_gene_top20 <- head(sig_gene,20)
sig_gene_tail20 <- tail(sig_gene,20)
sig_gene_top40 <- rbind(sig_gene_top20,sig_gene_tail20)


#sig_gene_top40 <- cbind(sig_gene_top40,xx)
#write.csv(sig_gene_top40,file = "sig_gene_top40_v2.csv")


sig_gene_v3 <- score_result_order[score_result_order$score>1.2,]
xx <- diff_protein
row.names(xx) <- xx$gene_name
xx <- xx[row.names(sig_gene_v3),]
sig_gene_v3 <- cbind(sig_gene_v3,xx)
write.csv(sig_gene_v3,file = "sig_gene_v3.csv")
