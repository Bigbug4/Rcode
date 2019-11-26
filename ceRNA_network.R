ceRNA network
```{r}
#mirna_target_gene <- read.csv("miRNAseq-2\\result_mirna_analysis\\starBase\\miRWalk_miRNA_Targets.csv")
#mirna_target_mrna <- read.csv("miRNAseq-2\\result_mirna_analysis\\starBase\\miRWalk_miRNA_Targets_gene2mirna_select.csv")
mirna_target_lncrna <- read.table("miRNAseq-2\\result_mirna_analysis\\starBase\\get_mirna_result_lncrna_v2.txt",sep = "\t",header = TRUE)
mirna_target_lncrna <- mirna_target_lncrna[,1:5]

mirna_target_mrna <- read.table("miRNAseq-2\\result_mirna_analysis\\starBase\\get_mirna_result_mrna_v2.txt",sep = "\t",header = TRUE)
mirna_target_mrna <- mirna_target_mrna[,1:5]

lncrna_target_rna <- read.table("miRNAseq-2\\result_mirna_analysis\\starBase\\get_rna_result_lncrna.txt",sep = "\t",header = TRUE)
lncrna_target_rna <- lncrna_target_rna[,1:5]
lncrna_target_rna <- lncrna_target_rna[-which(lncrna_target_rna$geneID=="geneID"),]

mrna_target_by_mirna <- read.table("miRNAseq-2\\result_mirna_analysis\\starBase\\get_mrna_result_by_mirna.txt",sep = "\t",header = TRUE)


diff_mrna <- read.table("diff_gene_protein_symbol_v2.txt")
diff_lncrna <- read.table("diff_gene_lncrna_symbol_v2.txt")
diff_mirna <- read.table("diff_gene_mirna.txt")
names(diff_mrna) <- "GeneName"
names(diff_lncrna) <- "GeneName"
names(diff_mirna) <- "GeneName"
diff_mirna <- gsub("-mir-",replacement = "-miR-",diff_mirna$GeneName)
write.table(diff_mirna,"diff_mirna_name.txt",row.names = FALSE,col.names = FALSE,quote = FALSE)

sig_gene <- read.csv("sig_gene_v3.csv")


mirna_total <- mirna_target_lncrna[which(mirna_target_lncrna$miRNAid%in%unique(mirna_target_mrna$miRNAid)),1:2]
mirna_total <- mirna_total[!duplicated(mirna_total),]
node_mirna2mrna <- mirna_target_mrna[which(mirna_target_mrna$miRNAid%in%mirna_total$miRNAid),]
node_mirna2mrna <- node_mirna2mrna[which(node_mirna2mrna$geneName%in%sig_gene$gene_name),]
mirna_total <- mirna_total[which(mirna_total$miRNAid%in%node_mirna2mrna$miRNAid),]
node_mirna2lncrna <- mirna_target_lncrna[which(mirna_target_lncrna$miRNAid%in%mirna_total$miRNAid),]


node_info_self <- rbind(node_mirna2mrna,node_mirna2lncrna)

ceRNA_mrna <- unique(node_mirna2mrna$geneName)
ceRNA_lncrna <- unique(node_mirna2lncrna$geneName)
ceRNA_mirna <- unique(node_info_self$miRNAname)
ceRNA_mirna <- data.frame(geneName=ceRNA_mirna,geneType=rep("miRNA",length(ceRNA_mirna)))
ceRNA_mrna <- data.frame(geneName=ceRNA_mrna,geneType=rep("protein_coding",length(ceRNA_mrna)))
ceRNA_lncrna <- data.frame(geneName=ceRNA_lncrna,geneType=rep("lncRNA",length(ceRNA_lncrna)))
ceRNA <- rbind(ceRNA_mrna,ceRNA_mirna,ceRNA_lncrna)
write.csv(ceRNA,file = "miRNAseq-2\\result_mirna_analysis\\starBase\\ceRNA.csv")


lncrna2mrna <- lncrna_target_rna[which(lncrna_target_rna$pairGeneName%in%diff_mrna$GeneName),]
lncrna2mrna <- lncrna2mrna[!duplicated(lncrna2mrna),]
node_trible <- lncrna2mrna[,c(5,2,3)]
names(node_trible) <- names(node_trible)[c(2,1,3)]

mirna2mrna_lncrna <- mirna_target_mrna[which(mirna_target_mrna$geneName%in%node_trible$geneName),]
mirna2mrna_lncrna <- mirna2mrna_lncrna[,c(4,2,5)]
names(mirna2mrna_lncrna)[2] <- "pairGeneName"
mirna2mrna_lncrna$geneType <- rep("miRNA",dim(mirna2rna_lncrna)[1])

node_trible <- rbind(node_trible,mirna2mrna_lncrna)



write.csv(node_info_self,file = "miRNAseq-2\\result_mirna_analysis\\starBase\\mirna_union_node_selected.csv",row.names = FALSE)
write.table(node_info_self[,c(2,4)],file = "miRNAseq-2\\result_mirna_analysis\\starBase\\mirna_union_node_v2.txt",row.names = FALSE,col.names = FALSE,quote = FALSE,sep = "\t")
write.table(node_trible,file = "miRNAseq-2\\result_mirna_analysis\\starBase\\node_trible_mrna2rna.txt",row.names = FALSE,col.names = FALSE,quote = FALSE,sep = "\t")
write.table(unique(node_trible$geneName),file = "miRNAseq-2\\result_mirna_analysis\\starBase\\node_trible_mrna.txt",row.names = FALSE,col.names = FALSE,quote = FALSE,sep = "\t")
write.table(unique(node_trible$pairGeneName),file = "miRNAseq-2\\result_mirna_analysis\\starBase\\node_trible_unmrna.txt",row.names = FALSE,col.names = FALSE,quote = FALSE,sep = "\t")



#
#
# ceRNA network verify
#

expr_data_choosed  <- read.csv("limma_result\\expr_data_choosed.csv",row.names = 1)
diff_lncrna <- result.lncrna[,c(1,9,10)]
expr_diff_lncrna <- expr_data_choosed[as.character(diff_lncrna$gene_id_v),]
row.names(expr_diff_lncrna) <- diff_lncrna$gene_name
expr_ce_lncrna <- expr_diff_lncrna[as.character(ceRNA_lncrna$geneName),]
expr_diff_mrna <- expr_data_choosed[as.character(diff_protein$gene_id_v),]
row.names(expr_diff_mrna) <- diff_protein$gene_name
expr_ce_mrna <- expr_diff_mrna[as.character(ceRNA_mrna$geneName),]

expr_cerna <- rbind(expr_ce_mrna,expr_ce_lncrna)

select_ceNetwork <- function(ce_lncrna=ceRNA_lncrna$geneName,ce_mrna=ceRNA_mrna$geneName){
  mrna_lncrna.mi_count <- matrix(nrow = length(ce_lncrna),ncol = length(ce_mrna))
  row.names(mrna_lncrna.mi_count) <- ce_lncrna
  colnames(mrna_lncrna.mi_count) <- ce_mrna
  mrna_lncrna.mi_name <- mrna_lncrna.mi_count
  mrna_lncrna.ceNet <- mrna_lncrna.mi_count
  for(mrna in 1:length(ce_mrna)){
    mirna_mrna <- node_info_self$miRNAname[which(as.character(node_info_self$geneName)==ce_mrna[mrna])]
    for (lncrna in 1:length(ce_lncrna)) {
      mirna_lncrna <- node_info_self$miRNAname[which(as.character(node_info_self$geneName)==ce_lncrna[lncrna])]
      inter_count <- length(intersect(mirna_mrna,mirna_lncrna))
      mrna_lncrna.mi_count[lncrna,mrna] <- inter_count
      mrna_lncrna.mi_name[lncrna,mrna] <- paste(intersect(mirna_mrna,mirna_lncrna),collapse = ",")
      correla <- cor(as.numeric(expr_ce_lncrna[as.character(ce_lncrna[lncrna]),]),as.numeric(expr_ce_mrna[as.character(ce_mrna[mrna]),]))
      mrna_lncrna.ceNet[lncrna,mrna] <- ifelse(!is.na(correla),correla,0)
      
      
    }
  }
  return(list(miRNA_counts=mrna_lncrna.mi_count,miRNA_name=mrna_lncrna.mi_name,ceNet = mrna_lncrna.ceNet))
}

ceNet_matrix <- select_ceNetwork(ce_lncrna=ceRNA_lncrna$geneName,ce_mrna=ceRNA_mrna$geneName)

ceNet_cor_matrix <- ceNet_matrix$ceNet

cut_null <- function(cor_matrix = ceNet_cor_matrix_cut,cutoff = 0.5){
  for(i in 1:dim(cor_matrix)[1]){
    for (j in 1:dim(cor_matrix)[2]) {
      cor_matrix[i,j] <- ifelse(cor_matrix[i,j]>cutoff,cor_matrix[i,j],0)
    }
  }
  for(i in dim(cor_matrix)[1]:1){
    if(max(abs(cor_matrix[i,]))==0){
      cor_matrix <- cor_matrix[-i,]
    }
  }
  for(j in dim(cor_matrix)[2]:1){
    if(max(abs(cor_matrix[,j]))==0){
      cor_matrix <- cor_matrix[,-j]
    }
  }
  return(cor_matrix)
}


ceRNA_lncrna_selected <- row.names(ceNet_cor_matrix_cut)
ceRNA_mrna_selected <- colnames(ceNet_cor_matrix_cut)
ceNet_symbol_select <- ceNet_matrix$miRNA_name[ceRNA_lncrna_selected,ceRNA_mrna_selected]


node_info <- read.csv("ceRNA\\mirna_union_node.csv")
node_info <- node_info[which(node_info$geneName%in%union(ceRNA_mrna_selected,ceRNA_lncrna_selected)),]


write.csv(ceNet_matrix$miRNA_counts,"ceRNA\\ce_matrix_mirna_count.csv")
write.table(ceNet_matrix$miRNA_name,"ceRNA\\ce_matrix_mirna_symbol.txt",quote = FALSE,sep = "\t")
write.csv(ceNet_matrix$ceNet,"ce_matrix_mirna_correlation.csv")


ceNet_cor_matrix_cut <- cut_null(cor_matrix = ceNet_cor_matrix,cutoff = 0.8)
ceNet_cor_matrix_cut_unregular <- ceNet_matrix$ceNet[ceRNA_lncrna_selected,ceRNA_mrna_selected]
write.csv(ceNet_cor_matrix_cut[ceRNA_lncrna_selected,ceRNA_mrna_selected],"ceRNA\\ce_matrix_mirna_correlation_cut_regular.csv",quote = FALSE)
write.csv(ceNet_cor_matrix_cut_unregular,"ceRNA\\ce_matrix_mirna_correlation_cut.csv")
write.csv(ceNet_symbol_select,"ceRNA\\ce_matrix_mirna_symbol_cut.csv",quote = FALSE)
## inner_mirna is the interact mirna between lncrna and mrna which were selected in ce_network_select.csv 

inner_mirna <- paste(ceNet_symbol_select,collapse = ",")
inner_mirna <- strsplit(inner_mirna,",")
inner_mirna <- inner_mirna[[1]]
inner_mirna <- inner_mirna[-which(inner_mirna=="")]
inner_mirna <- unique(inner_mirna)

write.table(inner_mirna,"ceRNA\\ce_network_mirna.txt",quote = FALSE)
inner_mirna <- read.table("ceRNA\\ce_network_mirna.txt")
ceNet_select_lncran_mrna <- read.csv("ceRNA\\ce_matrix_mirna_correlation_cut_regular.csv",row.names = 1)

node_info_select <- node_info_self[which(node_info_self$geneName%in%union(ceRNA_lncrna_selected,ceRNA_mrna_selected)&node_info_self$miRNAname%in%inner_mirna),]
ceRNA_lncrna_selected <- ceRNA_lncrna_selected[which(ceRNA_lncrna_selected%in%node_info_select$geneName)]
ceRNA_mrna_selected <- ceRNA_mrna_selected[which(ceRNA_mrna_selected%in%node_info_select$geneName)]


write.csv(node_info_select,"ceRNA\\ce_network_selected_node_info_all.csv",row.names = FALSE,quote = FALSE)


ceRNA_lncrna_selected_info <- diff_lncrna[which(diff_lncrna$gene_name%in%ceRNA_lncrna_selected),]
ceRNA_mrna_selected_info <- sig_gene[which(sig_gene$gene_name%in%ceRNA_mrna_selected),]
write.csv(ceRNA_lncrna_selected_info,"ceRNA\\ce_network_selected_lncRNA_info.csv")
write.csv(ceRNA_mrna_selected_info,"ceRNA\\ce_network_selected_mRNA_info.csv",row.names = FALSE)

