DIR=getwd()
setwd("E:\\Rcode\\data")

# Visualization

library(ggplot2)

result <- read.csv("result_annotation_v2.csv")
result.protein <- result[ which(result$gene_biotype=='protein_coding'),]
result.protein <- result.protein[,c(2,3,5,6,12)]
result.protein$gene_biotype <- rep('protein',nrow(result.protein))
result.lncRNA <- result[ which(result$gene_biotype%in%c("lincRNA","non coding","3prime_overlapping_ncRNA","antisense","retained_intron","sense_intronic","sense_overlapping","macro_lncRNA")),]
result.lncRNA <- result.lncRNA[,c(2,3,5,6,12)]
result.lncRNA$gene_biotype <- rep('lncRNA',nrow(result.lncRNA))
result.miRNA <- read.csv("result.mirna.csv")
result.miRNA <- result.miRNA[,c(2,3,5,6)]
result.miRNA$gene_biotype <- rep('miRNA',nrow(result.miRNA))

# different gene show
result.volcano <- rbind(result.protein,result.lncRNA,result.miRNA)
write.csv(result.volcano,'volcano.csv')
result.volcano <- read.csv('volcano.csv',row.names = 1)
significant <- ifelse(result.volcano$adj.P.Val<0.01 & abs(result.volcano$logFC)>= 1,ifelse(result.volcano$logFC > 1,'Up','Down'),'Nosignificant')
v_color <- c(Up = "red",Nosignificant = "gray",Down = "blue")
ggplot(result.volcano, aes(result.volcano$logFC, -log10(result.volcano$adj.P.Val), shape=factor(gene_biotype), col = significant)) +
  geom_point() +
  scale_color_manual(values = v_color) +
  labs(title = "Differential Expression",x="log2 (fold change)",y="-log10 (adj.P.Value)") +
  geom_hline(yintercept = -log10(0.01), lty=4,col="grey",lwd=0.6) +
  geom_vline(xintercept = c(-1, 1), lty=4,col="grey",lwd=0.6) + theme(plot.title = element_text(hjust = 0.5), legend.position="right", legend.title = element_blank(),panel.grid = element_blank())

# heatmap

# pheatmap::pheatmap(expr_diff_protein[1:100,sample_list],scale = "row",cluster_cols = FALSE)
expr_data_choosed <- read.csv("expr_gene_data.csv",row.names = 1)
names(expr_data_choosed) <- gsub(".","-",names(expr_data_choosed),fixed = TRUE)
# expr_data_choosed_cpm <- read.csv("expr_gene_data_cpm.csv",row.names = 1)
# names(expr_data_choosed_cpm) <- gsub(".","-",names(expr_data_choosed_cpm),fixed = TRUE)
sample_list <- read.csv("sample_list.csv")
sample_list <- as.vector(t(sample_list))
# heatmap(log10(as.matrix(expr_data_choosed[1:55, sample_list])))
# heatmap(log10(as.matrix(expr_data_choosed_cpm[1:55, sample_list])))
pheatmap::pheatmap(expr_data_choosed[1:55, sample_list],scale = "row",cluster_cols = FALSE,fontsize_row = 8, fontsize_col = 8)
# pheatmap::pheatmap(as.matrix(expr_data_choosed_cpm[1:55, sample_list]),scale = "row",cluster_cols = FALSE,fontsize_row = 8, fontsize_col = 8)

expr_stage_series <- read.csv("expr_stage_series_a3.csv",row.names = 1)
pheatmap::pheatmap(expr_stage_series[1:55, ],scale = "row",cluster_cols = FALSE,fontsize_row = 8, fontsize_col = 8)

expr_stage_series_var <- read.csv("expr_stage_series_var_a3.csv",row.names = 1)
pheatmap::pheatmap(expr_stage_series_var[1:55, ],scale = "row",cluster_rows = FALSE,cluster_cols = FALSE,fontsize_row = 8, fontsize_col = 8)

expr_stage_gene <- read.csv("expr_series_var_2.csv",row.names = 1)
pheatmap::pheatmap(expr_stage_gene[1:55, ],cluster_rows = FALSE,cluster_cols = FALSE,color = colorRampPalette(c("#2075FF","white"))(50),scale = "none",fontsize_row = 8, fontsize_col = 8)

signal_gene <- read.csv("signal_gene.csv",row.names = 1)
pheatmap::pheatmap(signal_gene[1:7, ],cluster_rows = FALSE,cluster_cols = FALSE,color = colorRampPalette(c("#2075FF","white"))(50),scale = "none",fontsize_row = 8, fontsize_col = 8)

# bar
sig_gene <- read.csv("sig_gene_v3.csv")
barplot(sig_gene$score[order(sig_gene$score,decreasing = TRUE)]-1,names.arg = sig_gene$geneName,
        horiz=TRUE, xlim = c(0,1),las=1,border="white",col = '#ff8800',xlab='SCORE - 1',cex.names=0.8)

grid(ny = NA)

# ROC
p <- read.csv("roc_test.csv")
attach(p)
# library(pROC)
# roc(roc_plot2, roc_plot1, plot=TRUE, print.thres=TRUE, print.auc=TRUE)

library(ROCR)
pred <- prediction(roc_plot2, roc_plot1) 
perf <- performance(pred,"tpr","fpr")
perf.auc<- performance(pred, measure = 'auc', x.measure = 'cutoff')
auc <- performance(pred,'auc')
auc = unlist(slot(auc,"y.values"))
df<- data.frame(x = attributes(perf)$x.values[[1]],y = attributes(perf)$y.values[[1]])

plot(perf,
     xlim=c(0,1), ylim=c(0,1),col='red', 
     main=paste("ROC curve (", "AUC = ",auc,")"),
     lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)

abline(0,1)

library(ggplot2)

ggplot(data = df) + geom_line(aes(x,y),colour = "#2E9FDF",size = 0.5) +   
  labs(title = paste("ROC curve (", "AUC = ",auc,")")) + 
  xlab("1 - Specificity") + 
  ylab("Sensitivity") + 
  theme(plot.title = element_text(size = 15)) +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), size = 0.5, color="grey", linetype="dashed") + theme_minimal()

load("df_RF.RData")
df1 <- df
names(df1)<- c("x1","y1")
load("df_SVM.RData")
df2 <- df
names(df2)<- c("x2","y2")
load("df_NB.RData")
df3 <- df
names(df3)<- c("x3","y3")

add_row <- function(df,k){
  for (i in 1:k){
    df <- rbind(df,c(1,1))
  }
  return(df)
}

df1 <- add_row(df1,3)
df3 <- add_row(df3,24)

df <- cbind(df1,df2,df3)

p1 <- ggplot(data = df) + geom_line(aes(x1,y1,colour = "#0F0"),size = 0.5,linetype=1) +
  labs(title = "ROC curve") + 
  xlab("1 - Specificity") + 
  ylab("Sensitivity") + 
  theme(plot.title = element_text(size = 15)) +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), size = 0.5, color="grey", linetype="dashed") + theme_minimal()

p2 <- geom_line(aes(x2,y2,colour = "#F00"),size = 0.5,linetype=1) 
p3 <- geom_line(aes(x3,y3,colour = "#00F"),size = 0.5,linetype=1) 

p1+p2+p3
save.image("visualization.RData")
# load("visualization.RData")
