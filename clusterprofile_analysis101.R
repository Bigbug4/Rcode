
DIR=getwd()
setwd("E:\\Rcode\\data")

# clusterprofile analysis

library(clusterProfiler)
library(org.Hs.eg.db)

result.sig <- read.csv("sig_gene_101.csv",row.names = 1)

enrichGO.sig_gene.CC <- enrichGO(gene = result.sig$gene_name,OrgDb="org.Hs.eg.db",keyType = "SYMBOL",pvalueCutoff = 0.05,qvalueCutoff = 0.05,ont = "CC",pAdjustMethod = "BH")
enrichGO.sig_gene <- enrichGO(gene = result.sig$gene_name,OrgDb="org.Hs.eg.db",keyType = "SYMBOL",pvalueCutoff = 0.05,qvalueCutoff = 0.2,ont = "all",pAdjustMethod = "BH")
enrichGO.sig_gene.BP <- enrichGO(gene = result.sig$gene_name,OrgDb="org.Hs.eg.db",keyType = "SYMBOL",pvalueCutoff = 0.05,qvalueCutoff = 0.05,ont = "BP",pAdjustMethod = "BH")
enrichGO.sig_gene.MF <- enrichGO(gene = result.sig$gene_name,OrgDb="org.Hs.eg.db",keyType = "SYMBOL",pvalueCutoff = 0.05,qvalueCutoff = 0.05,ont = "MF",pAdjustMethod = "BH")

dotplot(enrichGO.sig_gene,showCategory=10,title="Enrichment GO Top10") #泡泡图
barplot(enrichGO.sig_gene, showCategory=10,title="Enrichment GO Top10")  #柱状图
# plotGOgraph(enrichGO.sig_gene) 	#GO图，看不清楚可以尝试左上角另存为pdf
emapplot(enrichGO.sig_gene,showCategory = 10)
cnetplot(enrichGO.sig_gene,showCategory = 10,foldChange = result.sig$logFC,circular = TRUE,colorEdge = TRUE)

sig_gene_go_term <- enrichGO.sig_gene@result[,c("Description","geneID")]


# kegg
diff_gene_protein_entrezid <- bitr(result.sig$gene_name,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db")
diff_gene_protein_entrezid$logFC <- result.sig[which(result.sig$gene_name%in%diff_gene_protein_entrezid$SYMBOL),]$logFC
# result.sig_entrezid <- bitr(row.names(result.sig),fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db")

gene_list_protein <- diff_gene_protein_entrezid$logFC
names(gene_list_protein) <- diff_gene_protein_entrezid$ENTREZID
gene_list_protein <- sort(gene_list_protein,decreasing = TRUE)
# gseKEGG_protein <- gseKEGG(geneList = gene_list_protein)
# gseKEGG_protein_dose <- gseKEGG(geneList = gene_list_protein,by = "DOSE")

enrichKEGG.protein <- enrichKEGG(diff_gene_protein_entrezid$ENTREZID,pvalueCutoff = 1,qvalueCutoff = 1)

cnetplot(enrichKEGG.protein,showCategory = 10)
dotplot(enrichKEGG.protein,color = "pvalue",showCategory = 10,font.size=8)	# 画气泡图
barplot(enrichKEGG.protein,color = "pvalue",showCategory = 10,font.size=8)
# browseKEGG(enrichKEGG.protein,'hsa04666')	# 显示通路图

entrezid2name <- function(gene_entrezid,sig_gene_entrezid2name){
  namelist <- gene_entrezid
  print(sig_gene_entrezid2name)
  for (id in 1:length(gene_entrezid)){
    namelist[id] <- gene_entrezid[id]
    entrez <- strsplit(namelist[id],split = "/")[[1]]
    
    for (i in 1:length(entrez)){
      namelist[id] <- gsub(entrez[i],replacement = sig_gene_entrezid2name[entrez[i],"SYMBOL"],namelist[id])
      
    }
  }
  return(namelist)
}

diff_gene_protein_kegg_term <- enrichKEGG.protein@result[,c("Description","geneID")]
# diff_gene_protein_kegg_term$geneName <- entrezid2name(diff_gene_protein_kegg_term$geneID, diff_gene_protein_entrezid)


write.csv(x = enrichGO.sig_gene@result,file = "sig_gene_go101.csv")
write.csv(x = enrichKEGG.protein@result,file = "sig_gene_kegg101.csv")
