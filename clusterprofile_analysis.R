
DIR=getwd()
setwd("E:\\Rcode\\data")

# clusterprofile analysis

library(clusterProfiler)
library(org.Hs.eg.db)

diff_gene_protein <- read.csv("diff_expr_protein.csv",row.names = 1)
result.sig <- read.csv("sig_gene_v3.csv",row.names = 1)

enrichGO.protein.CC <- enrichGO(gene = diff_gene_protein$gene_name,OrgDb="org.Hs.eg.db",keyType = "SYMBOL",pvalueCutoff = 0.01,qvalueCutoff = 0.05,ont = "CC",pAdjustMethod = "BH")
enrichGO.protein <- enrichGO(gene = diff_gene_protein$gene_name,OrgDb="org.Hs.eg.db",keyType = "SYMBOL",pvalueCutoff = 0.01,qvalueCutoff = 0.05,ont = "all",pAdjustMethod = "BH")
enrichGO.protein.BP <- enrichGO(gene = diff_gene_protein$gene_name,OrgDb="org.Hs.eg.db",keyType = "SYMBOL",pvalueCutoff = 0.01,qvalueCutoff = 0.05,ont = "BP",pAdjustMethod = "BH")
enrichGO.protein.MF <- enrichGO(gene = diff_gene_protein$gene_name,OrgDb="org.Hs.eg.db",keyType = "SYMBOL",pvalueCutoff = 0.01,qvalueCutoff = 0.05,ont = "MF",pAdjustMethod = "BH")

dotplot(enrichGO.protein,showCategory=10,title="Enrichment GO Top10") #泡泡图
barplot(enrichGO.protein, showCategory=10,title="EnrichmentGO Top10")  #柱状图
# plotGOgraph(enrichGO.protein) 	#GO图，看不清楚可以尝试左上角另存为pdf
emapplot(enrichGO.protein)

cnetplot(enrichGO.protein,foldChange = diff_gene_protein$logFC,circular = TRUE,colorEdge = TRUE)


diff_gene_protein_entrezid <- bitr(diff_gene_protein$gene_name,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db")
diff_gene_protein_entrezid$logFC <- diff_gene_protein[which(diff_gene_protein$gene_name%in%diff_gene_protein_entrezid$SYMBOL),]$logFC
# result.sig_entrezid <- bitr(row.names(result.sig),fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db")


gene_list_protein <- diff_gene_protein_entrezid$logFC
names(gene_list_protein) <- diff_gene_protein_entrezid$ENTREZID
gene_list_protein <- sort(gene_list_protein,decreasing = TRUE)
# gseKEGG_protein <- gseKEGG(geneList = gene_list_protein)
# gseKEGG_protein_dose <- gseKEGG(geneList = gene_list_protein,by = "DOSE")

enrichKEGG.protein <- enrichKEGG(diff_gene_protein_entrezid$ENTREZID,pvalueCutoff = 0.1,qvalueCutoff = 1)

cnetplot(enrichKEGG.protein,foldChange = gene_list_protein)
dotplot(enrichKEGG.protein,font.size=8)	# 画气泡图
browseKEGG(enrichKEGG.protein,'mmu01100')	# 显示通路图

# sig_gene
sig_gene_entrezid <- diff_gene_protein_entrezid[which(diff_gene_protein_entrezid$SYMBOL%in%names(sig_gene)),]
row.names(sig_gene_entrezid) <- sig_gene_entrezid$SYMBOL
sig_gene_entrezid <- sig_gene_entrezid[names(sig_gene),]
row.names(sig_gene_entrezid) <- sig_gene_entrezid$ENTREZID
enrichKEGG.sig_gene <- enrichKEGG(sig_gene_entrezid$ENTREZID,qvalueCutoff = 0.2,pvalueCutoff = 0.2)
sig_gene_kegg_term <- enrichKEGG.sig_gene@result[,c("Description","geneID")]
sig_gene_kegg_term$geneName <- entrezid2name()

entrezid2name <- function(gene_entrezid=sig_gene_kegg_term$geneID,sig_gene_entrezid2name = sig_gene_entrezid){
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

write.csv(x = enrichKEGG.protein,file = "clusterprofile_result\\kegg_mrna.csv")
write.csv(x = enrichKEGG.sig_gene,file = "clusterprofile_result\\sig_gene_kegg_summary.csv")
write.csv(x = enrichKEGG.sig_gene@result,file = "clusterprofile_result\\sig_gene_kegg.csv")


enrichKEGG_result.sig <- enrichKEGG(gene =  result.sig_entrezid$ENTREZID,keyType = "kegg",use_internal_data = TRUE)


path_result.sig <- pathview(gene.data = diff_gene_protein_entrezid$ENTREZID,pathway.id = )
pathview_protein <-pathview(gene.data = diff_gene_protein_entrezid$ENTREZID,pathway.id = gseKEGG_protein@result$ID)



barplot(enrichGO.protein,showCategory = 30)
emapplot(enrichGO.protein,showCategory = 30)
cnetplot(enrichGO.protein,showCategory = 30)

barplot(gseKEGG_protein,showCategory = 30)