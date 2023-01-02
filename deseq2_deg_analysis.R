library( "DESeq2" )
cohort_df$group=as.vector(cohort_df$Cluster)
cohort_df$group[which(cohort_df$group == 'C1' | cohort_df$group == 'C2')]='C1_2'
cohort_df$group[which(cohort_df$group == 'C3' | cohort_df$group == 'C4')]='C3_4'
unique(cohort_df$group)
metadata=data.frame('srrid'=cohort_df$SRRid,'group'=cohort_df$group)

dds <- DESeqDataSetFromMatrix(countData=expression_counts, 
                                 colData=metadata, 
                                 design=~group, tidy = TRUE)
dds <- DESeq(dds)
res <- results(dds,contrast=c("group","C1_2", "C3_4"))
res_df=data.frame(res)
res_df$gene=rownames(res_df)
res_df=res_df %>% arrange(res_df$log2FoldChange*(-1))



library(topGO)
library(biomaRt)
library(org.Hs.eg.db)
library(fgsea)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(tibble)
library(reactome.db)
mSigDB.hallmark=gmtPathways('h.all.v7.5.1.symbols.gmt.txt')

markers_df=res_df

cluster.genes<- markers_df %>%
  arrange(-1*log2FoldChange) %>% 
  dplyr::select(gene, log2FoldChange)
ranks=cluster.genes$log2FoldChange
names(ranks)=cluster.genes$gene
fgseaRes <- fgsea(mSigDB.hallmark, ranks, maxSize=1000, nPermSimple = 1000)
message(sum(fgseaRes$padj<=0.05,na.rm=T))
fgseaRes.sig=fgseaRes[which(fgseaRes$padj<=0.05),]
fgseaRes.sig=fgseaRes.sig %>% arrange(-1*fgseaRes.sig$NES)
print(fgseaRes.sig[,c('pathway','NES','padj')])


fgseaRes_plot=fgseaRes.sig
fgseaRes_plot$FDR=fgseaRes_plot$padj
ggplot(fgseaRes_plot, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill= FDR)) +
  coord_flip() +
  scale_fill_distiller() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="") + 
  theme_minimal()
ggsave("gsea.png",height = 5,width = 7,dpi = 300)
