subset.T.raw=subset(filtered_seurat,subset=(cell.types == 'T.CD4' | cell.types == 'T.CD8'))
subset.T.umap <- SCTransform(subset.T.raw, verbose = TRUE)
DefaultAssay(subset.T.umap)='SCT'
features.include=setdiff(VariableFeatures(subset.T.umap),genes.exclude)
subset.T.umap <- ScaleData(subset.T.umap,features=features.include)
subset.T.umap <- RunPCA(subset.T.umap,features=features.include)
subset.T.umap <- FindNeighbors(subset.T.umap, dims = 1:30)
subset.T.umap <- FindClusters(subset.T.umap, resolution = c(0.7,0.6,0.5))
subset.T.umap <- RunUMAP(subset.T.umap, dims = 1:30)

Idents(subset.T.umap)='SCT_snn_res.0.7'
Idents(subset.T.umap)='cell.types' # from original author
Idents(subset.T.umap)='samples'
Idents(subset.T.umap)='Cluster'
Idents(subset.T.umap)='Cell.Type'

metadata=subset.T.umap@meta.data
metadata$Cluster=NA
for (patientId in patientId_filtered) {
  index=which(metadata$samples == patientId)
  metadata$Cluster[index]=pseudobulk_cluster_vec[patientId]
}
metadata$Cluster=factor(metadata$Cluster,c('C1','C2','C3','C4'))
subset.T.umap@meta.data=metadata

DimPlot(subset.T.umap, reduction = "umap",label = F, pt.size=0.01, raster = FALSE)
ggsave("UMAP1.png",height = 3,width =6,dpi = 300) 

DefaultAssay(subset.T.umap) <- "SCT"
FeaturePlot(subset.T.umap, cols = c("gray88", "forestgreen"), features = c("CD8A",'GZMK','ITGAE','IFNG','TBX21',
                                                                           'CD4','FOXP3','CXCL13','STAT4','IL21',
                                                                           "TCF7","SELL",'CCR7','CD69','MKI67',
                                                                           'PDCD1','CTLA4','HAVCR2','LAG3',"ENTPD1"), ncol = 5, order = T, pt.size = 0.2) & NoLegend() & NoAxes()
ggsave("UMAP2.png",height = 7,width =10,dpi = 300)



# compute DEGs across UMAP clusters
DefaultAssay(subset.T.umap)='SCT'
markers_subset_T <- FindAllMarkers(subset.T.umap, only.pos = F, min.pct = 0.1, logfc.threshold = 0)
markers_subset_T <- markers_subset_T %>%
  dplyr::arrange(cluster, avg_log2FC*(-1))
markers_subset_T <- markers_subset_T[ , c(6, 7, 2:4, 1, 5)]
markers_subset_T_filtered=subset(markers_subset_T,markers_subset_T$p_val_adj<=0.05)


subset.T.umap <- RenameIdents(object = subset.T.umap, 
                              '0' = 'CD8-Trm',
                              '1' = 'Naïve-like',
                              '2' = 'Undefined', 
                              '3' = 'CD8-Effector', 
                              '4' = 'CD8-Proliferating', 
                              '5' = 'CD8-T Cells', 
                              '6' = 'CD4-Treg', 
                              '7' = 'CD4-Tfh',
                              '8' = 'Undefined'
)
subset.T.umap$Cell.Type=Idents(subset.T.umap)
metadata=subset.T.umap@meta.data
metadata$Cell.Type=factor(metadata$Cell.Type, c('CD8-Trm','CD8-Effector','CD8-Proliferating','CD8-T Cells',
                                                'CD4-Treg','CD4-Tfh','Naïve-like','Undefined'))
subset.T.umap@meta.data=metadata

library(viridis)
ggplot(data=metadata, aes(x=Cluster,fill=Cell.Type)) +
  geom_bar(aes(fill=Cell.Type),position="fill") + xlab('') + ylab('Fraction') +
  scale_fill_manual(values=c(rocket(8)[4:7],mako(8)[3:5],cividis(5)[3])) +
  theme_bw() + labs(fill = "Cell Type")
ggsave("bar.png",height = 5,width =4.5,dpi = 300)

ggplot(data=metadata, aes(x=Cell.Type,fill=Cluster)) +
  geom_bar(aes(fill=Cluster),position="fill") + xlab('') + ylab('Fraction') +
  scale_fill_manual(values=c(mako(8)[c(5,4,2,1)])) +
  theme_bw()
ggsave("bar.png",height = 2,width =5,dpi = 300)



# fishers
cellType_vec=rev(c('CD8-Trm','CD8-Effector','CD8-Proliferating','CD8-T Cells',
                   'CD4-Treg','CD4-Tfh','Naïve-like','Undefined'))
cluster_vec=c('C1','C2','C3','C4')

fisher_stats=data.frame()
fisher_pVal=data.frame()
for (cellType in cellType_vec) {
  pVal_vec=c()
  odds_vec=c()
  for (cluster_num in cluster_vec) {
    n1=sum(metadata$Cell.Type == cellType & metadata$Cluster == cluster_num)
    n2=sum(metadata$Cell.Type != cellType & metadata$Cluster == cluster_num)
    n3=sum(metadata$Cell.Type == cellType & metadata$Cluster != cluster_num)
    n4=sum(metadata$Cell.Type != cellType & metadata$Cluster != cluster_num)
    df_combined=data.frame('cluster1'=c(n1,n2),'cluster2'=c(n3,n4))
    fisher_results=fisher.test(df_combined,alternative = "two.sided")
    pVal_vec=append(pVal_vec,fisher_results$p.value)
    odds_vec=append(odds_vec,fisher_results$estimate)
  }
  temp_stats=data.frame('Cell.Type'=cellType,'Cluster'=cluster_vec,'Odds'=odds_vec)
  temp_pVal=data.frame('Cell.Type'=cellType,'Cluster'=cluster_vec,'pval'=pVal_vec)
  
  fisher_stats=rbind(fisher_stats,temp_stats)
  fisher_pVal=rbind(fisher_pVal,temp_pVal)
}

plot.data <- fisher_stats
plot.data$log2OR=log2(plot.data$Odds+0.1)
plot.pVal <- fisher_pVal
qVal=p.adjust(plot.pVal$pval,method="BH") 
plot.pVal$p.adj=qVal
plot.data$stars <- cut(plot.pVal$p.adj, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))  # Create column of significance labels

plot.data$Cell.Type=factor(plot.data$Cell.Type,cellType_vec)
p <- ggplot(aes(x=Cell.Type, y=Cluster, fill=log2OR), data=plot.data)
p + geom_tile() + scale_fill_gradient2(low='dodgerblue4',mid='gray97',high='gold3') + 
  geom_text(aes(label=stars), color="black", size=5) + 
  labs(y=NULL, x=NULL, fill="log2(OR)") + coord_flip() + theme_bw() +
  theme(axis.text.x=element_text(angle = -45, hjust = 0),legend.title=element_text(size=8)) +
  guides(fill = guide_colourbar(barwidth = 1,
                                barheight = 10))
ggsave("OR.png",height = 3.5,width =4,dpi = 300)





# dot plot
EXHAUSTION.GENES=read.table('exhaustion.txt',header=T)$EXHAUSTION
IL7.GENES=read.table('IL7-Induced.txt',header=T)$IL7.INDUCED
TRM.GENES=read.table('TRM.txt',header=T)$TRM
STEMNESS.GENES=read.table('stemness.txt')$V1
CYTOTOXICITY.GENES=read.table('cytotoxicity.txt')$V1
CYTOKINES.GENES=read.table('cytokine.txt')$V1

immune.functions.pathways=list('Immune Exhaustion'=EXHAUSTION.GENES,'Trm Signature'=TRM.GENES,
                               'Immune Stemness'=STEMNESS.GENES,'Immune Cytotoxicity'=CYTOTOXICITY.GENES, 'Cytokines Signature'=CYTOKINES.GENES)
mSigDB.hallmark=gmtPathways('h.all.v7.5.1.symbols.gmt.txt')


clusters=unique(markers_df$cluster)
combined_df=data.frame()
for (cluster.group in clusters) {
  cluster.genes<- markers_df %>%
    dplyr::filter(cluster == cluster.group) %>%
    arrange(desc(avg_log2FC)) %>% 
    dplyr::select(gene, avg_log2FC)
  ranks=cluster.genes$avg_log2FC
  names(ranks)=cluster.genes$gene
  fgseaRes <- fgsea(immune.functions.pathways, ranks, maxSize=1000, nPermSimple = 10000)
  message(cluster.group,': ',sum(fgseaRes$padj<0.25,na.rm=T)) 
  
  if (sum(fgseaRes$padj<0.25,na.rm=T)>0) {
    fgseaRes.sig=fgseaRes[which(fgseaRes$padj<0.25),]
    fgseaRes.sig=fgseaRes.sig %>% arrange(-1*fgseaRes.sig$NES)
    print(fgseaRes.sig[,c('pathway','NES','padj')])
    
    combined_df=rbind(combined_df,data.frame('Cell Type'=cluster.group,fgseaRes.sig[,c('pathway','NES','padj')]))
  }
}
combined_df$log10.PVal=-log10(combined_df$padj)
ggplot(data=combined_df, aes(Cell.Type, pathway)) +
  geom_point(aes(size = log10.PVal, color=NES)) + 
  geom_point(data=combined_df %>% filter(padj<0.05), aes(Cell.Type, pathway, size = log10.PVal*3),pch=21)+
  scale_colour_gradient2(low="deepskyblue4", mid="gray80", high="firebrick4") +
  theme_classic() + scale_size_continuous(range = c(2.5, 10)) +
  theme(text = element_text(size=17),axis.text.x=element_text(angle = -30, hjust = 0)) + 
  xlab('') + ylab('')
ggsave("celltype_NES.png",height = 3,width =8,dpi = 300)




