# cohortPender_df: data frame of transcriptome (in TPM) & clinical data in the Liu cohort. 

# Fig. 3a
features=c('PSME1','PSMB10','TAP2','B2M','CIITA','HLA.E','HLA.G','HLA.F','HLA.DRB6','HLA.DQA2','HLA.DQB2')
combo_vec=unlist(features) 
df=subset(cohortPender_df,select=(c('patientId','best_response','tumorType','HLA.A','HLA.B','HLA.C',"HLA.DRA", "HLA.DRB1","HLA.DQA1","HLA.DQB1","HLA.DPA1","HLA.DPB1","HLA.DMA",'HLA.DMB',"HLA.DOA",'HLA.DOB',combo_vec)))
CRPR_df=subset(df,df$best_response == 'Partial response' | df$best_response == 'Physician assessed PR' | df$best_response == 'Partial Response')
SD_df=subset(df,df$best_response == 'Stable disease' | df$best_response == 'Physician assessed SD' | df$best_response == 'Physician Assessed SD')
PD_df=subset(df,df$best_response == 'Progressive disease' | df$best_response == 'Physician assessed PD' | df$best_response == 'Physician assessed PD ')
MR_df=subset(df,df$best_response == 'Mixed')
NE_df=subset(df,df$best_response == 'Not evaluable')
df=rbind(CRPR_df,SD_df,PD_df,MR_df,NE_df)
patientList=df$patientId
rownames(df)=patientList
tumorType=df$tumorType
df=df[4:ncol(df)]
df=log2(df+1)
colnames(df)=c('HLA-A','HLA-B','HLA-C',"HLA-DRA","HLA-DRB1","HLA-DQA1","HLA-DQB1","HLA-DPA1","HLA-DPB1","HLA-DMA",'HLA-DMB',"HLA-DOA",'HLA-DOB','PSME1','PSMB10','TAP2','B2M','CIITA','HLA-E','HLA-G','HLA-F','HLA-DRB6','HLA-DQA2','HLA-DQB2')
data=as.matrix(df)

annotated_clusters=cutree(out_Pan$tree_row, k=num_clusters)
clusters_df=data.frame('orig'=annotated_clusters)
clusters_df$new=NA
clusters_df$new[clusters_df$orig==4]=1
clusters_df$new[clusters_df$orig==1]=2
clusters_df$new[clusters_df$orig==2]=3
clusters_df$new[clusters_df$orig==3]=4
clusters_df$new=as.character(clusters_df$new)

annotation_row=data.frame(row.names=patientList, 
                          'Tumor Type' = tumorType,
                          'BOR' = c(rep("CRPR", nrow(CRPR_df)),rep("SD", nrow(SD_df)), rep("PD", nrow(PD_df)), rep('MR',nrow(MR_df)), rep('NE',nrow(NE_df))),
                          'Cluster'=clusters_df$new
)
colnames(annotation_row)=c('Tumor Type','BOR','Cluster')

# annotation_col <- data.frame(row.names = c('PSME1','PSMB10','TAP2',"HLA-A","HLA-B","HLA-C","HLA-E","HLA-G",'HLA-F','B2M','CIITA','HLA-DRB6','HLA-DQA2','HLA-DQB2',"HLA-DRA", "HLA-DRB1","HLA-DQA1","HLA-DQB1","HLA-DPA1","HLA-DPB1","HLA-DMA",'HLA-DMB',"HLA-DOA",'HLA-DOB'), 
#                             Genes = c(rep('Proteosome-related',3),rep("HLA Class I", 7), rep("HLA Class II", 14)))

colors=c(magma(10),viridis(10))
names(colors)=unique(annotation_row$`Tumor Type`)
annotation_colors = list('Tumor Type' = colors, 
                         'BOR' = c("CRPR" = "#0073C2FF","SD" = '#EFC000FF','MR' = 'snow2',"PD" = 'gray45','NE' = 'black'),
                         'Cluster' = c('1' = '#CD9B1D','2'='#EAD6A3','3'='#BFB1D6', '4'='#5E3C99'))

# Image: 590 * 570
my_palette <- colorRampPalette(c("deepskyblue4","white","firebrick3"))(200)
# heatmap: 590 * 650
# surv: 420 * 650
# DCB bar: 370 * 500
out_Pan=pheatmap::pheatmap(data, color = my_palette, cluster_cols=T,cluster_rows=T, angle_col = "315",
                           clustering_distance_rows = "euclidean", clustering_method = 'complete',
                           annotation_names_row = F,annotation_names_col = F, show_rownames = F, 
                           annotation_row=annotation_row,annotation_colors=annotation_colors,
                           legend=T, fontsize = 8,main='Pender Pan-cancer')

# Extended Data Fig. 3a
p_vec=c()
for (clusters in (2:8)) {
  num_clusters=clusters
  cohortPender_df$cluster=NA
  pre_cluster=sort(cutree(out_Pan$tree_row, k=num_clusters))
  for (i in (1:length(names(pre_cluster)))) {
    patientId=names(pre_cluster)[i]
    index=which(cohortPender_df$patientId == patientId)
    if (length(index) == 0) {
      message('Not Found: ', patientId)
    }
    cohortPender_df$cluster[index]=pre_cluster[i]
  }
  cohortPender_df$cluster=factor(cohortPender_df$cluster)
  coxfit=coxph(Surv(PFS_time,PFS) ~ cluster, data = cohortPender_df)
  adjusted_p=summary(coxfit)$sctest[3] * choose(clusters,2)
  p_vec=append(p_vec,adjusted_p)
}
chisq_vec=-log10(p_vec)
chisq_cohortPender_df=data.frame('clusters'=2:8,'pVal'=chisq_vec)
# 600 * 430
ggplot(chisq_cohortPender_df) +
  geom_step(aes(x=seq_along(chisq_cohortPender_df$clusters)-0.5, y=chisq_cohortPender_df$pVal)) +
  xlab("Clusters") +
  ylab("-log10(adjusted p value)") + theme_bw() +
  ggtitle('Pender Pan-cancer') +
  geom_hline(aes(yintercept=-log10(0.05)),linetype=2)+
  scale_x_discrete(name ="Clusters", 
                   limits=c("2","3","4",'5','6','7','8'='')) +
  annotate(geom = "text", y = c(-log10(0.05))+0.08, x = c(5.8), label = c('Adjusted p=0.05'), size = 3.5, hjust = 0)



num_clusters=4
cohortPender_df$cluster=NA
pre_cluster=sort(cutree(out_Pan$tree_row, k=num_clusters))
for (i in (1:length(names(pre_cluster)))) {
  patientId=names(pre_cluster)[i]
  index=which(cohortPender_df$patientId == patientId)
  if (length(index) == 0) {
    message('Not Found: ', patientId)
  }
  k=pre_cluster[i]
  if (k == 4) {
    k=1
  } else if (k==1) {
    k=2
  } else if (k==2) {
    k=3
  } else if (k==3) {
    k=4
  }
  cohortPender_df$cluster[index]=k
}
cohortPender_df$cluster=factor(cohortPender_df$cluster)
#write.csv(cohortPender_df,'/Users/katana/HLA-Evolution/CTLA4-science-dataset/HLA-clustering/cohort122_clusters.csv')

# Fig. 3b
cohortPender_df$Cluster=paste('C',cohortPender_df$cluster,sep='')
cohortPender_df$Cluster=factor(cohortPender_df$Cluster)
fit <- survfit(Surv(PFS_time,PFS) ~ Cluster,data = cohortPender_df)
# 530 * 560
ggsurvplot(fit, data = cohortPender_df,palette = c('#CD9B1D','#EAD6A3','#BFB1D6','#5E3C99'),
           conf.int = F, 
           pval = T,  
           risk.table = TRUE,
           xlab = '', ylab = 'Survival',
           legend.labs = levels(factor(cohortPender_df$Cluster)),
           ggtheme = theme_bw()+theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                      panel.background = element_blank()))
coxph(Surv(PFS_time,PFS) ~ Cluster, data = cohortPender_df) %>% gtsummary::tbl_regression(exp = TRUE)


# Fig. 3c
# Tumor types among clusters
data=subset(cohortPender_df,select=c('Cluster','tumorType'))
data$`Tumor Type`=data$tumorType
colors=length(unique(data$`Tumor Type`))
# 420 * 500
ggplot(data, aes(Cluster)) + 
  geom_bar(aes(fill=`Tumor Type`)) +
  scale_fill_manual(values=c(magma(colors/2),viridis(colors/2))) +
  xlab("") + ylab('Counts') + theme_bw()



# Fig. 3d
C1=subset(cohortPender_df,cohortPender_df$cluster == 1)
C2=subset(cohortPender_df,cohortPender_df$cluster == 2)
C3=subset(cohortPender_df,cohortPender_df$cluster == 3)
C4=subset(cohortPender_df,cohortPender_df$cluster == 4)
dfList=list(C1, C2, C3, C4)
HLA_supertypes=lapply(dfList, function(x) {
  supertypes = data.frame(mean(log2(x$HLA.A+1),na.rm=T), mean(log2(x$HLA.B+1),na.rm=T), mean(log2(x$HLA.C+1),na.rm=T), mean(log2(x$HLA.E+1),na.rm=T), mean(log2(x$HLA.F+1),na.rm=T), mean(log2(x$HLA.G+1),na.rm=T), mean(log2(x$B2M+1),na.rm=T), mean(log2(x$HLA.DRA+1),na.rm=T), mean(log2(x$HLA.DRB1+1),na.rm=T), mean(log2(x$HLA.DQA1+1),na.rm=T), mean(log2(x$HLA.DQB1+1),na.rm=T), mean(log2(x$HLA.DPA1+1),na.rm=T), mean(log2(x$HLA.DPB1+1),na.rm=T),mean(log2(x$HLA.DMA+1),na.rm=T), mean(log2(x$HLA.DMB+1),na.rm=T), mean(log2(x$HLA.DOA+1),na.rm=T), mean(log2(x$HLA.DOB+1),na.rm=T),            
                          mean(log2(x$PSME1+1),na.rm=T), mean(log2(x$PSME2+1),na.rm=T), mean(log2(x$PSMB8+1),na.rm=T), mean(log2(x$PSMB9+1),na.rm=T), mean(log2(x$PSMB10+1),na.rm=T))          
  
  colnames(supertypes)=c("HLA-A","HLA-B","HLA-C","HLA-E","HLA-F","HLA-G",'B2M',"HLA-DRA", "HLA-DRB1","HLA-DQA1","HLA-DQB1","HLA-DPA1","HLA-DPB1","HLA-DMA","HLA-DMB","HLA-DOA","HLA-DOB",
                         'PSME1','PSME2','PSMB8','PSMB9','PSMB10')
  return(supertypes)
} )

HLA_matrix <- data.frame(matrix(unlist(HLA_supertypes), nrow=length(HLA_supertypes), byrow=TRUE))
colnames(HLA_matrix)=c("HLA-A","HLA-B","HLA-C","HLA-E","HLA-F","HLA-G",'B2M',"HLA-DRA", "HLA-DRB1","HLA-DQA1","HLA-DQB1","HLA-DPA1","HLA-DPB1","HLA-DMA","HLA-DMB","HLA-DOA","HLA-DOB",
                       'PSME1','PSME2','PSMB8','PSMB9','PSMB10')
rownames(HLA_matrix)=c('C1','C2','C3','C4')
data=as.matrix(HLA_matrix)

annotation_col <- data.frame(row.names = c("HLA-A","HLA-B","HLA-C","HLA-E","HLA-F","HLA-G",'B2M',
                                           "HLA-DRA", "HLA-DRB1","HLA-DQA1","HLA-DQB1","HLA-DPA1","HLA-DPB1","HLA-DMA","HLA-DMB","HLA-DOA","HLA-DOB",
                                           'PSME1','PSME2','PSMB8','PSMB9','PSMB10'), 
                             Genes = c(rep("HLA Class I", 7), rep("HLA Class II", 10),rep("Proteasome", 5)))

annotation_colors = list(
  Genes = c('Proteasome' = 'dodgerblue4',"HLA Class I" = "forestgreen","HLA Class II" = 'thistle2')
)

my_palette <- colorRampPalette(c("deepskyblue4","white","goldenrod2"))(200)
# 600 * 500
pheatmap::pheatmap(data, color = my_palette, cluster_cols=T,cluster_rows=F, annotation_colors = annotation_colors,
                   annotation_col=annotation_col, annotation_legend = T, 
                   annotation_names_row = F,annotation_names_col = F, show_rownames = T,
                   border_color=NA, legend=T, angle_col = "315")



# Fig. 3e
# immune signatures
data=cohortPender_df
data$`CD8 T cells`=data$CD8T
data$`CD4 Tmem`=data$CD4ActTmem
data$`M1`=data$macrophage1
data$`M2`=data$macrophage2
data$`Tfh`=data$Th17
data$`NK cells`=data$NKcells
data$Tregs=data$Tregs

dat.m1 <- melt(data %>% filter(!is.na(cluster)),
               id.vars='cluster', 
               measure.vars=c('CD8 T cells','CD4 Tmem',
                              'M1','M2',
                              'Tfh','NK cells','Tregs'))
my_comparisons=list(c('C1','C2'),c('C3','C4'),c('C1','C3'))
dat.m1$Cluster=NA
dat.m1$Cluster[which(dat.m1$cluster == 1)]='C1'
dat.m1$Cluster[which(dat.m1$cluster == 2)]='C2'
dat.m1$Cluster[which(dat.m1$cluster == 3)]='C3'
dat.m1$Cluster[which(dat.m1$cluster == 4)]='C4'

# 620 * 360
ggplot(dat.m1, aes(x=Cluster, y=value, fill=Cluster)) +
  geom_violin() +
  stat_compare_means(comparisons = my_comparisons, size=3) + 
  geom_boxplot(width=0.1, fill="white") + 
  guides(fill=FALSE) +
  facet_grid(cols = vars(variable),scales = "free") + xlab('') + ylab('Immune Infiltration') +
  scale_fill_manual(values = c('#CD9B1D',"#EAD6A3", "#BFB1D6",'#5E3C99')) +
  theme_bw()


# Fig. 3f
# exhaustion signatures
data=cohortPender_df
data$LAG3=log2(data$LAG3+1)
data$PDCD1=log2(data$PDCD1+1)
data$CTLA4=log2(data$CTLA4+1)
data$HAVCR2=log2(data$HAVCR2+1)
data$TIGIT=log2(data$TIGIT+1)
data$CD274=log2(data$CD274+1)
data$PDCD1LG2=log2(data$PDCD1LG2+1)

dat.m1 <- melt(data %>% filter(!is.na(cluster)),
               id.vars='cluster', 
               measure.vars=c('LAG3','PDCD1','CTLA4',
                              'HAVCR2','TIGIT','CD274','PDCD1LG2'))

my_comparisons=list(c('C1','C2'),c('C3','C4'),c('C1','C3'))
dat.m1$Cluster=NA
dat.m1$Cluster[which(dat.m1$cluster == 1)]='C1'
dat.m1$Cluster[which(dat.m1$cluster == 2)]='C2'
dat.m1$Cluster[which(dat.m1$cluster == 3)]='C3'
dat.m1$Cluster[which(dat.m1$cluster == 4)]='C4'
# 620 * 360
ggplot(dat.m1, aes(x=Cluster, y=value, fill=Cluster)) +
  geom_violin() +
  stat_compare_means(comparisons = my_comparisons, size=3) + 
  geom_boxplot(width=0.1, fill="white") + xlab('') + ylab('Gene Expression') +
  guides(fill=FALSE) +
  facet_grid(cols = vars(variable),scales = "free") +
  scale_fill_manual(values = c('#CD9B1D',"#EAD6A3", "#BFB1D6",'#5E3C99')) +
  theme_bw()



# Extended Data Fig. 3b
facet_df=cohortPender_df
facet_df$`clinical benefit`=facet_df$clinical_benefit

C1=subset(facet_df,facet_df$cluster== 1)
C1_percent=paste(round(length(which(C1$`clinical benefit` == 'DCB')) / nrow(C1) * 100,0), '%',sep='')
C2=subset(facet_df,facet_df$cluster== 2)
C2_percent=paste(round(length(which(C2$`clinical benefit` == 'DCB')) / nrow(C2) * 100,0), '%',sep='')
C3=subset(facet_df,facet_df$cluster== 3)
C3_percent=paste(round(length(which(C3$`clinical benefit` == 'DCB')) / nrow(C3) * 100,0), '%',sep='')
C4=subset(facet_df,facet_df$cluster== 4)
C4_percent=paste(round(length(which(C4$`clinical benefit` == 'DCB')) / nrow(C4) * 100,0), '%',sep='')
data=subset(facet_df,select=c('Cluster','clinical benefit'))
# 370 * 500
ggplot(data, aes(Cluster,fill=`clinical benefit`)) + 
  geom_bar(position = "fill") +
  scale_fill_manual(values=c("gold2", "gray24")) + ylim(c(0,1.1)) +
  xlab("") + ylab('Fraction') + theme_bw() + theme(text = element_text(size=8)) +
  annotate(geom = "text", y = c(1,1,1,1) - 0.1, x = c(0.95,2,3,4)-0.2, label = c(C1_percent,C2_percent,C3_percent,C4_percent), size = 3.5, hjust = 0)


# Fisher's
df=facet_df
C1_DCB=subset(df, df$`clinical benefit` == 'DCB' & df$cluster == 1)
C1_NCB=subset(df, df$`clinical benefit` == 'NCB' & df$cluster == 1)
C2_DCB=subset(df, df$`clinical benefit` == 'DCB' & df$cluster == 2)
C2_NCB=subset(df, df$`clinical benefit` == 'NCB' & df$cluster == 2)
C3_DCB=subset(df, df$`clinical benefit` == 'DCB' & df$cluster == 3)
C3_NCB=subset(df, df$`clinical benefit` == 'NCB' & df$cluster == 3)
C4_DCB=subset(df, df$`clinical benefit` == 'DCB' & df$cluster == 4)
C4_NCB=subset(df, df$`clinical benefit` == 'NCB' & df$cluster == 4)

C1=data.frame('DCB'= nrow(C1_DCB),'NCB'= nrow(C1_NCB))
C2=data.frame('DCB'= nrow(C2_DCB),'NCB'= nrow(C2_NCB))
C3=data.frame('DCB'= nrow(C3_DCB),'NCB'= nrow(C3_NCB))
C4=data.frame('DCB'= nrow(C4_DCB),'NCB'= nrow(C4_NCB))
df_combined=rbind(C2,C3)
fisher.test(df_combined,alternative = "two.sided")




# Extended Data Fig. 3c
# proinflammatory
data=cohortPender_df
data$IFNG=log2(data$IFNG+1)
data$IDO1=log2(data$IDO1+1)
data$STAT1=log2(data$STAT1+1)
data$CXCL9=log2(data$CXCL9+1)
data$CXCL10=log2(data$CXCL10+1)
dat.m1 <- melt(data %>% filter(!is.na(cluster)),
               id.vars='cluster', 
               measure.vars=c('IFNG','IDO1','STAT1','CXCL9',
                              'CXCL10'))

my_comparisons=list(c('C1','C2'),c('C3','C4'),c('C1','C3'))
dat.m1$Cluster=NA
dat.m1$Cluster[which(dat.m1$cluster == 1)]='C1'
dat.m1$Cluster[which(dat.m1$cluster == 2)]='C2'
dat.m1$Cluster[which(dat.m1$cluster == 3)]='C3'
dat.m1$Cluster[which(dat.m1$cluster == 4)]='C4'
# 620 * 360
ggplot(dat.m1, aes(x=Cluster, y=value, fill=Cluster)) +
  geom_violin() +
  stat_compare_means(comparisons = my_comparisons, size=3) + 
  geom_boxplot(width=0.1, fill="white") + xlab('') + ylab('Gene Expression') +
  guides(fill=FALSE) +
  facet_grid(cols = vars(variable),scales = "free") +
  scale_fill_manual(values = c('#EAD6A3',"#CD9B1D", "#BFB1D6",'#5E3C99')) +
  theme_bw()




# 650 * 400
dat.m1 <- melt(cohortPender_df %>% filter(!is.na(cluster)),
               id.vars='cluster', 
               measure.vars=c('SSGSEA_PROLIFERATION','SSGSEA_ANGIOGENESIS'))

dat.m1$metrics=NA
dat.m1$metrics[which(dat.m1$variable == 'SSGSEA_PROLIFERATION')]='Proliferation Score'
dat.m1$metrics[which(dat.m1$variable == 'SSGSEA_ANGIOGENESIS')]='Angiogenesis Score'

my_comparisons=list(c('C1','C2'),c('C3','C4'),c('C1','C3'))
dat.m1$Cluster=NA
dat.m1$Cluster[which(dat.m1$cluster == 1)]='C1'
dat.m1$Cluster[which(dat.m1$cluster == 2)]='C2'
dat.m1$Cluster[which(dat.m1$cluster == 3)]='C3'
dat.m1$Cluster[which(dat.m1$cluster == 4)]='C4'
dat.m1$metrics=factor(dat.m1$metrics,c('Proliferation Score','Angiogenesis Score'))

# Extended Data Fig. 3d
# 400 * 500
ggplot(cohortPan_df, aes(x=Cluster, y=SSGSEA_PROLIFERATION, fill=Cluster)) +
  geom_violin() +
  stat_compare_means(comparisons = my_comparisons, size=3) + 
  geom_boxplot(width=0.1, fill="white") + xlab('') + ylab('Proliferation Enrichment Score') +
  guides(fill=FALSE) +
  scale_fill_manual(values = c('#CD9B1D',"#EAD6A3", "#BFB1D6",'#5E3C99')) +
  theme_bw()
# Extended Data Fig. 3e
ggplot(cohortPan_df, aes(x=Cluster, y=SSGSEA_ANGIOGENESIS, fill=Cluster)) +
  geom_violin() +
  stat_compare_means(comparisons = my_comparisons, size=3) + 
  geom_boxplot(width=0.1, fill="white") + xlab('') + ylab('Angiogenesis Enrichment Score') +
  guides(fill=FALSE) +
  scale_fill_manual(values = c('#CD9B1D',"#EAD6A3", "#BFB1D6",'#5E3C99')) +
  theme_bw()







