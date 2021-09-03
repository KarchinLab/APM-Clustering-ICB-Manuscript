# cohortVanAllen_df: data frame of transcriptome (in TPM) & clinical data in the Van Allen cohort. 

# housekeeping gene score: for normalization
cohortVanAllen_df$ACTB_log=log2(cohortVanAllen_df$ACTB+1)
cohortVanAllen_df$GAPDH_log=log2(cohortVanAllen_df$GAPDH+1)
cohortVanAllen_df$SDHA_log=log2(cohortVanAllen_df$SDHA+1)
cohortVanAllen_df$HPRT1_log=log2(cohortVanAllen_df$HPRT1+1)
cohortVanAllen_df$HBS1L_log=log2(cohortVanAllen_df$HBS1L+1)
cohortVanAllen_df$AHSP_log=log2(cohortVanAllen_df$AHSP+1)
cohortVanAllen_df$housekeeping=(cohortVanAllen_df$ACTB_log+cohortVanAllen_df$GAPDH_log+cohortVanAllen_df$SDHA_log+cohortVanAllen_df$HPRT1_log+cohortVanAllen_df$HBS1L_log+cohortVanAllen_df$AHSP_log)/6
housekeeping_mean=mean(cohortVanAllen_df$housekeeping,na.rm=T)


# Fig 1b
features=c('PSME1','PSMB10','TAP2','B2M','CIITA','HLA_E','HLA_G','HLA_F','HLA_DRB6','HLA_DQA2','HLA_DQB2')
combo_vec=unlist(features)
df=subset(cohortVanAllen_df,select=(c('patientId','BOR',"progression_free","progression",'HLA_A','HLA_B','HLA_C',"HLA_DRA", "HLA_DRB1","HLA_DQA1","HLA_DQB1","HLA_DPA1","HLA_DPB1","HLA_DMA",'HLA_DMB',"HLA_DOA",'HLA_DOB',combo_vec)))

CRPR_df=subset(df,df$BOR == 'CR' | df$BOR == 'PR')
SD_df=subset(df,df$BOR == 'SD')
PD_df=subset(df,df$BOR == 'PD')
X_df=subset(df,df$BOR == 'X')

df=rbind(CRPR_df,SD_df,PD_df, X_df)
patientList=df$patientId
rownames(df)=patientList
df[5:ncol(df)]=log2(df[5:ncol(df)]+1)
df=df[5:ncol(df)]
colnames(df)=c('HLA-A','HLA-B','HLA-C',"HLA-DRA", "HLA-DRB1","HLA-DQA1","HLA-DQB1","HLA-DPA1","HLA-DPB1","HLA-DMA",'HLA-DMB',"HLA-DOA",'HLA-DOB','PSME1','PSMB10','TAP2','B2M','CIITA','HLA-E','HLA-G','HLA-F','HLA-DRB6','HLA-DQA2','HLA-DQB2')
data=as.matrix(df)


annotated_clusters=cutree(out_38$tree_row, k=num_clusters)
clusters_df=data.frame('orig'=annotated_clusters)
clusters_df$new=clusters_df$orig
clusters_df$new=as.character(clusters_df$new)
annotation_row <- data.frame(row.names = patientList, 
                             'BOR' = c(rep("CRPR", nrow(CRPR_df)),rep("SD", nrow(SD_df)), rep("PD", nrow(PD_df)), rep('NE', nrow(X_df))),
                             'Cluster'=clusters_df$new)
annotation_colors = list(
  'Cluster' = c('1' = '#CD9B1D','2'='#EAD6A3','3'='#BFB1D6', '4'='#5E3C99'),
  'BOR' = c("CRPR" = "#0073C2FF","SD" = '#EFC000FF',"PD" = 'gray35','NE' = 'black')
)

# Image: 590 * 570
my_palette <- colorRampPalette(c("deepskyblue4","white","firebrick3"))(200)

# heatmap: 590 * 650
# surv: 420 * 650
# DCB bar: 370 * 500
out_38=pheatmap::pheatmap(data, color = my_palette, cluster_cols=T,cluster_rows=T, angle_col = "315",
                          annotation_names_row = F,annotation_names_col = F, show_rownames = F, 
                          annotation_colors=annotation_colors,annotation_row=annotation_row,
                          legend=T, fontsize = 8,main='Van Allen Melanoma anti-CTLA4')


num_clusters=4
cohortVanAllen_df$cluster=NA
pre_cluster=sort(cutree(out_38$tree_row, k=num_clusters))
for (i in (1:length(names(pre_cluster)))) {
  patientId=names(pre_cluster)[i]
  index=which(cohortVanAllen_df$patientId == patientId)
  if (length(index) == 0) {
    message('Not Found: ', patientId)
  }
  k=pre_cluster[i]
  # attention
  if (k==1) {
    k=1
  } else if (k==2) {
    k=2
  } else if (k==3) {
    k=3
  } else if (k==4) {
    k=4
  }
  cohortVanAllen_df$cluster[index]=k
}
cohortVanAllen_df$cluster=factor(cohortVanAllen_df$cluster)

cohortVanAllen_df$Cluster=paste('C',cohortVanAllen_df$cluster,sep='')
cohortVanAllen_df$Cluster=factor(cohortVanAllen_df$Cluster)

# Fig. 1d
fit <- survfit(Surv(progression_free,progression) ~ Cluster,data = cohortVanAllen_df)
# 430 * 550
ggsurvplot(fit, data = cohortVanAllen_df, palette = c('#CD9B1D',"#EAD6A3", "#BFB1D6",'#5E3C99'),
           conf.int = F, 
           pval = T,  
           risk.table = TRUE,
           xlab = '', ylab = 'Survival', 
           legend.labs = levels(factor(cohortVanAllen_df$Cluster)),
           ggtheme = theme_bw()+theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                      panel.background = element_blank()))
coxph(Surv(progression_free,progression) ~ Cluster, data = cohortVanAllen_df) %>% gtsummary::tbl_regression(exp = TRUE)


# Fig. 1f
# 420 * 500
facet_df=cohortVanAllen_df
facet_df$`clinical benefit`=NA
facet_df$`clinical benefit`[which(facet_df$BOR == 'PD')]='NCB'
facet_df$`clinical benefit`[which(facet_df$BOR != 'PD')]='DCB'
facet_df$`clinical benefit`[which(facet_df$BOR == 'X' & facet_df$OS < 180)]='NCB'

C1=subset(facet_df,facet_df$cluster== 1)
C1_percent=paste(round(length(which(C1$`clinical benefit` == 'DCB')) / nrow(C1) * 100,0), '%',sep='')
C2=subset(facet_df,facet_df$cluster== 2)
C2_percent=paste(round(length(which(C2$`clinical benefit` == 'DCB')) / nrow(C2) * 100,0), '%',sep='')
C3=subset(facet_df,facet_df$cluster== 3)
C3_percent=paste(round(length(which(C3$`clinical benefit` == 'DCB')) / nrow(C3) * 100,0), '%',sep='')
C4=subset(facet_df,facet_df$cluster== 4)
C4_percent=paste(round(length(which(C4$`clinical benefit` == 'DCB')) / nrow(C4) * 100,0), '%',sep='')

# 420 * 500
ggplot(facet_df, aes(Cluster,fill=`clinical benefit`)) + 
  geom_bar(position = "fill") +
  scale_fill_manual(values=c("gold2", "gray24")) + ylim(c(0,1.1)) +
  xlab("") + ylab('Fraction') + theme_bw() + theme(text = element_text(size=8)) +
  annotate(geom = "text", y = c(1,1,1) - 0.1, x = c(1,2,3)-0.2, label = c(C1_percent,C2_percent,C3_percent), size = 3, hjust = 0)

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
df_combined=rbind(C1,C4)
fisher.test(df_combined,alternative = "two.sided")


# Extended Fig. 1b
surv_df=cohortVanAllen_df
surv_df$B2M=surv_df$B2M > unname(quantile(surv_df$B2M,na.rm=TRUE)["50%"])
surv_df$CIITA=surv_df$CIITA > unname(quantile(surv_df$CIITA,na.rm=TRUE)["50%"])
surv_df$`HLA-A`=surv_df$HLA_A > unname(quantile(surv_df$HLA_A,na.rm=TRUE)["50%"])
surv_df$`HLA-B`=surv_df$HLA_B > unname(quantile(surv_df$HLA_B,na.rm=TRUE)["50%"])
surv_df$`HLA-C`=surv_df$HLA_C > unname(quantile(surv_df$HLA_C,na.rm=TRUE)["50%"])
surv_df$`HLA-E`=surv_df$HLA_E > unname(quantile(surv_df$HLA_E,na.rm=TRUE)["50%"])
surv_df$`HLA-G`=surv_df$HLA_G > unname(quantile(surv_df$HLA_G,na.rm=TRUE)["50%"])
surv_df$`HLA-F`=surv_df$HLA_F > unname(quantile(surv_df$HLA_F,na.rm=TRUE)["50%"])
surv_df$`HLA-DRA`=surv_df$HLA_DRA > unname(quantile(surv_df$HLA_DRA,na.rm=TRUE)["50%"])
surv_df$`HLA-DRB1`=surv_df$HLA_DRB1 > unname(quantile(surv_df$HLA_DRB1,na.rm=TRUE)["50%"])
surv_df$`HLA-DQA1`=surv_df$HLA_DQA1 > unname(quantile(surv_df$HLA_DQA1,na.rm=TRUE)["50%"])
surv_df$`HLA-DQB1`=surv_df$HLA_DQB1 > unname(quantile(surv_df$HLA_DQB1,na.rm=TRUE)["50%"])
surv_df$`HLA-DPA1`=surv_df$HLA_DPA1 > unname(quantile(surv_df$HLA_DPA1,na.rm=TRUE)["50%"])
surv_df$`HLA-DPB1`=surv_df$HLA_DPB1 > unname(quantile(surv_df$HLA_DPB1,na.rm=TRUE)["50%"])
surv_df$`HLA-DMA`=surv_df$HLA_DMA > unname(quantile(surv_df$HLA_DMA,na.rm=TRUE)["50%"])
surv_df$`HLA-DMB`=surv_df$HLA_DMB > unname(quantile(surv_df$HLA_DMB,na.rm=TRUE)["50%"])
surv_df$`HLA-DOA`=surv_df$HLA_DOA > unname(quantile(surv_df$HLA_DOA,na.rm=TRUE)["50%"])
surv_df$`HLA-DOB`=surv_df$HLA_DOB > unname(quantile(surv_df$HLA_DOB,na.rm=TRUE)["50%"])
surv_df$PSME1=surv_df$PSME1 > unname(quantile(surv_df$PSME1,na.rm=TRUE)["50%"])
surv_df$PSME2=surv_df$PSME2 > unname(quantile(surv_df$PSME2,na.rm=TRUE)["50%"])
surv_df$ERAP1=surv_df$ERAP1 > unname(quantile(surv_df$ERAP1,na.rm=TRUE)["50%"])
surv_df$TAPBP=surv_df$TAPBP > unname(quantile(surv_df$TAPBP,na.rm=TRUE)["50%"])
surv_df$NLRC5=surv_df$NLRC5 > unname(quantile(surv_df$NLRC5,na.rm=TRUE)["50%"])
surv_df$PSMB8=surv_df$PSMB8 > unname(quantile(surv_df$PSMB8,na.rm=TRUE)["50%"])
surv_df$PSMB9=surv_df$PSMB9 > unname(quantile(surv_df$PSMB9,na.rm=TRUE)["50%"])
surv_df$PSMB10=surv_df$PSMB10 > unname(quantile(surv_df$PSMB10,na.rm=TRUE)["50%"])
surv_df$TAP1=surv_df$TAP1 > unname(quantile(surv_df$TAP1,na.rm=TRUE)["50%"])
surv_df$TAP2=surv_df$TAP2 > unname(quantile(surv_df$TAP2,na.rm=TRUE)["50%"])
surv_df$`HLA-DRB5`=surv_df$HLA_DRB5 > unname(quantile(surv_df$HLA_DRB5,na.rm=TRUE)["50%"])
surv_df$`HLA-DRB6`=surv_df$HLA_DRB6 > unname(quantile(surv_df$HLA_DRB6,na.rm=TRUE)["50%"])
surv_df$`HLA-DQA2`=surv_df$HLA_DQA2 > unname(quantile(surv_df$HLA_DQA2,na.rm=TRUE)["50%"])
surv_df$`HLA-DQB2`=surv_df$HLA_DQB2 > unname(quantile(surv_df$HLA_DQB2,na.rm=TRUE)["50%"])
genes=c('B2M',"HLA-A","HLA-B","HLA-C","HLA-E","HLA-G",'HLA-F','CIITA',"HLA-DRA", "HLA-DRB1",
        "HLA-DQA1","HLA-DQB1","HLA-DPA1","HLA-DPB1","HLA-DMA",'HLA-DMB',"HLA-DOA",'HLA-DOB',
        'PSME1','PSME2','ERAP1','TAPBP','NLRC5','PSMB8','PSMB9','PSMB10','TAP1','TAP2',
        'HLA-DRB5','HLA-DRB6','HLA-DQA2','HLA-DQB2')

plot_list = list()
for (i in (1:length(genes))) {
  gene=levels(factor(genes))[i]
  dat.m1 <- melt(surv_df,
                 id.vars=c('progression_free','progression'), 
                 measure.vars=c(gene))
  # 2. Create the plot
  fit <- survfit( Surv(progression_free,progression) ~ value, data = dat.m1 )
  ggsurv=ggsurvplot_facet(fit, dat.m1, facet.by = c("variable"),legend.labs=c('low','high'),short.panel.labs=T,
                          palette = "jco", pval = T)
  plot_list[[i]]=ggsurv
  message(gene)
}

# Save plots to tiff. Makes a separate file for each plot.
for (i in (1:length(genes))) {
  ggsave(plot_list[[i]], file=paste0('VanAllen',levels(factor(genes))[i],".png"), width = 10, height = 8, units = "cm")
}




# Fig. 2b
C1=subset(cohortVanAllen_df,cohortVanAllen_df$cluster == 1)
C2=subset(cohortVanAllen_df,cohortVanAllen_df$cluster == 2)
C3=subset(cohortVanAllen_df,cohortVanAllen_df$cluster == 3)
C4=subset(cohortVanAllen_df,cohortVanAllen_df$cluster == 4)
dfList=list(C1, C2, C3, C4)
HLA_supertypes=lapply(dfList, function(x) {
  supertypes = data.frame(mean(log2(x$HLA_A+1),na_rm=T), mean(log2(x$HLA_B+1),na_rm=T), mean(log2(x$HLA_C+1),na_rm=T), mean(log2(x$HLA_E+1),na_rm=T), mean(log2(x$HLA_F+1),na_rm=T), mean(log2(x$HLA_G+1),na_rm=T), mean(log2(x$B2M+1),na_rm=T), mean(log2(x$HLA_DRA+1),na_rm=T), mean(log2(x$HLA_DRB1+1),na_rm=T), mean(log2(x$HLA_DQA1+1),na_rm=T), mean(log2(x$HLA_DQB1+1),na_rm=T), mean(log2(x$HLA_DPA1+1),na_rm=T), mean(log2(x$HLA_DPB1+1),na_rm=T),mean(log2(x$HLA_DMA+1),na_rm=T), mean(log2(x$HLA_DMB+1),na_rm=T), mean(log2(x$HLA_DOA+1),na_rm=T), mean(log2(x$HLA_DOB+1),na_rm=T),            
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
                   border_color=NA, legend=T, angle_col = "315",main='Van Allen Melanoma Anti-CTLA4')


# Fig. 2d
# immune signatures
data=cohortVanAllen_df
data$`CD8 T cells`=data$T.cells.CD8
data$`CD4 Tmem`=data$T.cells.CD4.memory.activated
data$`M1`=data$Macrophages.M1
data$`M2`=data$Macrophages.M2
data$`Tfh`=data$T.cells.follicular.helper
data$`NK cells`=data$NK.cells.activated
data$Tregs=data$T.cells.regulatory

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
  theme_bw() + ggtitle('Van Allen Melanoma Anti-CTLA4')


# Fig. 2f
# exhaustion signatures
data=cohortVanAllen_df
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



# # Fig. 2h hazard ratio
cohortVanAllen_df$CYT=(log2(cohortVanAllen_df$GZMA + 1) + log2(cohortVanAllen_df$PRF1 + 1)) / 2
cohortVanAllen_df$highCYT=cohortVanAllen_df$CYT > unname(quantile(cohortVanAllen_df$CYT,na.rm=TRUE)["75%"])
cohortVanAllen_df$highTMB=cohortVanAllen_df$nonsyn_muts > 100
cohortVanAllen_df$C1=F
cohortVanAllen_df$C1[which(cohortVanAllen_df$cluster == 1)]=T
cohortVanAllen_df$C2=F
cohortVanAllen_df$C2[which(cohortVanAllen_df$cluster == 2)]=T
explanatory = c('highCYT','highTMB',"C1",'C2')
dependent = "Surv(progression_free,progression)"
# 710 * 310
cohortVanAllen_df %>%
  hr_plot(dependent, explanatory, remove_ref = T,breaks=c(0.1,1,5),
          column_space = c(-0.4, -0.2, 0.8, 1), dependent_label = "", prefix = "", suffix = "", table_text_size = 3, title_text_size = 18)




# Extended Data Fig. 2b
# IFN response
dat.m1 <- melt(cohortVanAllen_df %>% filter(!is.na(cluster)),
               id.vars='cluster', 
               measure.vars=c('SSGSEA_IFNA_RESPONSE','SSGSEA_IFNG_RESPONSE'))
dat.m1$label=NA
dat.m1$label[which(dat.m1$variable == 'SSGSEA_IFNA_RESPONSE')]='IFN-alpha Response'
dat.m1$label[which(dat.m1$variable == 'SSGSEA_IFNG_RESPONSE')]='IFN-gamma Response'

my_comparisons=list(c('C1','C2'),c('C3','C4'),c('C1','C3'))
dat.m1$Cluster=NA
dat.m1$Cluster[which(dat.m1$cluster == 1)]='C1'
dat.m1$Cluster[which(dat.m1$cluster == 2)]='C2'
dat.m1$Cluster[which(dat.m1$cluster == 3)]='C3'
dat.m1$Cluster[which(dat.m1$cluster == 4)]='C4'
# 400 * 500
ggplot(dat.m1, aes(x=Cluster, y=value, fill=Cluster)) +
  geom_violin() +
  stat_compare_means(comparisons = my_comparisons, size=3) + 
  geom_boxplot(width=0.1, fill="white") + xlab('') + ylab('Enrichment Score') +
  guides(fill=FALSE) +
  facet_grid(cols = vars(label),scales = "free") +
  scale_fill_manual(values = c('#CD9B1D',"#EAD6A3", "#BFB1D6",'#5E3C99')) +
  theme_bw() + ggtitle('Van Allen Melanoma Anti-CTLA4')



# Extended Data Fig. 2d
# purity
dat.m1 <- melt(cohortVanAllen_df %>% filter(!is.na(cluster)),
               id.vars='cluster', 
               measure.vars=c('purity'))

my_comparisons=list(c('C1','C2'),c('C3','C4'),c('C1','C3'))
dat.m1$Cluster=NA
dat.m1$Cluster[which(dat.m1$cluster == 1)]='C1'
dat.m1$Cluster[which(dat.m1$cluster == 2)]='C2'
dat.m1$Cluster[which(dat.m1$cluster == 3)]='C3'
dat.m1$Cluster[which(dat.m1$cluster == 4)]='C4'
dat.m1$variable='Tumor Purity'
# 350 * 500
ggplot(cohortVanAllen_df, aes(x=Cluster, y=purity, fill=Cluster)) +
  geom_violin() +
  stat_compare_means(comparisons = my_comparisons, size=3) + 
  geom_boxplot(width=0.1, fill="white") + xlab('') + ylab('Tumor Purity') +
  guides(fill=FALSE) +
  scale_fill_manual(values = c('#CD9B1D',"#EAD6A3", "#BFB1D6",'#5E3C99')) +
  theme_bw() + ggtitle('Van Allen Melanoma Anti-CTLA4')



# Extended Data Fig. 2f
# mutation
dat.m1 <- melt(cohortVanAllen_df %>% filter(!is.na(cluster)),
               id.vars='cluster', 
               measure.vars=c('nonsyn_muts'))
dat.m1$label='TMB'
dat.m1$value=log10(dat.m1$value)
my_comparisons=list(c('C1','C2'),c('C3','C4'),c('C1','C3'))
dat.m1$Cluster=NA
dat.m1$Cluster[which(dat.m1$cluster == 1)]='C1'
dat.m1$Cluster[which(dat.m1$cluster == 2)]='C2'
dat.m1$Cluster[which(dat.m1$cluster == 3)]='C3'
dat.m1$Cluster[which(dat.m1$cluster == 4)]='C4'
# 350 * 500
cohortVanAllen_df$TMB=log10(cohortVanAllen_df$nonsyn_muts)
ggplot(cohortVanAllen_df, aes(x=Cluster, y=TMB, fill=Cluster)) +
  geom_violin() +
  stat_compare_means(label.x=1,label.y=4) + 
  geom_boxplot(width=0.1, fill="white") + xlab('') + ylab('TMB') +
  guides(fill=FALSE) +
  scale_fill_manual(values = c('#CD9B1D',"#EAD6A3", "#BFB1D6",'#5E3C99')) +
  theme_bw() + ggtitle('Van Allen Melanoma Anti-CTLA4')





### Figs. for Methods
# Extended Data Fig. 6b
df=subset(cohortVanAllen_df,select=(c('patientId','BOR',"progression_free","progression",'B2M','HLA_A','HLA_B','HLA_C','HLA_E','HLA_G','HLA_F','CIITA',"HLA_DRA", "HLA_DRB1","HLA_DQA1","HLA_DQB1","HLA_DPA1","HLA_DPB1","HLA_DMA",'HLA_DMB',"HLA_DOA",'HLA_DOB')))

CRPR_df=subset(df,df$BOR == 'CR' | df$BOR == 'PR')
SD_df=subset(df,df$BOR == 'SD')
PD_df=subset(df,df$BOR == 'PD')
X_df=subset(df,df$BOR == 'X')

df=rbind(CRPR_df,SD_df,PD_df, X_df)
patientList=df$patientId
rownames(df)=patientList
df[5:ncol(df)]=log2(df[5:ncol(df)]+1)
df=df[5:ncol(df)]
data=as.matrix(df)

# Compute hierarchical clustering
res.hc <- df %>%
  dist(method = "euclidean") %>% # Compute dissimilarity matrix
  hclust(method = "complete")     # Compute hierachical clustering

# Visualize using factoextra
# Cut in 4 groups and color by groups
# 1200 * 300
fviz_dend(res.hc, k = 4, # Cut in four groups
          cex = 0.5, # label size
          k_colors = c("#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07"),
          color_labels_by_k = TRUE, # color labels by groups
          rect = TRUE # Add rectangle around groups
)
# 450 * 470
out_38=pheatmap::pheatmap(data, color = my_palette, cluster_cols=T,cluster_rows=T, angle_col = "315",
                          annotation_names_row = F,annotation_names_col = F, show_rownames = F, 
                          legend=T, fontsize = 8,main='Base Model: Van Allen Cohort')





