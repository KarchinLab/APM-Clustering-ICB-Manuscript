# cohortRiaz_df: data frame of transcriptome (in log2(TPM+1)) & clinical data in the Liu cohort. 

# housekeeping gene score: for normalization
cohortRiaz_df$housekeeping=(cohortRiaz_df$ACTB+cohortRiaz_df$GAPDH+cohortRiaz_df$SDHA+cohortRiaz_df$HPRT1+cohortRiaz_df$HBS1L+cohortRiaz_df$AHSP)/6
housekeeping_mean=mean(cohortRiaz_df$housekeeping,na.rm=T)

# Fig. 4a
cohortRiaz_df$Cluster=paste('C',cohortRiaz_df$cluster,sep='')
cohortRiaz_df$Cluster=factor(cohortRiaz_df$Cluster)
fit <- survfit(Surv(OS,OS_status) ~ Cluster,data = cohortRiaz_df)
ggsurvplot(fit, data = cohortRiaz_df, palette = c('#CD9B1D','#EAD6A3','#BFB1D6','#5E3C99'),
           conf.int = F, 
           pval = T,  
           risk.table = TRUE,
           legend.labs = levels(factor(cohortRiaz_df$Cluster)),
           xlab = '', ylab = 'Survival',title='Riaz anti-PD1 (validation)',
           ggtheme = theme_bw()+theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                      panel.background = element_blank()))
coxph(Surv(OS,OS_status) ~ cluster, data = cohortRiaz_df) %>% gtsummary::tbl_regression(exp = TRUE)



# Fig. 4b: hazard ratio
cohortRiaz_df$PD1_PDL1=(cohortRiaz_df$PDCD1 + cohortRiaz_df$CD274) / 2
cohortRiaz_df$highPD1_PDL1=cohortRiaz_df$PD1_PDL1 > unname(quantile(cohortRiaz_df$PD1_PDL1,na.rm=TRUE)["75%"])
cohortRiaz_df$highTMB=cohortRiaz_df$Mutation.Load > 100
cohortRiaz_df$C1=F
cohortRiaz_df$C1[which(cohortRiaz_df$cluster == 1)]=T
cohortRiaz_df$C2=F
cohortRiaz_df$C2[which(cohortRiaz_df$cluster == 2)]=T
explanatory = c('highTMB','highPD1_PDL1',"C1",'C2')
dependent = "Surv(OS,OS_status)"
# 710 * 310
cohortRiaz_df %>%
  hr_plot(dependent, explanatory, remove_ref = T,
          column_space = c(-0.4, -0.2, 0.8, 1), dependent_label = "", 
          prefix = "", suffix = "", table_text_size = 3, title_text_size = 18, 
          breaks=c(0.2,0.5,1,1.5,2))



# Fig. 4c-e (ex. HLA.A: pre-treatment HLA-A expression; A_on: on-treatment HLA-A expression)
cohortRiaz_df$HLAI_pre=(cohortRiaz_df$HLA.A + cohortRiaz_df$HLA.B + cohortRiaz_df$HLA.C + cohortRiaz_df$B2M) / 4
cohortRiaz_df$HLAI_on=(cohortRiaz_df$A_on + cohortRiaz_df$B_on + cohortRiaz_df$C_on + cohortRiaz_df$B2M_on) / 4
cohortRiaz_df$HLAII_pre=(cohortRiaz_df$HLA.DRA + cohortRiaz_df$HLA.DRB1 + cohortRiaz_df$HLA.DQA1 + cohortRiaz_df$HLA.DQB1 + cohortRiaz_df$HLA.DPA1 + cohortRiaz_df$HLA.DPB1) / 6
cohortRiaz_df$HLAII_on=(cohortRiaz_df$DRA_on + cohortRiaz_df$DRB1_on + cohortRiaz_df$DQA1_on + cohortRiaz_df$DQB1_on + cohortRiaz_df$DPA1_on + cohortRiaz_df$DPB1_on) / 6
cohortRiaz_df$proteasome = (cohortRiaz_df$PSME1 + cohortRiaz_df$PSME2 + cohortRiaz_df$PSMB8 + cohortRiaz_df$PSMB9 + cohortRiaz_df$PSMB10) / 5
cohortRiaz_df$proteasome_on = (cohortRiaz_df$PSME1_on + cohortRiaz_df$PSME2_on + cohortRiaz_df$PSMB8_on + cohortRiaz_df$PSMB9_on + cohortRiaz_df$PSMB10_on) / 5
cohortRiaz_df$exhaustion_pre=(cohortRiaz_df$PDCD1 + cohortRiaz_df$CTLA4 + cohortRiaz_df$LAG3 + cohortRiaz_df$TIGIT + cohortRiaz_df$HAVCR2 + cohortRiaz_df$PDCD1LG2 + cohortRiaz_df$TIGIT) / 7
cohortRiaz_df$exhaustion_on=(cohortRiaz_df$PDCD1_on + cohortRiaz_df$CTLA4_on + cohortRiaz_df$LAG3_on + cohortRiaz_df$TIGIT_on + cohortRiaz_df$HAVCR2_on + cohortRiaz_df$PDCD1LG2_on + cohortRiaz_df$TIGIT_on) / 7

paired=subset(cohortRiaz_df, !is.na(cohortRiaz_df$A_on)) # filter out patients with no paired data
paired$Cluster=NA
paired$Cluster[which(paired$cluster == 1)] = 'C1'
paired$Cluster[which(paired$cluster == 2)] = 'C2'
paired$Cluster[which(paired$cluster == 3)] = 'C3'
paired$Cluster[which(paired$cluster == 4)] = 'C4'

paired$pre=paired$proteasome
paired$on=paired$proteasome_on
# APM differential expression between pre- and post-ICI treatment
# 450 * 500
ggpaired(paired, cond1 = "pre", cond2 = "on",
         fill = "condition", palette = "jco", facet.by='Cluster', xlab = '', ylab = 'Proteasome Expression') + 
  stat_compare_means(label.x = 1, label.y = 10, size=3, paired = TRUE) + ylim(c(4.5,10.5))


# Fig. 4f-j, Extended Data Fig. 4f-g (template: features can be substituted)
# Comparison of immune cell infiltration between pre- and post-ICI treatment
paired$pre=paired$T.cells.CD8
paired$on=paired$T.cells.CD8_on
# APM differential expression between pre- and post-ICI treatment
ggpaired(paired, cond1 = "pre", cond2 = "on",
         fill = "condition", palette = "jco", facet.by='Cluster', xlab = '', ylab = 'CD8 T Cell Infiltration') + 
  stat_compare_means(label.x = 1, label.y = 12, size=3, paired = TRUE) + coord_cartesian(ylim = c(0,15))


# Fig. 4k
# immunoinhibitory genes
paired$pre=paired$exhaustion_pre
paired$on=paired$exhaustion_on
ggpaired(paired, cond1 = "pre", cond2 = "on",
         fill = "condition", palette = "jco", facet.by='Cluster', xlab = '', ylab = 'Immunoinhibitor Expression') + 
  stat_compare_means(label.x = 1, label.y = 7.5, size=3, paired = TRUE) + coord_cartesian(ylim = c(0,8))



# Extended Data Fig. 4a
# 420 * 500
facet_df=cohortRiaz_df
facet_df$`clinical benefit`='NCB'
facet_df$`clinical benefit`[which(facet_df$BOR != 'PD')]='DCB'
facet_df$`clinical benefit`[which(facet_df$BOR == 'NE' & facet_df$OS < 180)]='NCB'

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
  xlab("") + ylab('Fraction') + theme_bw() + ggtitle('Riaz Cohort (validation)') + theme(text = element_text(size=8)) +
  annotate(geom = "text", y = c(1,1,1) - 0.1, x = c(1,2,3)-0.2, label = c(C1_percent,C2_percent,C3_percent), size = 3.5, hjust = 0)

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



# Extended Data Fig. 4b
C1=subset(cohortRiaz_df,cohortRiaz_df$cluster == 1)
C2=subset(cohortRiaz_df,cohortRiaz_df$cluster == 2)
C3=subset(cohortRiaz_df,cohortRiaz_df$cluster == 3)
C4=subset(cohortRiaz_df,cohortRiaz_df$cluster == 4)
dfList=list(C1, C2, C3, C4)
HLA_supertypes=lapply(dfList, function(x) {
  supertypes = data.frame(mean(x$HLA.A,na.rm=T), mean(x$HLA.B,na.rm=T), mean(x$HLA.C,na.rm=T), mean(x$HLA.E,na.rm=T), mean(x$HLA.F,na.rm=T), mean(x$HLA.G,na.rm=T), mean(x$B2M,na.rm=T), mean(x$HLA.DRA,na.rm=T), mean(x$HLA.DRB1,na.rm=T), mean(x$HLA.DQA1,na.rm=T), mean(x$HLA.DQB1,na.rm=T), mean(x$HLA.DPA1,na.rm=T), mean(x$HLA.DPB1,na.rm=T),mean(x$HLA.DMA,na.rm=T), mean(x$HLA.DMB,na.rm=T), mean(x$HLA.DOA,na.rm=T), mean(x$HLA.DOB,na.rm=T),            
                          mean(x$PSME1,na.rm=T), mean(x$PSME2,na.rm=T), mean(x$PSMB8,na.rm=T), mean(x$PSMB9,na.rm=T), mean(x$PSMB10,na.rm=T))          
  
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
# heatmap.2(data,col= my_palette, scale = 'none', trace = "none", density.info = "none",key.xlab="Expression",keysize = 1)
# 600 * 500
pheatmap::pheatmap(data, color = my_palette, cluster_cols=T,cluster_rows=F, annotation_colors = annotation_colors,
                   annotation_col=annotation_col, annotation_legend = T, 
                   annotation_names_row = F,annotation_names_col = F, show_rownames = T,
                   border_color=NA, legend=T, angle_col = "315",main='Riaz cohort (validation)')




# Extended Data Fig. 4c
# immune signatures
data=cohortRiaz_df
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
  scale_fill_manual(values = c('#EAD6A3',"#CD9B1D", "#BFB1D6",'#5E3C99')) +
  theme_bw() + ggtitle('Riaz cohort (validation)')


# Extended Data Fig. 4d
# proinflammatory
data=cohortRiaz_df
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


# Extended Data Fig. 4e
# exhaustion signatures
data=cohortRiaz_df
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
  scale_fill_manual(values = c('#EAD6A3',"#CD9B1D", "#BFB1D6",'#5E3C99')) +
  theme_bw()




# Extended Data Fig. 5a-o (template: features can be substituted)
# BOR
paired=subset(cohortRiaz_df_pre, !is.na(cohortRiaz_df_pre$A_on))
paired$clinical='NCB'
paired$clinical[which(paired$BOR != 'PD')]='DCB'
paired$clinical[which(paired$BOR == 'NE' & paired$OS < 180)]='NCB'

paired$pre=paired$proteasome
paired$on=paired$proteasome_on
# 400 * 500
ggpaired(paired, cond1 = "pre", cond2 = "on",
         fill = "condition", palette = "jco", facet.by='clinical', xlab = '', ylab='Proteasome Expression') +
  stat_compare_means(label.x = 1, label.y = 10, size=3.5, paired = TRUE)



