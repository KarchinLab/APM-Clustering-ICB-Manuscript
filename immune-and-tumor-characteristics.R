# immune infiltration
data=cohort_df
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

ggplot(dat.m1, aes(x=Cluster, y=value, fill=Cluster)) +
  geom_violin() +
  stat_compare_means(comparisons = my_comparisons, size=3) + 
  geom_boxplot(width=0.1, fill="white") + 
  guides(fill=FALSE) +
  facet_grid(cols = vars(variable),scales = "free") + xlab('') + ylab('Immune Infiltration') +
  scale_fill_manual(values = c('#CD9B1D',"#EAD6A3", "#BFB1D6",'#5E3C99')) +
  theme_bw()
ggsave('Immune-Infiltration.png',height = 5.5,width =8, dpi=300)




# TMB
dat.m1 <- melt(cohort_df,
               id.vars='Cluster', 
               measure.vars=c('nonsyn_muts'))
dat.m1$value=log10(dat.m1$value)
my_comparisons=list(c('C1','C2'),c('C3','C4'),c('C1','C3'))

cohort_df$TMB=log10(cohort_df$nonsyn_muts)
ggplot(cohort_df, aes(x=Cluster, y=TMB, fill=Cluster)) +
  geom_violin() +
  stat_compare_means(label.x=1,label.y=4, size=4) + 
  geom_boxplot(width=0.1, fill="white") + xlab('') + ylab('TMB') +
  guides(fill=FALSE) +
  scale_fill_manual(values = c('#CD9B1D',"#EAD6A3", "#BFB1D6",'#5E3C99')) +
  theme_bw()  + theme(axis.text = element_text(size = 12)) + 
  ggtitle('Liu Melanoma Anti-PD1')
ggsave('tmb.png',height = 4,width =5.5, dpi=300)




# tumor-intrinsic characteristics
dat.m1 <- melt(cohort_df %>% filter(!is.na(Cluster)),
               id.vars='Cluster', 
               measure.vars=c("SSGSEA_PROLIFERATION","SSGSEA_ANGIOGENESIS","SSGSEA_TGFB_SIGNALING",
                              "SSGSEA_HYPOXIA","SSGSEA_ROS"))
dat.m1$label=NA
dat.m1$label[which(dat.m1$variable == 'SSGSEA_TGFB_SIGNALING')]='TGF-β Signaling'
dat.m1$label[which(dat.m1$variable == 'SSGSEA_HYPOXIA')]='Hypoxia'
dat.m1$label[which(dat.m1$variable == 'SSGSEA_PROLIFERATION')]='Proliferation'
dat.m1$label[which(dat.m1$variable == 'SSGSEA_ANGIOGENESIS')]='Angiogenesis'
dat.m1$label[which(dat.m1$variable == 'SSGSEA_ROS')]='ROS Pathway'

my_comparisons=list(c('C1','C2'),c('C3','C4'),c('C1','C3'))
ggplot(dat.m1, aes(x=Cluster, y=value, fill=Cluster)) +
  geom_violin() +
  stat_compare_means(comparisons = my_comparisons, size=3) + 
  geom_boxplot(width=0.1, fill="white") + xlab('') + ylab('Enrichment Score') +
  guides(fill=FALSE) +
  facet_grid(cols = vars(label),scales = "free") +
  scale_fill_manual(values = c('#CD9B1D',"#EAD6A3", "#BFB1D6",'#5E3C99')) +
  theme_bw() + ggtitle('Liu Melanoma Anti-PD1')
ggsave('ssgsea.png',height = 4,width =6, dpi=300)

