# KM survival across APM clusters
cohort_df$Cluster=paste('C',cohort_df$cluster,sep='')
cohort_df$Cluster=factor(cohort_df$Cluster)
fit <- survfit(Surv(progression_free,progression) ~ Cluster,data = cohort_df)
ggsurvplot(fit, data = cohort_df, palette = c('#CD9B1D',"#EAD6A3", "#BFB1D6",'#5E3C99'),
           conf.int = F, 
           pval = T,  
           risk.table = F,
           xlab = '', ylab = 'Survival', 
           legend.labs = levels(factor(cohort_df$Cluster)),
           ggtheme = theme_bw()+theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                      panel.background = element_blank()))
ggsave('survival.png',height = 6,width =4.5, dpi=300)



# clinical benefit barplots across APM clusters
facet_df=cohort_df
facet_df$`clinical benefit`=NA
facet_df$`clinical benefit`[which(facet_df$BOR == 'PD')]='NCB'
facet_df$`clinical benefit`[which(facet_df$BOR != 'PD')]='DCB'
facet_df$`clinical benefit`[which(facet_df$BOR == 'MR' & facet_df$progression_free < 180)]='NCB'

C1=subset(facet_df,facet_df$cluster== 1)
C1_percent=paste(round(length(which(C1$`clinical benefit` == 'DCB')) / nrow(C1) * 100,0), '%',sep='')
C2=subset(facet_df,facet_df$cluster== 2)
C2_percent=paste(round(length(which(C2$`clinical benefit` == 'DCB')) / nrow(C2) * 100,0), '%',sep='')
C3=subset(facet_df,facet_df$cluster== 3)
C3_percent=paste(round(length(which(C3$`clinical benefit` == 'DCB')) / nrow(C3) * 100,0), '%',sep='')
C4=subset(facet_df,facet_df$cluster== 4)
C4_percent=paste(round(length(which(C4$`clinical benefit` == 'DCB')) / nrow(C4) * 100,0), '%',sep='')

theme_set(theme_classic((base_size=15)))
ggplot(facet_df, aes(Cluster,fill=`clinical benefit`)) + 
  geom_bar(position = "fill") +
  scale_fill_manual(values=c("gold2", "gray24")) + ylim(c(0,1.1)) +
  xlab("") + ylab('Fraction') + theme_bw() + theme(text = element_text(size=12)) +
  annotate(geom = "text", y = c(1,1,1,1) - 0.1, x = c(1,2,3,4)-0.2, label = c(C1_percent,C2_percent,C3_percent,C4_percent), 
           size = 3, hjust = 0)
ggsave('bar.png',height = 6,width =4.7, dpi=300)


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
df_combined=rbind(C1,C3)
fisher.test(df_combined,alternative = "two.sided")




# survival by PD1 expression in APM C1 and C2
c1_c2=subset(cohort_df,cohort_df$cluster == 1 | cohort_df$cluster == 2)
c1_c2$highPD1=c1_c2$PDCD1 > unname(quantile(c1_c2$PDCD1,na.rm=TRUE)["50%"])
fit <- survfit(Surv(progression_free,progression) ~ highPD1,data = c1_c2)
# 430 * 550
ggsurvplot(fit, data = c1_c2, palette = c("deepskyblue3", "firebrick3"),
           conf.int = F, 
           title = "Liu Cohort C1 + C2",
           pval = T,  
           risk.table = F,
           xlab = '', ylab = 'Survival', 
           legend.labs = c('Low PDCD1','High PDCD1'),
           ggtheme = theme_bw()+theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                      panel.background = element_blank()))
ggsave('survival.png',height = 5,width =4.5, dpi=300)


c3_c4=subset(cohort_df,cohort_df$cluster == 3 | cohort_df$cluster == 4)
c3_c4$highPD1=c3_c4$PDCD1 > unname(quantile(c3_c4$PDCD1,na.rm=TRUE)["50%"])
fit <- survfit(Surv(progression_free,progression) ~ highPD1,data = c3_c4)
# 430 * 550
ggsurvplot(fit, data = c3_c4, palette = c("deepskyblue3", "firebrick3"),
           conf.int = F, 
           title = "Liu Cohort C3 + C4",
           pval = T,  
           risk.table = F,
           xlab = '', ylab = 'Survival', 
           legend.labs = c('Low PDCD1','High PDCD1'),
           ggtheme = theme_bw()+theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                      panel.background = element_blank()))
ggsave('survival.png',height = 5,width =4.5, dpi=300)


cohort_df$meta.cluster=NA
cohort_df$meta.cluster[which(cohort_df$Cluster == 'C1' | cohort_df$Cluster == 'C2')]='C1+C2'
cohort_df$meta.cluster[which(cohort_df$Cluster == 'C3' | cohort_df$Cluster == 'C4')]='C3+C4'
cohort_df$PDCD1_log2=log2(cohort_df$PDCD1+1)
ggplot(cohort_df, aes(x=meta.cluster, y=PDCD1_log2, fill=meta.cluster)) +
  geom_violin() +
  stat_compare_means(size=4, label.x=1.4,label.y=3.9) + 
  geom_boxplot(width=0.1, fill="white") + 
  guides(fill=FALSE) + xlab('') + ylab('PDCD1 Expression') +
  theme_bw() +
  scale_fill_manual(values = c('firebrick3','deepskyblue3'))
ggsave('pd1.png',height = 4,width =3.5, dpi=300)

