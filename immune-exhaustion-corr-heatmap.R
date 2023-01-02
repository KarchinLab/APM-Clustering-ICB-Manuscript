exhaustion.signature=rev(c('PDCD1','CTLA4','HAVCR2','TIGIT','LAG3','CD96','BTLA'))
cluster_vec=c('C1','C2','C3','C4')

rho_df=data.frame()
pVal_df=data.frame()
for (exhaust.gene in exhaustion.signature) {
  pVal_vec=c()
  rho_vec=c()
  for (cluster_num in cluster_vec) {
    cluster_subset=expression_tpm_matrix[which(cohort_df$Cluster == cluster_num),]
    spearman_stats=cor.test(cluster_subset[,exhaust.gene], cluster_subset$T.cells.CD8)
    
    pVal_vec=append(pVal_vec,spearman_stats$p.value)
    rho_vec=append(rho_vec,spearman_stats$estimate)
  }
  temp_stats=data.frame('Exhaustion Signature'=exhaust.gene,'Cluster'=cluster_vec,'rho'=rho_vec)
  temp_pVal=data.frame('Exhaustion Signature'=exhaust.gene,'Cluster'=cluster_vec,'pval'=pVal_vec)
  
  rho_df=rbind(rho_df,temp_stats)
  pVal_df=rbind(pVal_df,temp_pVal)
}

plot.data <- rho_df
plot.pVal <- pVal_df
qVal=p.adjust(plot.pVal$pval,method="BH") 
plot.pVal$p.adj=qVal
plot.data$stars <- cut(plot.pVal$p.adj, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))  # Create column of significance labels

plot.data$Cluster=factor(plot.data$Cluster,rev(cluster_vec))
p <- ggplot(aes(x=Cluster, y=Exhaustion.Signature, fill=rho), data=plot.data)
p + geom_tile() + scale_fill_gradient2(low='dodgerblue4',mid='gray90',high='firebrick3') + 
  geom_text(aes(label=stars), color="black", size=5) + 
  labs(y='Exhaustion Signature', x=NULL, fill="Spearman's Rho") + coord_flip() + theme_bw() +
  theme(axis.text.x=element_text(angle = -45, hjust = 0),legend.title=element_text(size=8)) +
  guides(fill = guide_colourbar(barwidth = 1,
                                barheight = 10))
ggsave("corr.png",height = 3,width =6,dpi = 300)
