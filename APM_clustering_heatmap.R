df=subset(cohort_df,select=(c('patientId','BOR',"progression_free","progression","HLA_A","HLA_B","HLA_C","HLA_DRA", "HLA_DRB1","HLA_DQA1","HLA_DQB1","HLA_DPA1","HLA_DPB1",'PSME1','TAPBP','NLRC5','PSMB10','TAP2','HLA_DRB6','HLA_DQA2','HLA_DQB2','CIITA','HLA_E','HLA_G','HLA_F','HLA_DMA','HLA_DOB')))

CRPR_df=subset(df,df$BOR == 'CR' | df$BOR == 'PR')
SD_df=subset(df,df$BOR == 'SD')
PD_df=subset(df,df$BOR == 'PD')
MR_df=subset(df,df$BOR == 'MR')

df=rbind(CRPR_df,SD_df,PD_df,MR_df)
patientList=df$patientId
rownames(df)=patientList

df[5:ncol(df)]=log2(df[5:ncol(df)]+1)
df=df[5:27]
colnames(df)=c("HLA-A","HLA-B","HLA-C","HLA-DRA", "HLA-DRB1","HLA-DQA1","HLA-DQB1","HLA-DPA1","HLA-DPB1",'PSME1','TAPBP','NLRC5','PSMB10','TAP2','HLA-DRB6','HLA-DQA2','HLA-DQB2','CIITA','HLA-E','HLA-G','HLA-F','HLA-DMA','HLA-DOB')
data=as.matrix(df)


annotation_row <- data.frame(row.names = patientList, 
                             'BOR' = c(rep("CRPR", nrow(CRPR_df)),rep("SD", nrow(SD_df)), rep("PD", nrow(PD_df)), rep('MR', nrow(MR_df))))


annotation_row$patientId=rownames(annotation_row)
rownames(annotation_row)=annotation_row$patientId
annotation_row=annotation_row[1:2]

annotation_colors = list(
  'BOR' = c("CRPR" = "#0073C2FF","SD" = '#EFC000FF','MR' = 'snow2',"PD" = 'gray35')
)
my_palette <- colorRampPalette(c("deepskyblue4","gray85","firebrick3"))(200)

png(filename = "heatmap.png",width = 1900, height = 2000,res = 300)
out_122=pheatmap::pheatmap(data, color = my_palette, cluster_cols=T,cluster_rows=T, clustering_distance_rows = "euclidean",
                           clustering_method = 'complete', show_rownames = F,
                           annotation_row=annotation_row,annotation_colors=annotation_colors,
                           legend=T, fontsize = 8, angle_col = "315",annotation_names_row = F)
dev.off()





# APM base model
df=subset(cohort_df,select=(c('patientId','BOR',"progression_free","progression","HLA_A","HLA_B","HLA_C",'HLA_E','HLA_G','HLA_F','NLRC5',"HLA_DRA", "HLA_DRB1","HLA_DQA1","HLA_DQB1","HLA_DQA2","HLA_DQB2","HLA_DPA1","HLA_DPB1",'CIITA')))

CRPR_df=subset(df,df$BOR == 'CR' | df$BOR == 'PR')
SD_df=subset(df,df$BOR == 'SD')
PD_df=subset(df,df$BOR == 'PD')
MR_df=subset(df,df$BOR == 'MR')

df=rbind(CRPR_df,SD_df,PD_df,MR_df)
patientList=df$patientId
rownames(df)=patientList

df[5:ncol(df)]=log2(df[5:ncol(df)]+1)
df=df[5:ncol(df)]
colnames(df)=c("HLA-A","HLA-B","HLA-C",'HLA-E','HLA-G','HLA-F','NLRC5',"HLA-DRA", "HLA-DRB1","HLA-DQA1","HLA-DQB1","HLA-DQA2","HLA-DQB2","HLA-DPA1","HLA-DPB1",'CIITA')
data=as.matrix(df)


my_palette <- colorRampPalette(c("deepskyblue4","gray85","firebrick3"))(200)

png(filename = "heatmap.png",width = 1300, height = 1700,res = 300)
out_122=pheatmap::pheatmap(data, color = my_palette, cluster_cols=T,cluster_rows=T, clustering_distance_rows = "euclidean",
                           clustering_method = 'complete', show_rownames = F,
                           legend=T, fontsize = 8, angle_col = "315",annotation_names_row = F,
                           main='Base Model: Liu Cohort')
dev.off()

