cohort_df$CYT=(log2(cohort_df$GZMA + 1) + log2(cohort_df$PRF1 + 1)) / 2

coxph_df=cohort_df
coxph_df$`*APM Clusters`=as.numeric(coxph_df$cluster)
coxph_df$`Mutational Subtype`=coxph_df$Mutational..Subtype
coxph_df$Stage=coxph_df$M.Stage
coxph_df$Stage[which(is.na(coxph_df$Stage) | coxph_df$Stage == 'UNKNOWN')]='Unknown or NA'
coxph_df$TMB=coxph_df$Mutation.Load
coxph_df$`Prior ipi`=0
coxph_df$`Prior ipi`[which(coxph_df$pretreatment == T)]=1
coxph_df$`CD8 T Cells`=coxph_df$T.cells.CD8
coxph_df$M2=coxph_df$Macrophages.M2
coxph_df$`*CYT`=(log2(coxph_df$GZMA + 1) + log2(coxph_df$PRF1 + 1)) / 2
coxph_df$`*TIDE Score`=coxph_df$tide_score
model <- coxph( Surv(OS,OS_status) ~ `*APM Clusters` + TMB + PDCD1 + `CD8 T Cells` + M2 + `Mutational Subtype` + Stage + `Prior ipi` +
                  `*TIDE Score` + `*CYT`,
                data = coxph_df )
ggforest(model, main = "",
         fontsize = 0.8,
         noDigits = 2)
ggsave("coxph.png",height = 6.5,width =7,dpi = 300)


coxph_df$`APM Clusters`=as.numeric(coxph_df$cluster)
coxph_df$`CYT`=(log2(coxph_df$GZMA + 1) + log2(coxph_df$PRF1 + 1)) / 2
coxph_df$`TIDE Score`=coxph_df$tide_score
model <- coxph( Surv(OS,OS_status) ~ `APM Clusters` + `TIDE Score` + `CYT`,
                data = coxph_df )
ggforest(model, main = "",
         fontsize = 0.8,
         noDigits = 2)
ggsave("coxph.png",height = 3,width =5,dpi = 300)
