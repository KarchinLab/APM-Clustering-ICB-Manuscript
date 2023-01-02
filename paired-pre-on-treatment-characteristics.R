# paired pre- and on-treatment characteristics across predicted APM clusters
cohort_riaz$HLAI_pre=(log2(cohort_riaz$HLA.A+1) + log2(cohort_riaz$HLA.B+1) + log2(cohort_riaz$HLA.C+1) + log2(cohort_riaz$HLA.E+1) + log2(cohort_riaz$HLA.G+1) + log2(cohort_riaz$HLA.F+1)) / 6
cohort_riaz$HLAI_on=(log2(cohort_riaz$A_on+1) + log2(cohort_riaz$B_on+1) + log2(cohort_riaz$C_on+1) + log2(cohort_riaz$E_on+1) + log2(cohort_riaz$G_on+1) + log2(cohort_riaz$F_on+1)) / 6
cohort_riaz$HLAII_pre=(log2(cohort_riaz$HLA.DRA+1) + log2(cohort_riaz$HLA.DRB1+1) + log2(cohort_riaz$HLA.DQA1+1) + log2(cohort_riaz$HLA.DQB1+1) + log2(cohort_riaz$HLA.DPA1+1) + log2(cohort_riaz$HLA.DPB1+1) + log2(cohort_riaz$HLA.DQA2+1) + log2(cohort_riaz$HLA.DQB2+1)) / 8
cohort_riaz$HLAII_on=(log2(cohort_riaz$DRA_on+1) + log2(cohort_riaz$DRB1_on+1) + log2(cohort_riaz$DQA1_on+1) + log2(cohort_riaz$DQB1_on+1) + log2(cohort_riaz$DPA1_on+1) + log2(cohort_riaz$DPB1_on+1) + log2(cohort_riaz$DQA2_on+1) + log2(cohort_riaz$DQB2_on+1)) / 8
cohort_riaz$proteasome = (log2(cohort_riaz$PSME1+1) + log2(cohort_riaz$PSME2+1) + log2(cohort_riaz$PSMB8+1) + log2(cohort_riaz$PSMB9+1) + log2(cohort_riaz$PSMB10+1)) / 5
cohort_riaz$proteasome_on = (log2(cohort_riaz$PSME1_on+1) + log2(cohort_riaz$PSME2_on+1) + log2(cohort_riaz$PSMB8_on+1) + log2(cohort_riaz$PSMB9_on+1) + log2(cohort_riaz$PSMB10_on+1)) / 5

paired=subset(cohort_riaz, !is.na(cohort_riaz$A_on))
ggpaired(paired, cond1 = "pre", cond2 = "on",
         fill = "condition", palette = "jco", facet.by='Cluster', xlab = '', ylab = 'HLA I Expression') + 
  stat_compare_means(size=3.5, paired = TRUE)
ggsave('paired.png',height = 4.5,width =3.7, dpi=300)



# paired pre- and on-treatment characteristics across clinically-defined patient responses
paired=subset(cohortRiaz_df_pre, !is.na(cohortRiaz_df_pre$A_on))
paired$clinical='NCB'
paired$clinical[which(paired$BOR != 'PD')]='DCB'
paired$clinical[which(paired$BOR == 'NE' & paired$OS < 180)]='NCB'

ggpaired(paired, cond1 = "pre", cond2 = "on",
         fill = "condition", palette = "jco", facet.by='clinical', xlab = '', ylab = 'HLA I Expression') + 
  stat_compare_means(size=3.5, paired = TRUE, label.x = 1, label.y = 13.2) + coord_cartesian(ylim = c(7,13.3))
ggsave('paired.png',height = 4.2,width =4, dpi=300)
