#dis expo formated~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list=ls())
library(TwoSampleMR)
library(data.table)
fre_info <- fread('/gpfs/lab/groupYU/members/guolining/GLN/MR_ana/variants.txt')
head(fre_info)
fre_info <- fre_info[,c(-1,-3,-4,-5,-7)]

dis <- fread('/gpfs/lab/groupYU/members/guolining/GLN/dise_gwas/SCZ/PGC3_SCZ_wave3.european.autosome.public.v3.vcf1.tsv')
head(dis)
#dis$beta <- log(dis$OR)
#dis$ncas <- 15212
#dis$ncon <- 29677
dis_s <- subset(dis,dis$`PVAL`<5e-8)
expo_format <- format_data(
  dis_s,
  type='exposure',
  snp_col = "ID",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col ="A1",
  other_allele_col = "A2",
  pval_col = "PVAL",
  eaf_col = "FCAS",
  #samplesize_col = "N",
  ncase_col = "NCAS",
  ncontrol_col = "NCON"
)
expo_sig <- clump_data(expo_format,clump_kb = 10000,clump_r2=0.001,pop = "EUR")
write.table(expo_sig,row.names = F,quote = F,file = "/gpfs/lab/groupYU/members/guolining/GLN/MR_ana/dis_expo_format/SCZ_format.txt")

