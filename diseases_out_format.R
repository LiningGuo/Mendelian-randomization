#dis out formated~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list=ls())
library(TwoSampleMR)
library(data.table)
fre_info <- fread('/gpfs/lab/groupYU/members/guolining/GLN/MR_ana/variants.txt')
head(fre_info)
fre_info <- fre_info[,c(-1,-3,-4,-5,-7)]

dis <- fread('/gpfs/lab/groupYU/members/guolining/GLN/dise_gwas/AD/2021/PGCALZ2ExcludingUKBand23andME_METALInverseVariance_MetaAnalysis_SNP.txt')
head(dis)
#dis$beta <- log(dis$OR)
#dis$ncas <- 15212
#dis$ncon <- 29677
out_format <- format_data(
  dis,
  type='outcome',
  snp_col = "rsid",
  beta_col = "beta",
  se_col = "standard_error",
  effect_allele_col ="effect_allele",
  other_allele_col = "other_allele",
  pval_col = "p_value",
  eaf_col = "effect_allele_frequency",
  samplesize_col = "N"
  #ncase_col = "N_cases",
  #ncontrol_col = "N_controls"
)

write.table(out_format,row.names = F,quote = F,file = "/gpfs/lab/groupYU/members/guolining/GLN/MR_ana/dis_out_format/AD_format1.txt")

