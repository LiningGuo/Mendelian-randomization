#hipps outa formated~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list=ls())
library(TwoSampleMR)
library(data.table)
fre_info <- fread('/gpfs/lab/groupYU/members/guolining/GLN/MR_ana/variants.txt')
head(fre_info)
fre_info <- fre_info[,c(-1,-3,-4,-5,-7)]

folder_path <- "/gpfs/lab/groupYU/members/guolining/GLN/hipp_sub"
txt_files <- list.files(path = folder_path, pattern = "\\.txt$", full.names = TRUE)
txt_files1 <- txt_files[c(13:15)]
txt_files1 <- txt_files[c(21)]

for (file_path in txt_files1) {
  outa_dat <- fread(file_path)
  outa_dat$N <- 33224
  M <- merge(fre_info,outa_dat,by.x="rsid",by.y="rsid")
  outa_format <- format_data(
    M,
    type='outcome',
    snp_col = "rsid",
    beta_col = "beta",
    se_col = "se",
    effect_allele_col ="A1",
    other_allele_col = "A2",
    pval_col = "p",
    eaf_col = "af",
    samplesize_col = "N"
  )
  filename <- substring(file_path, first = 50, last = 54)
  write.table(outa_format,row.names = F,quote = F,file = paste0("/gpfs/lab/groupYU/members/guolining/GLN/MR_ana/hipps_out_format/",filename,'format1.txt'))
  rm('outa_dat','M','outa_format')
}
