#hipps expo formated~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list=ls())
library(TwoSampleMR)
library(data.table)
fre_info <- fread('/gpfs/lab/groupYU/members/guolining/GLN/MR_ana/variants.txt')
head(fre_info)
fre_info <- fre_info[,c(-1,-3,-4,-5,-7)]

folder_path <- "/gpfs/lab/groupYU/members/guolining/GLN/hipp_sub"
txt_files <- list.files(path = folder_path, pattern = "\\.txt$", full.names = TRUE)
txt_files1 <- txt_files[c(43:44)]

for (file_path in txt_files1) {
  expo_dat <- fread(file_path)
  expo_s <- subset(expo_dat,expo_dat$p<5e-8)
  expo_s$N <- 33224
  M <- merge(fre_info,expo_s,by.x="rsid",by.y="rsid")
  expo_format <- format_data(
    M,
    type='exposure',
    snp_col = "rsid",
    beta_col = "beta",
    se_col = "se",
    effect_allele_col ="A1",
    other_allele_col = "A2",
    pval_col = "p",
    eaf_col = "af",
    samplesize_col = "N"
  )
  expo_sig <- clump_data(expo_format,clump_kb = 10000,clump_r2=0.001,pop = "EUR")
  
  filename <- substring(file_path, first = 50, last = 54)
  write.table(expo_sig,row.names = F,quote = F,file = paste0("/gpfs/lab/groupYU/members/guolining/GLN/MR_ana/hipps_expo_format/",filename,'format.txt'))
  rm('expo_sig')
  }




