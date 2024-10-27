#hipps_dis~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list=ls())
library(TwoSampleMR)
library(data.table)

expo_sig <- fread('/gpfs/lab/groupYU/members/guolining/GLN/MR_ana/hipps_expo_format/0270_format.txt')

outa_dat <- fread('/gpfs/lab/groupYU/members/guolining/GLN/MR_ana/dis_out_format/AD_format1.txt')

#head(outa_dat)
mydata <- harmonise_data(expo_sig,outa_dat,action= 2)
#check procy SNP
miss_SNP_index <- which(expo_sig$SNP %in% mydata$SNP == F)
miss_SNP <- expo_sig$SNP[miss_SNP_index]
rm('miss_SNP_index')
# find proxy by way1
library(httr)
library(jsonlite)
out_proxy1 <- list()
server <- "http://grch37.rest.ensembl.org"
for (k in 1:dim(as.data.frame(miss_SNP))[1]) {
  na_1 <- sprintf("/ld/human/%s/1000GENOMES:phase_3:EUR?r2=0.8;window_size=500",miss_SNP[k])
  proxy1 <- GET(paste(server, na_1, sep = ""), content_type("application/json"))
  stop_for_status(proxy1)
  out_proxy1[[k]] <- fromJSON(toJSON(content(proxy1)))
}
out_proxy1_1 <- do.call('rbind',out_proxy1)
out_proxy1_1$variation1 <- unlist(out_proxy1_1$variation1)
out_proxy1_1$variation2 <- unlist(out_proxy1_1$variation2)
out_proxy1_1$population_name <- unlist(out_proxy1_1$population_name)
out_proxy1_1$d_prime <- unlist(out_proxy1_1$d_prime)
out_proxy1_1$r2 <- unlist(out_proxy1_1$r2)

snp_info <- data.frame()
for (i in 1:dim(as.data.frame(miss_SNP))[1]) {
  index_snp_miss <- which(out_proxy1_1$variation1 %in% miss_SNP[i])
  snp_info[i,1] <- miss_SNP[i]
  index_snp_in_outdata <- which(outa_dat$SNP %in% out_proxy1_1$variation2[index_snp_miss])
  if (length(index_snp_in_outdata) == 0){
    next
  }
  snp_in <- outa_dat$SNP[index_snp_in_outdata]
  index_in_proxy <- which(out_proxy1_1$variation2 %in% snp_in)
  index_r2max <- which(out_proxy1_1$r2[index_in_proxy] == max(out_proxy1_1$r2[index_in_proxy]))
  index_yes <- index_r2max[1]
  snp_yes <- out_proxy1_1$variation2[index_in_proxy][index_yes]
  index_in_outdata <- which(outa_dat$SNP %in% snp_yes)
  outa_dat[index_in_outdata,]$SNP <- miss_SNP[i]
  snp_info[i,2] <- snp_yes
}
colnames(snp_info) <- c('miss.snp','agency.snp')
mydata2 <- harmonise_data(expo_sig,outa_dat,action= 2)
rm('i','out_proxy1','proxy1','k','na_1','server','outa_dat')
rm(list=ls(pattern = 'index'))

mydata2 <- mydata

# 3. detect outliers
# outlier_det_loo <- mr_leaveoneout(har_dat)
library(dplyr)
library(RadialMR)
outlier_det_ste <- steiger_filtering(mydata2)
outlier_det_egg <- egger_radial(mydata2, 0.05, 1, TRUE)
outlier_det_ivw <- ivw_radial(mydata2, 0.05, 1, 0.0001, TRUE)

snp_outl_pal <- data.frame(mydata2$SNP[!mydata2$mr_keep]) %>% rename_with(~ 'SNP', 1)
snp_outl_ste <- data.frame(mydata2$SNP[!outlier_det_ste$steiger_dir]) %>% rename_with(~ 'SNP', 1)
snp_outl_egg <- data.frame(outlier_det_egg$data$SNP[which(outlier_det_egg$data$Outliers=='Outlier')]) %>% rename_with(~ 'SNP', 1)
snp_outl_ivw <- data.frame(outlier_det_ivw$data$SNP[which(outlier_det_ivw$data$Outliers=='Outlier')]) %>% rename_with(~ 'SNP', 1)

snp_outl_com <- unique(rbind(snp_outl_pal, snp_outl_ste, snp_outl_ivw, snp_outl_egg))
#snp_outl_com <- unique(rbind(snp_outl_pal, snp_outl_ste, snp_outl_egg))

har_dat <- anti_join(mydata2, snp_outl_com, by = 'SNP')

rm(list=ls(pattern = 'outlier'))
# Power Analysis
MAF <- pmin(har_dat$eaf.exposure,1-har_dat$eaf.exposure)
beta <- har_dat$beta.exposure
se <- har_dat$se.exposure
N <- har_dat$samplesize.exposure
R2 <- (2*MAF*(1-MAF)*beta^2)/(2*MAF*(1-MAF)*beta^2+2*MAF*(1-MAF)*N*se^2)
F <- (R2*(N-2))/(1-R2)
min(F)
rm('MAF','beta','se','N')

sum(R2)

res <- mr(har_dat,method_list = c('mr_two_sample_ml','mr_egger_regression',
                                 'mr_egger_regression_bootstrap','mr_simple_median',
                                 'mr_weighted_median','mr_penalised_weighted_median',
                                 'mr_ivw','mr_ivw_radial',
                                 'mr_ivw_mre','mr_ivw_fe',
                                 'mr_simple_mode','mr_weighted_mode',
                                 'mr_weighted_mode_nome','mr_simple_mode_nome',
                                 'mr_sign','mr_uwr'))
res <- mr(har_dat)
mr_or <- generate_odds_ratios(res)
rm('res')
# sensitivity analysis2 with robust methods2
library(mr.raps)
mr_sen_raps <- data.frame(mr.raps(har_dat$beta.exposure,
                                  har_dat$beta.outcome,
                                  har_dat$se.exposure,
                                  har_dat$se.outcome))
colnames(mr_sen_raps) <- c('b','se','pval','naive.se','chi.sq.test')
mr_sen_raps_or <- generate_odds_ratios(mr_sen_raps)
rm('mr_sen_raps')
# MRMix analysis (unstandardize)

mr_het <- mr_heterogeneity(har_dat,method_list = c('mr_two_sample_ml','mr_egger_regression',
                                                   'mr_ivw','mr_ivw_radial',
                                                   'mr_ivw_mre','mr_ivw_fe',
                                                   'mr_uwr'))

mr_het$I2 <- (mr_het$Q-mr_het$Q_df)/mr_het$Q

res_pleiotropy <- mr_pleiotropy_test(har_dat)

rm('outa_dat')

## MR-PRESSO
library(MRPRESSO)
#min(NbDistribution) = dim(har_data)[1]/0.05
result_presso <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", 
                           SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                           OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = har_dat, 
                           NbDistribution = 1000, SignifThreshold = 0.05)

presso_pval <- result_presso$`MR-PRESSO results`$`Global Test`$Pvalue
presso_outlier_SNP <- result_presso[["MR-PRESSO results"]][["Distortion Test"]][["Outliers Indices"]]

save.image("/gpfs/lab/groupYU/members/guolining/GLN/MR_ana/MR_result/0270_ad.RData")
rm(list=ls())
#plot
res <- mr(har_dat)
p1 <- mr_scatter_plot(res,har_dat)
print(p1)

#Leave-one-out
res_loo <- mr_leaveoneout(har_dat)
p2 <- mr_leaveoneout_plot(res_loo)
p2[[1]]

## MR-PRESSO
library(MRPRESSO)
result_presso <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", 
                           SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                           OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = har_dat, 
                           NbDistribution = 5000, SignifThreshold = 0.05)

presso_pval <- result_presso$`MR-PRESSO results`$`Global Test`$Pvalue
presso_outlier_SNP <- result_presso[["MR-PRESSO results"]][["Distortion Test"]][["Outliers Indices"]]

rm('outa_dat')



