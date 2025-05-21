library(dplyr)
library(readr)
library(TwoSampleMR)
library(tidyverse)
library(ieugwasr)
library(MRPRESSO)

LADA.data <- read_tsv("D:/Others/药物和LADA/data/LADACTRL_MA_filtered_2018_DLC.txt", col_types=cols(.default="c"))
LADA.data[,5:10] <- sapply(LADA.data[,5:10], as.numeric)
LADA.data$beta <- log(LADA.data$OR)

LTL <- read_tsv(gzfile("D:/Others/药物和LADA/data/UKB_telomere_gwas_summarystats.tsv.gz"), col_types=cols(.default="c"))
LTL <- LTL[,c('variant_id', 'chromosome', 'base_pair_location')]
LTL$rs_number <- paste0(LTL$chromosome, ':', LTL$base_pair_location)
LTL <- LTL[,c('rs_number', 'variant_id')]

LADA.data <- merge(LADA.data, LTL, by='rs_number', all.x = T)
LADA.data <- LADA.data[LADA.data$P<5e-8, ]
dat <- data.frame(rsid=LADA.data$variant_id, pval=LADA.data$P)
retained_SNPs <- ld_clump(dat, clump_kb = 10000, clump_r2 = 0.001, bfile = 'D:/LD_reference/EUR',
                          plink_bin = 'D:/LD_reference/plink.exe')
LADA.data <- LADA.data[LADA.data$variant_id%in%retained_SNPs$rsid, ]

write.csv(LADA.data, 'LADA_IVs_5e-8.csv', row.names = F, quote = F)

LADA <- read.csv('LADA_IVs_5e-8.csv', stringsAsFactors=F, header=T)

LADA <- format_data(LADA, snp_col = 'variant_id', effect_allele_col = 'reference_allele',
                    other_allele_col = 'other_allele', eaf_col = 'eaf',
                    beta_col = 'beta', se_col = 'OR_se', pval_col = 'P')

mr_method <- c('mr_ivw', 'mr_ivw_radial', 'mr_weighted_mode', 'mr_raps', 
               'Cochran Q test', 'I^2', 'mr_egger_intercept_test')
exposure <- NULL
outcome <- NULL
beta <- NULL
beta.CI <- NULL
OR <- NULL
OR.CI <- NULL
p.value <- NULL
method <- NULL

#### SAA1 ####
SAA1 <- read_tsv(gzfile('E:/pQTL/LADA/15515_2_SAA1_SAA.txt.gz'), col_types=cols(.default="c"))
SAA1 <- SAA1[SAA1$rsids%in%LADA$SNP, ]
SAA1[,c('Beta', 'SE', 'Pval')] <- sapply(SAA1[,c('Beta', 'SE', 'Pval')], as.numeric)
SAA1 <- format_data(SAA1, type='outcome', snp_col = 'rsids', effect_allele_col = 'effectAllele',
                    other_allele_col = 'otherAllele', beta_col = 'Beta',
                    se_col = 'SE', pval_col = 'Pval')
LADA.SAA1 <- harmonise_data(LADA, SAA1)
LADA.SAA1.data <- data.frame(SNP=LADA.SAA1$SNP, EA=LADA.SAA1$effect_allele.exposure,
                             OA=LADA.SAA1$other_allele.exposure, EAF=LADA.SAA1$eaf.exposure,
                             beta.exposure=LADA.SAA1$beta.exposure, se.exposure=LADA.SAA1$se.exposure,
                             pval.exposure=LADA.SAA1$pval.exposure, beta.outcome=LADA.SAA1$beta.outcome,
                             se.outcome=LADA.SAA1$se.outcome, pval.outcome=LADA.SAA1$pval.outcome)
write.csv(LADA.SAA1.data, 'LADA_to_SAA1_data.csv', quote = F, row.names = F)


res <- mr(LADA.SAA1, method_list = c('mr_ivw', 'mr_ivw_fe', 'mr_weighted_median', 'mr_raps'));res
or.res <- generate_odds_ratios(res)

beta <- c(beta, round(-res$b,3))
beta.CI <- c(beta.CI, paste0(round(-or.res$lo_ci,3), ' to ', round(-or.res$up_ci,3)))
OR <- c(OR, round(exp(-res$b), 3))
OR.CI <- c(OR.CI, paste0(round(exp(-or.res$up_ci),3), ' to ', round(exp(-or.res$lo_ci),3)))
p.value <- c(p.value, res$pval)

LADA.SAA1 <- LADA.SAA1[!LADA.SAA1$SNP%in%'rs2516494', ]
radial_data <- RadialMR::format_radial(BXG=LADA.SAA1$beta.exposure, BYG = LADA.SAA1$beta.outcome,
                                       seBXG=LADA.SAA1$se.exposure, seBYG = LADA.SAA1$se.outcome,
                                       RSID=LADA.SAA1$SNP)
radial_ivw <- RadialMR::ivw_radial(radial_data, alpha=0.05, weights = 3, summary = T);radial_ivw

### heterogeneity test
heter <- mr_heterogeneity(LADA.SAA1, method_list = 'mr_ivw'); heter

beta <- c(beta, round(heter$Q, 3))
beta.CI <- c(beta.CI, '')
OR <- c(OR, '')
OR.CI <- c(OR.CI, '')
p.value <- c(p.value, heter$Q_pval)

### I^2
Q <- heter$Q
Q_df <- heter$Q_df
I.sqr <- max(0, (Q-Q_df)/Q*100);I.sqr

beta <- c(beta, paste0(round(I.sqr, 2), '%'))
beta.CI <- c(beta.CI, '')
OR <- c(OR, '')
OR.CI <- c(OR.CI, '')
p.value <- c(p.value, '')

### mr-egger intercept ###

intercept <- mr_pleiotropy_test(LADA.SAA1); intercept

beta <- c(beta, intercept$egger_intercept)
beta.CI <- c(beta.CI, '')
OR <- c(OR, '')
OR.CI <- c(OR.CI, '')
p.value <- c(p.value, intercept$pval)

exposure <- c(exposure, rep('LADA', 7))
outcome <- c(outcome, rep('SAA1', 7))
method <- c(method, mr_method)

#### SAA2 ####
SAA2 <- read_tsv(gzfile('E:/pQTL/LADA/18832_65_SAA2_SAA2.txt.gz'), col_types=cols(.default="c"))
SAA2 <- SAA2[SAA2$rsids%in%LADA$SNP, ]
SAA2[,c('Beta', 'SE', 'Pval')] <- sapply(SAA2[,c('Beta', 'SE', 'Pval')], as.numeric)
SAA2 <- format_data(SAA2, type='outcome', snp_col = 'rsids', effect_allele_col = 'effectAllele',
                    other_allele_col = 'otherAllele', beta_col = 'Beta',
                    se_col = 'SE', pval_col = 'Pval')
LADA.SAA2 <- harmonise_data(LADA, SAA2)
LADA.SAA2.data <- data.frame(SNP=LADA.SAA2$SNP, EA=LADA.SAA2$effect_allele.exposure,
                             OA=LADA.SAA2$other_allele.exposure, EAF=LADA.SAA2$eaf.exposure,
                             beta.exposure=LADA.SAA2$beta.exposure, se.exposure=LADA.SAA2$se.exposure,
                             pval.exposure=LADA.SAA2$pval.exposure, beta.outcome=LADA.SAA2$beta.outcome,
                             se.outcome=LADA.SAA2$se.outcome, pval.outcome=LADA.SAA2$pval.outcome)
write.csv(LADA.SAA2.data, 'LADA_to_SAA2_data.csv', quote = F, row.names = F)

res <- mr(LADA.SAA2, method_list = c('mr_ivw', 'mr_ivw_radial', 'mr_weighted_median', 'mr_raps'));res
or.res <- generate_odds_ratios(res)

LADA.SAA2 <- LADA.SAA2[!LADA.SAA2$SNP%in%'rs2516494', ]
radial_data <- RadialMR::format_radial(BXG=LADA.SAA2$beta.exposure, BYG = LADA.SAA2$beta.outcome,
                                       seBXG=LADA.SAA2$se.exposure, seBYG = LADA.SAA2$se.outcome,
                                       RSID=LADA.SAA2$SNP)
radial_ivw <- RadialMR::ivw_radial(radial_data, alpha=0.05, weights = 3, summary = T);radial_ivw

beta <- c(beta, round(-res$b,3))
beta.CI <- c(beta.CI, paste0(round(-or.res$lo_ci,3), ' to ', round(-or.res$up_ci,3)))
OR <- c(OR, round(exp(-res$b), 3))
OR.CI <- c(OR.CI, paste0(round(exp(-or.res$up_ci),3), ' to ', round(exp(-or.res$lo_ci),3)))
p.value <- c(p.value, res$pval)

### heterogeneity test
heter <- mr_heterogeneity(LADA.SAA2, method_list = 'mr_ivw'); heter

beta <- c(beta, round(heter$Q, 3))
beta.CI <- c(beta.CI, '')
OR <- c(OR, '')
OR.CI <- c(OR.CI, '')
p.value <- c(p.value, heter$Q_pval)

### I^2
Q <- heter$Q
Q_df <- heter$Q_df
I.sqr <- max(0, (Q-Q_df)/Q*100);I.sqr

beta <- c(beta, paste0(round(I.sqr, 2), '%'))
beta.CI <- c(beta.CI, '')
OR <- c(OR, '')
OR.CI <- c(OR.CI, '')
p.value <- c(p.value, '')

### mr-egger intercept ###

intercept <- mr_pleiotropy_test(LADA.SAA2); intercept

beta <- c(beta, intercept$egger_intercept)
beta.CI <- c(beta.CI, '')
OR <- c(OR, '')
OR.CI <- c(OR.CI, '')
p.value <- c(p.value, intercept$pval)

exposure <- c(exposure, rep('LADA', 7))
outcome <- c(outcome, rep('SAA2', 7))
method <- c(method, mr_method)


#### CXCL10 ####
CXCL10 <- read_tsv(gzfile('E:/pQTL/LADA/4141_79_CXCL10_IP_10.txt.gz'), col_types=cols(.default="c"))
CXCL10 <- CXCL10[CXCL10$rsids%in%LADA$SNP, ]
CXCL10[,c('Beta', 'SE', 'Pval')] <- sapply(CXCL10[,c('Beta', 'SE', 'Pval')], as.numeric)
CXCL10 <- format_data(CXCL10, type='outcome', snp_col = 'rsids', effect_allele_col = 'effectAllele',
                      other_allele_col = 'otherAllele', beta_col = 'Beta',
                      se_col = 'SE', pval_col = 'Pval')
LADA.CXCL10 <- harmonise_data(LADA, CXCL10)
LADA.CXCL10.data <- data.frame(SNP=LADA.CXCL10$SNP, EA=LADA.CXCL10$effect_allele.exposure,
                               OA=LADA.CXCL10$other_allele.exposure, EAF=LADA.CXCL10$eaf.exposure,
                               beta.exposure=LADA.CXCL10$beta.exposure, se.exposure=LADA.CXCL10$se.exposure,
                               pval.exposure=LADA.CXCL10$pval.exposure, beta.outcome=LADA.CXCL10$beta.outcome,
                               se.outcome=LADA.CXCL10$se.outcome, pval.outcome=LADA.CXCL10$pval.outcome)
write.csv(LADA.CXCL10.data, 'LADA_to_CXCL10_data.csv', quote = F, row.names = F)

res <- mr(LADA.CXCL10, method_list = c('mr_ivw', 'mr_ivw_radial', 'mr_weighted_median', 'mr_raps'));res
or.res <- generate_odds_ratios(res)

LADA.CXCL10 <- LADA.CXCL10[!LADA.CXCL10$SNP%in%'rs2516494', ]
radial_data <- RadialMR::format_radial(BXG=LADA.CXCL10$beta.exposure, BYG = LADA.CXCL10$beta.outcome,
                                       seBXG=LADA.CXCL10$se.exposure, seBYG = LADA.CXCL10$se.outcome,
                                       RSID=LADA.CXCL10$SNP)
radial_ivw <- RadialMR::ivw_radial(radial_data, alpha=0.05, weights = 3, summary = T);radial_ivw

beta <- c(beta, round(-res$b,3))
beta.CI <- c(beta.CI, paste0(round(-or.res$lo_ci,3), ' to ', round(-or.res$up_ci,3)))
OR <- c(OR, round(exp(-res$b), 3))
OR.CI <- c(OR.CI, paste0(round(exp(-or.res$up_ci),3), ' to ', round(exp(-or.res$lo_ci),3)))
p.value <- c(p.value, res$pval)

### heterogeneity test
heter <- mr_heterogeneity(LADA.CXCL10, method_list = 'mr_ivw'); heter

beta <- c(beta, round(heter$Q, 3))
beta.CI <- c(beta.CI, '')
OR <- c(OR, '')
OR.CI <- c(OR.CI, '')
p.value <- c(p.value, heter$Q_pval)

### I^2
Q <- heter$Q
Q_df <- heter$Q_df
I.sqr <- max(0, (Q-Q_df)/Q*100);I.sqr

beta <- c(beta, paste0(round(I.sqr, 2), '%'))
beta.CI <- c(beta.CI, '')
OR <- c(OR, '')
OR.CI <- c(OR.CI, '')
p.value <- c(p.value, '')

### mr-egger intercept ###

intercept <- mr_pleiotropy_test(LADA.CXCL10); intercept

beta <- c(beta, intercept$egger_intercept)
beta.CI <- c(beta.CI, '')
OR <- c(OR, '')
OR.CI <- c(OR.CI, '')
p.value <- c(p.value, intercept$pval)

exposure <- c(exposure, rep('LADA', 7))
outcome <- c(outcome, rep('CXCL10', 7))
method <- c(method, mr_method)


result <- data.frame(exposure=exposure, outcome=outcome, method=method, beta=beta, CI.95=beta.CI, 
                     OR=OR, CI.95=OR.CI, p.value=p.value)
write.csv(result, 'result_LADA_proteins.csv', quote = F, row.names = F)
