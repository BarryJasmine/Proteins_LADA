library(TwoSampleMR)
library(stringr)
library(readr)
library(ieugwasr)
library(echoLD)
library(data.table)   # For efficient data handling
library(Matrix)       # For sparse matrix operations
library(reticulate)   # To interface with Python

##### read LADA dataset
LADA.data <- read_tsv("D:/Others/药物和LADA/data/LADACTRL_MA_filtered_2018_DLC.txt", col_types=cols(.default="c"))
LADA.data[,5:10] <- sapply(LADA.data[,5:10], as.numeric)
LADA.data$beta <- log(LADA.data$OR)

LTL <- read_tsv(gzfile("D:/Others/药物和LADA/data/UKB_telomere_gwas_summarystats.tsv.gz"), col_types=cols(.default="c"))
LTL <- LTL[,c('variant_id', 'chromosome', 'base_pair_location')]
LTL$rs_number <- paste0(LTL$chromosome, ':', LTL$base_pair_location)
LTL <- LTL[,c('rs_number', 'variant_id')]

LADA.data <- merge(LADA.data, LTL, by='rs_number', all.x = T)

#### CXCL10 ####

CXCL10 <- read.table('D:/Others/药物和LADA/data/GTEx_CXCL10/Cells_EBV-transformed_lymphocytes.lite.CXCL10.txt',
                     sep='\t', stringsAsFactors = F, header = T)
CXCL10$N <- 147
dat <- data.frame(rsid=CXCL10$SNP, pval=CXCL10$p)
retained_snps <- ld_clump(dat, clump_kb = 1000, clump_p = 0.1, 
                          bfile = 'D:/LD_reference/EUR', plink_bin = 'D:/LD_reference/plink/plink.exe')
CXCL10 <- CXCL10[CXCL10$SNP%in%retained_snps$rsid, ]
CXCL10 <- format_data(CXCL10, snp_col = 'SNP', other_allele_col = 'A2',
                      effect_allele_col = 'A1', beta_col = 'b', se_col = 'SE', pval_col = 'p',
                      chr_col = 'Chr', pos_col = 'BP', samplesize_col = 'N')
LADA <- LADA.data[LADA.data$variant_id%in%CXCL10$SNP, ]
LADA <- format_data(LADA, type='outcome', snp_col = 'variant_id', effect_allele_col = 'reference_allele',
                    other_allele_col = 'other_allele', eaf_col = 'eaf',
                    beta_col = 'beta', se_col = 'OR_se', pval_col = 'P', samplesize_col = 'n_samples')
CXCL10.LADA <- harmonise_data(CXCL10, LADA)
steiger_filtering(CXCL10.LADA)
content <- data.frame(SNP=CXCL10.LADA$SNP, EA=CXCL10.LADA$effect_allele.exposure, OA=CXCL10.LADA$other_allele.exposure,
                      EAF=CXCL10.LADA$eaf.outcome, beta.exposure=CXCL10.LADA$beta.exposure, se.exposure=CXCL10.LADA$se.exposure,
                      pval.exposure=CXCL10.LADA$pval.exposure, beta.outcome=CXCL10.LADA$beta.outcome, se.outcome=CXCL10.LADA$se.outcome,
                      pval.outcome=CXCL10.LADA$pval.outcome, F_stat=(CXCL10.LADA$beta.exposure/CXCL10.LADA$se.exposure)^2)
write.csv(content, 'CXCL10_to_LADA_EBV.csv', quote = F, row.names = F)

generate_odds_ratios(mr(CXCL10.LADA))

### LCL ###
CXCL10 <- read.table('D:/Others/药物和LADA/data/LCL_CXCL10.txt',
                     sep='\t', stringsAsFactors = F, header = T)
CXCL10$N <- 373
LADA <- LADA.data[LADA.data$variant_id%in%CXCL10$SNP, ]
dat <- data.frame(rsid=LADA$variant_id, pval=LADA$P)
retained_snps <- ld_clump(dat, clump_kb = 1000, clump_p = 0.001, 
                          bfile = 'D:/LD_reference/EUR', plink_bin = 'D:/LD_reference/plink/plink.exe')
CXCL10 <- CXCL10[CXCL10$SNP%in%c(retained_snps$rsid), ]
CXCL10 <- format_data(CXCL10, snp_col = 'SNP', other_allele_col = 'A1',
                      effect_allele_col = 'A2', beta_col = 'b', se_col = 'SE', pval_col = 'p',
                      chr_col = 'Chr', pos_col = 'BP', eaf_col = 'Freq', samplesize_col = 'N')
LADA <- format_data(LADA, type='outcome', snp_col = 'variant_id', effect_allele_col = 'reference_allele',
                    other_allele_col = 'other_allele', eaf_col = 'eaf',
                    beta_col = 'beta', se_col = 'OR_se', pval_col = 'P', samplesize_col = 'n_samples')

CXCL10.LADA <- harmonise_data(CXCL10, LADA)
steiger_filtering(CXCL10.LADA)
content <- data.frame(SNP=CXCL10.LADA$SNP, EA=CXCL10.LADA$effect_allele.exposure, OA=CXCL10.LADA$other_allele.exposure,
                      EAF=CXCL10.LADA$eaf.outcome, beta.exposure=CXCL10.LADA$beta.exposure, se.exposure=CXCL10.LADA$se.exposure,
                      pval.exposure=CXCL10.LADA$pval.exposure, beta.outcome=CXCL10.LADA$beta.outcome, se.outcome=CXCL10.LADA$se.outcome,
                      pval.outcome=CXCL10.LADA$pval.outcome, F_stat=(CXCL10.LADA$beta.exposure/CXCL10.LADA$se.exposure)^2)
write.csv(content, 'CXCL10_to_LADA_LCL.csv', quote = F, row.names = F)

generate_odds_ratios(mr(CXCL10.LADA))


#### SAA1 ####

#### Artery Tibial ###
SAA1 <- read.table('D:/Others/药物和LADA/data/GTEx_SAA1/Artery_Tibial.lite.SAA1.txt',
                   sep='\t', stringsAsFactors = F, header = T)
dat <- data.frame(rsid=SAA1$SNP, pval=SAA1$p)
retained_snps <- ld_clump(dat, bfile = 'D:/LD_reference/EUR', plink_bin = 'D:/LD_reference/plink/plink.exe')
SAA1 <- SAA1[SAA1$SNP%in%retained_snps$rsid, ]
SAA1$N <- 584
SAA1 <- format_data(SAA1, snp_col = 'SNP', other_allele_col = 'A2',
                    effect_allele_col = 'A1', beta_col = 'b', se_col = 'SE', pval_col = 'p',
                    chr_col = 'Chr', pos_col = 'BP',samplesize_col = 'N')

LADA <- LADA.data[LADA.data$variant_id%in%SAA1$SNP, ]
LADA <- format_data(LADA, type='outcome', snp_col = 'variant_id', effect_allele_col = 'reference_allele',
                    other_allele_col = 'other_allele', eaf_col = 'eaf',
                    beta_col = 'beta', se_col = 'OR_se', pval_col = 'P', samplesize_col = 'n_samples')
SAA1.LADA <- harmonise_data(SAA1, LADA)
steiger_filtering(SAA1.LADA)
content <- data.frame(SNP=SAA1.LADA$SNP, EA=SAA1.LADA$effect_allele.exposure, OA=SAA1.LADA$other_allele.exposure,
                      EAF=SAA1.LADA$eaf.outcome, beta.exposure=SAA1.LADA$beta.exposure, se.exposure=SAA1.LADA$se.exposure,
                      pval.exposure=SAA1.LADA$pval.exposure, beta.outcome=SAA1.LADA$beta.outcome, se.outcome=SAA1.LADA$se.outcome,
                      pval.outcome=SAA1.LADA$pval.outcome, F_stat=(SAA1.LADA$beta.exposure/SAA1.LADA$se.exposure)^2)
write.csv(content, 'SAA1_to_LADA_artery.csv', quote = F, row.names = F)

generate_odds_ratios(mr(SAA1.LADA))

#### Lung ###
SAA1 <- read.table('D:/Others/药物和LADA/data/GTEx_SAA1/Lung.lite.SAA1.txt',
                   sep='\t', stringsAsFactors = F, header = T)
dat <- data.frame(rsid=SAA1$SNP, pval=SAA1$p)
retained_snps <- ld_clump(dat, clump_kb = 1000, clump_r2 = 0.01,
                          bfile = 'D:/LD_reference/EUR', plink_bin = 'D:/LD_reference/plink/plink.exe')
SAA1 <- SAA1[SAA1$SNP%in%retained_snps$rsid, ]
SAA1$N <- 515
SAA1 <- format_data(SAA1, snp_col = 'SNP', other_allele_col = 'A2',
                    effect_allele_col = 'A1', beta_col = 'b', se_col = 'SE', pval_col = 'p',
                    chr_col = 'Chr', pos_col = 'BP', samplesize_col = 'N')
LADA <- LADA.data[LADA.data$variant_id%in%SAA1$SNP, ]
LADA <- format_data(LADA, type='outcome', snp_col = 'variant_id', effect_allele_col = 'reference_allele',
                    other_allele_col = 'other_allele', eaf_col = 'eaf',
                    beta_col = 'beta', se_col = 'OR_se', pval_col = 'P', samplesize_col = 'n_samples')
SAA1.LADA <- harmonise_data(SAA1, LADA)
steiger_filtering(SAA1.LADA)
content <- data.frame(SNP=SAA1.LADA$SNP, EA=SAA1.LADA$effect_allele.exposure, OA=SAA1.LADA$other_allele.exposure,
                      EAF=SAA1.LADA$eaf.outcome, beta.exposure=SAA1.LADA$beta.exposure, se.exposure=SAA1.LADA$se.exposure,
                      pval.exposure=SAA1.LADA$pval.exposure, beta.outcome=SAA1.LADA$beta.outcome, se.outcome=SAA1.LADA$se.outcome,
                      pval.outcome=SAA1.LADA$pval.outcome, F_stat=(SAA1.LADA$beta.exposure/SAA1.LADA$se.exposure)^2)
write.csv(content, 'SAA1_to_LADA_lung.csv', quote = F, row.names = F)

generate_odds_ratios(mr(SAA1.LADA))

#### Nerve Tibial ###
SAA1 <- read.table('D:/Others/药物和LADA/data/GTEx_SAA1/Nerve_Tibial.lite.SAA1.txt',
                   sep='\t', stringsAsFactors = F, header = T)
dat <- data.frame(rsid=SAA1$SNP, pval=SAA1$p)
retained_snps <- ld_clump(dat, clump_kb = 1000, clump_r2 = 0.01,
                          bfile = 'D:/LD_reference/EUR', plink_bin = 'D:/LD_reference/plink/plink.exe')
SAA1 <- SAA1[SAA1$SNP%in%retained_snps$rsid, ]
SAA1$N <- 532
SAA1 <- format_data(SAA1, snp_col = 'SNP', other_allele_col = 'A2',
                    effect_allele_col = 'A1', beta_col = 'b', se_col = 'SE', pval_col = 'p',
                    chr_col = 'Chr', pos_col = 'BP', samplesize_col = 'N')
LADA <- LADA.data[LADA.data$variant_id%in%SAA1$SNP, ]
LADA <- format_data(LADA, type='outcome', snp_col = 'variant_id', effect_allele_col = 'reference_allele',
                    other_allele_col = 'other_allele', eaf_col = 'eaf',
                    beta_col = 'beta', se_col = 'OR_se', pval_col = 'P', samplesize_col = 'n_samples')
SAA1.LADA <- harmonise_data(SAA1, LADA)
SAA1.LADA$mr_keep <- T
steiger_filtering(SAA1.LADA)
content <- data.frame(SNP=SAA1.LADA$SNP, EA=SAA1.LADA$effect_allele.exposure, OA=SAA1.LADA$other_allele.exposure,
                      EAF=SAA1.LADA$eaf.outcome, beta.exposure=SAA1.LADA$beta.exposure, se.exposure=SAA1.LADA$se.exposure,
                      pval.exposure=SAA1.LADA$pval.exposure, beta.outcome=SAA1.LADA$beta.outcome, se.outcome=SAA1.LADA$se.outcome,
                      pval.outcome=SAA1.LADA$pval.outcome, F_stat=(SAA1.LADA$beta.exposure/SAA1.LADA$se.exposure)^2)
write.csv(content, 'SAA1_to_LADA_nerve_tibial.csv', quote = F, row.names = F)

generate_odds_ratios(mr(SAA1.LADA))

#### SAA2 ####

#### Adipose visceral omentum ###
SAA2 <- read.table('D:/Others/药物和LADA/data/GTEx_SAA2/Adipose_Visceral_Omentum.lite.SAA2.txt',
                   sep='\t', stringsAsFactors = F, header = T)
dat <- data.frame(rsid=SAA2$SNP, pval=SAA2$p)
retained_snps <- ld_clump(dat, bfile = 'D:/LD_reference/EUR', plink_bin = 'D:/LD_reference/plink/plink.exe')
SAA2 <- SAA2[SAA2$SNP%in%retained_snps$rsid, ]
SAA2$N <- 469
SAA2 <- format_data(SAA2, snp_col = 'SNP', other_allele_col = 'A2',
                    effect_allele_col = 'A1', beta_col = 'b', se_col = 'SE', pval_col = 'p',
                    chr_col = 'Chr', pos_col = 'BP', samplesize_col = 'N')
LADA <- LADA.data[LADA.data$variant_id%in%SAA2$SNP, ]
LADA <- format_data(LADA, type='outcome', snp_col = 'variant_id', effect_allele_col = 'reference_allele',
                    other_allele_col = 'other_allele', eaf_col = 'eaf',
                    beta_col = 'beta', se_col = 'OR_se', pval_col = 'P', samplesize_col = 'n_samples')
SAA2.LADA <- harmonise_data(SAA2, LADA)
steiger_filtering(SAA2.LADA)
content <- data.frame(SNP=SAA2.LADA$SNP, EA=SAA2.LADA$effect_allele.exposure, OA=SAA2.LADA$other_allele.exposure,
                      EAF=SAA2.LADA$eaf.outcome, beta.exposure=SAA2.LADA$beta.exposure, se.exposure=SAA2.LADA$se.exposure,
                      pval.exposure=SAA2.LADA$pval.exposure, beta.outcome=SAA2.LADA$beta.outcome, se.outcome=SAA2.LADA$se.outcome,
                      pval.outcome=SAA2.LADA$pval.outcome, F_stat=(SAA2.LADA$beta.exposure/SAA2.LADA$se.exposure)^2)
write.csv(content, 'SAA2_to_LADA_adipose.csv', quote = F, row.names = F)
generate_odds_ratios(mr(SAA2.LADA))

#### Artery_Tibial ###
SAA2 <- read.table('D:/Others/药物和LADA/data/GTEx_SAA2/Artery_Tibial.lite.SAA2.txt',
                   sep='\t', stringsAsFactors = F, header = T)
dat <- data.frame(rsid=SAA2$SNP, pval=SAA2$p)
retained_snps <- ld_clump(dat, clump_kb = 1000, clump_r2 = 0.01,
                          bfile = 'D:/LD_reference/EUR', plink_bin = 'D:/LD_reference/plink/plink.exe')
SAA2 <- SAA2[SAA2$SNP%in%retained_snps$rsid, ]
SAA2$N <- 584
SAA2 <- format_data(SAA2, snp_col = 'SNP', other_allele_col = 'A2',
                    effect_allele_col = 'A1', beta_col = 'b', se_col = 'SE', pval_col = 'p',
                    chr_col = 'Chr', pos_col = 'BP', samplesize_col = 'N')
LADA <- LADA.data[LADA.data$variant_id%in%SAA2$SNP, ]
LADA <- format_data(LADA, type='outcome', snp_col = 'variant_id', effect_allele_col = 'reference_allele',
                    other_allele_col = 'other_allele', eaf_col = 'eaf',
                    beta_col = 'beta', se_col = 'OR_se', pval_col = 'P', samplesize_col = 'n_samples')
SAA2.LADA <- harmonise_data(SAA2, LADA)
steiger_filtering(SAA2.LADA)
content <- data.frame(SNP=SAA2.LADA$SNP, EA=SAA2.LADA$effect_allele.exposure, OA=SAA2.LADA$other_allele.exposure,
                      EAF=SAA2.LADA$eaf.outcome, beta.exposure=SAA2.LADA$beta.exposure, se.exposure=SAA2.LADA$se.exposure,
                      pval.exposure=SAA2.LADA$pval.exposure, beta.outcome=SAA2.LADA$beta.outcome, se.outcome=SAA2.LADA$se.outcome,
                      pval.outcome=SAA2.LADA$pval.outcome, F_stat=(SAA2.LADA$beta.exposure/SAA2.LADA$se.exposure)^2)
write.csv(content, 'SAA2_to_LADA_artery.csv', quote = F, row.names = F)
generate_odds_ratios(mr(SAA2.LADA))

#### Nerve Tibial ###
SAA2 <- read.table('D:/Others/药物和LADA/data/GTEx_SAA2/Nerve_Tibial.lite.SAA2.txt',
                   sep='\t', stringsAsFactors = F, header = T)
dat <- data.frame(rsid=SAA2$SNP, pval=SAA2$p)
retained_snps <- ld_clump(dat, clump_kb = 1000, clump_r2 = 0.01,
                          bfile = 'D:/LD_reference/EUR', plink_bin = 'D:/LD_reference/plink/plink.exe')
SAA2 <- SAA2[SAA2$SNP%in%retained_snps$rsid, ]
SAA2$N <- 532
SAA2 <- format_data(SAA2, snp_col = 'SNP', other_allele_col = 'A2',
                    effect_allele_col = 'A1', beta_col = 'b', se_col = 'SE', pval_col = 'p',
                    chr_col = 'Chr', pos_col = 'BP', samplesize_col = 'N')
LADA <- LADA.data[LADA.data$variant_id%in%SAA2$SNP, ]
LADA <- format_data(LADA, type='outcome', snp_col = 'variant_id', effect_allele_col = 'reference_allele',
                    other_allele_col = 'other_allele', eaf_col = 'eaf',
                    beta_col = 'beta', se_col = 'OR_se', pval_col = 'P', samplesize_col = 'n_samples')
SAA2.LADA <- harmonise_data(SAA2, LADA)
SAA2.LADA$mr_keep <- T
steiger_filtering(SAA2.LADA)
content <- data.frame(SNP=SAA2.LADA$SNP, EA=SAA2.LADA$effect_allele.exposure, OA=SAA2.LADA$other_allele.exposure,
                      EAF=SAA2.LADA$eaf.outcome, beta.exposure=SAA2.LADA$beta.exposure, se.exposure=SAA2.LADA$se.exposure,
                      pval.exposure=SAA2.LADA$pval.exposure, beta.outcome=SAA2.LADA$beta.outcome, se.outcome=SAA2.LADA$se.outcome,
                      pval.outcome=SAA2.LADA$pval.outcome, F_stat=(SAA2.LADA$beta.exposure/SAA2.LADA$se.exposure)^2)
write.csv(content, 'SAA2_to_LADA_nerve.csv', quote = F, row.names = F)
generate_odds_ratios(mr(SAA2.LADA))
