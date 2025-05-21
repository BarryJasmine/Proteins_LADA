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

#### UKB-PPP Analysis ####

### CXCL10 chr4:76942271-76944650###

fp <- read.table(gzfile('D:/Others/药物和LADA/data/discovery_chr4_CXCL10_P02778_OID20697_v1_Inflammation.gz'), 
                 sep=' ', header = T, stringsAsFactors = F)
annotation <- read.table(gzfile('D:/Others/药物和LADA/data/olink_rsid_map_mac5_info03_b0_7_chr4_patched_v2.tsv.gz'), 
                         sep='\t', header = T, stringsAsFactors = F)
annotation <- annotation[,c('ID', 'rsid', 'POS19')]
fp <- merge(fp, annotation, by='ID', all.x = T)
fp$Pval <- 10^(-fp$LOG10P)
fp <- fp[fp$Pval<5e-8, ]
upper <- 76944650
down <- 76942271
region <- 1000000
fp.ivs <- fp[(fp$POS19<upper+region)&(fp$POS19>down-region), ]

protein <- format_data(fp.ivs, snp_col = 'rsid', effect_allele_col = 'ALLELE1',
                       other_allele_col = 'ALLELE0', beta_col = 'BETA',
                       se_col = 'SE', pval_col = 'Pval', eaf_col = 'A1FREQ',
                       samplesize_col = 'N')
LADA <- LADA.data[LADA.data$variant_id%in%protein$SNP, ]
LADA <- format_data(LADA, type='outcome', snp_col = 'variant_id', effect_allele_col = 'reference_allele',
                    other_allele_col = 'other_allele', eaf_col = 'eaf',
                    beta_col = 'beta', se_col = 'OR_se', pval_col = 'P',
                    samplesize_col = 'n_samples')
protein.LADA <- harmonise_data(protein, LADA)
dat <- data.frame(rsid=protein.LADA$SNP, pval=protein.LADA$pval.exposure)
retained_snps <- ld_clump(dat, clump_kb = 10000, clump_r2 = 0.001, bfile = 'D:/LD_reference/EUR',
                          plink_bin = 'D:/LD_reference/plink/plink.exe')
protein.LADA <- protein.LADA[protein.LADA$SNP%in%retained_snps$rsid, ]
steiger_filtering(protein.LADA)

res <- generate_odds_ratios(mr(protein.LADA));res


dat <- data.frame(rsid=protein$SNP, pval=protein$pval.exposure)
retained_snps <- ld_clump(dat, clump_kb = 10000, clump_r2 = 0.001, bfile = 'D:/LD_reference/EUR',
                          plink_bin = 'D:/LD_reference/plink/plink.exe')
