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

#### Steiger filtering ####
Protein <- NULL
SNP <- NULL
EA <- NULL
OA <- NULL
EAF <- NULL
Beta <- NULL
SE <- NULL
Pval <- NULL

### CXCL10 ###
filename <- '4141_79_CXCL10_IP_10.txt'
#### cis-acting pQTL ####
input_file <- paste0('D:/Others/中介-蛋白组-MM/data/pval_5e-08_cis-acting/', filename)

fp <- read.table(input_file, sep=' ', header = T, stringsAsFactors = F)
protein <- format_data(fp, snp_col = 'rsids', effect_allele_col = 'effectAllele',
                       other_allele_col = 'otherAllele', beta_col = 'Beta',
                       se_col = 'SE', pval_col = 'Pval', eaf_col = 'effectAlleleFreq',
                       samplesize_col = 'N')
LADA <- LADA.data[LADA.data$variant_id%in%protein$SNP, ]

LADA <- format_data(LADA, type='outcome', snp_col = 'variant_id', effect_allele_col = 'reference_allele',
                    other_allele_col = 'other_allele', eaf_col = 'eaf',
                    beta_col = 'beta', se_col = 'OR_se', pval_col = 'P',
                    samplesize_col = 'n_samples')
protein.LADA <- harmonise_data(protein, LADA)

steiger_filtering(protein.LADA)

### SAA1 ###

filename <- '15515_2_SAA1_SAA.txt'
#### cis- and trans-acting pQTL ####
input_file <- paste0('D:/Others/中介-蛋白组-MM/data/pval_5e-08_cis-acting/', filename)

fp <- read.table(input_file, sep=' ', header = T, stringsAsFactors = F)
protein <- format_data(fp, snp_col = 'rsids', effect_allele_col = 'effectAllele',
                       other_allele_col = 'otherAllele', beta_col = 'Beta',
                       se_col = 'SE', pval_col = 'Pval', eaf_col = 'effectAlleleFreq',
                       samplesize_col = 'N')
LADA <- LADA.data[LADA.data$variant_id%in%protein$SNP, ]

LADA <- format_data(LADA, type='outcome', snp_col = 'variant_id', effect_allele_col = 'reference_allele',
                    other_allele_col = 'other_allele', eaf_col = 'eaf',
                    beta_col = 'beta', se_col = 'OR_se', pval_col = 'P',
                    samplesize_col = 'n_samples')
protein.LADA <- harmonise_data(protein, LADA)

steiger_filtering(protein.LADA)

### SAA2 ###

filename <- '18832_65_SAA2_SAA2.txt'
#### cis- and trans-acting pQTL ####
input_file <- paste0('D:/Others/中介-蛋白组-MM/data/pval_5e-08_cis-acting/', filename)

fp <- read.table(input_file, sep=' ', header = T, stringsAsFactors = F)
protein <- format_data(fp, snp_col = 'rsids', effect_allele_col = 'effectAllele',
                       other_allele_col = 'otherAllele', beta_col = 'Beta',
                       se_col = 'SE', pval_col = 'Pval', eaf_col = 'effectAlleleFreq',
                       samplesize_col = 'N')
LADA <- LADA.data[LADA.data$variant_id%in%protein$SNP, ]

LADA <- format_data(LADA, type='outcome', snp_col = 'variant_id', effect_allele_col = 'reference_allele',
                    other_allele_col = 'other_allele', eaf_col = 'eaf',
                    beta_col = 'beta', se_col = 'OR_se', pval_col = 'P',
                    samplesize_col = 'n_samples')
protein.LADA <- harmonise_data(protein, LADA)

steiger_filtering(protein.LADA)


