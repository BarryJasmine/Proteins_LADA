library(MendelianRandomization)
library(MVMR)
library(TwoSampleMR)
library(readr)
library(ieugwasr)

### CXCL10: chr4: 76021118-76023497
### read CXCL10 data ####
CXCL10 <- read.table(gzfile('E:/pQTL/LADA/4141_79_CXCL10_IP_10.txt.gz'), sep='\t', header=T, stringsAsFactors = F)
head(CXCL10)

### CXCL9:chr4: 76001275-76007509
### read CXCL9 data ###
CXCL9 <- read.table(gzfile('E:/pQTL/LADA/9188_119_CXCL9_MIG.txt.gz'), sep='\t', header=T, stringsAsFactors = F)
head(CXCL9)

### CXCL11:chr4: 76033682-76036070
### read CXCL9 data ###
CXCL11 <- read.table(gzfile('E:/pQTL/LADA/3038_9_CXCL11_I_TAC.txt.gz'), sep='\t', header=T, stringsAsFactors = F)
head(CXCL11)

### read LADA summary data
LADA.data <- read_tsv("D:/Others/药物和LADA/data/LADACTRL_MA_filtered_2018_DLC.txt", col_types=cols(.default="c"))
LADA.data[,5:10] <- sapply(LADA.data[,5:10], as.numeric)
LADA.data$beta <- log(LADA.data$OR)

LTL <- read_tsv(gzfile("D:/Others/药物和LADA/data/UKB_telomere_gwas_summarystats.tsv.gz"), col_types=cols(.default="c"))
LTL <- LTL[,c('variant_id', 'chromosome', 'base_pair_location')]
LTL$rs_number <- paste0(LTL$chromosome, ':', LTL$base_pair_location)
LTL <- LTL[,c('rs_number', 'variant_id')]

LADA.data <- merge(LADA.data, LTL, by='rs_number', all.x = T)

### cis-pqtl of CXCL10
CXCL10.ivs <- CXCL10[CXCL10$Pval<5e-8, ]
# upper <- 76023497
# down <- 76021118
# region <- 1000000
# CXCL10.ivs <- CXCL10.ivs[CXCL10.ivs$Chrom=='chr4',]
# CXCL10.ivs <- CXCL10.ivs[((CXCL10.ivs$Pos)<=upper+region)&((CXCL10.ivs$Pos)>=down-region), ]
dat <- data.frame(rsid=CXCL10.ivs$rsids, pval=CXCL10.ivs$Pval)
remained_snps <- ld_clump(dat, clump_kb = 10000, clump_r2 = 0.01, bfile = 'D:/LD_reference/EUR',
                          plink_bin = 'D:/LD_reference/plink/plink.exe')
CXCL10.ivs <- CXCL10.ivs[CXCL10.ivs$rsids%in%remained_snps$rsid, ]

### cis-pqtl of CXCL9 
CXCL9.ivs <- CXCL9[CXCL9$Pval<5e-8, ]
# upper <- 76007509
# down <- 76001275
# CXCL9.ivs <- CXCL9.ivs[CXCL9.ivs$Chrom=='chr4',]
# CXCL9.ivs <- CXCL9.ivs[((CXCL9.ivs$Pos)<=upper+region)&((CXCL9.ivs$Pos)>=down-region), ]
dat <- data.frame(rsid=CXCL9.ivs$rsids, pval=CXCL9.ivs$Pval)
remained_snps <- ld_clump(dat, clump_kb = 10000, clump_r2 = 0.001, bfile = 'D:/LD_reference/EUR',
                          plink_bin = 'D:/LD_reference/plink/plink.exe')
CXCL9.ivs <- CXCL9.ivs[CXCL9.ivs$rsids%in%remained_snps$rsid, ]

### cis-pqtl of CXCL11 
CXCL11.ivs <- CXCL11[CXCL11$Pval<5e-8, ]
# upper <- 76036070
# down <- 76033682
# CXCL11.ivs <- CXCL11.ivs[CXCL11.ivs$Chrom=='chr4',]
# CXCL11.ivs <- CXCL11.ivs[((CXCL11.ivs$Pos)<=upper+region)&((CXCL11.ivs$Pos)>=down-region), ]
dat <- data.frame(rsid=CXCL11.ivs$rsids, pval=CXCL11.ivs$Pval)
remained_snps <- ld_clump(dat, clump_kb = 10000, clump_r2 = 0.001, bfile = 'D:/LD_reference/EUR',
                          plink_bin = 'D:/LD_reference/plink/plink.exe')
CXCL11.ivs <- CXCL11.ivs[CXCL11.ivs$rsids%in%remained_snps$rsid, ]


### CXCL11, CXCL10 and CXCL9 significant IVs
ivs <- union(CXCL10.ivs$rsids, CXCL9.ivs$rsids)
ivs <- union(ivs, CXCL11.ivs$rsids)
### clump MVMR IVs
CXCL10.data <- CXCL10[CXCL10$rsids%in%ivs, ]
dat <- data.frame(rsid=CXCL10.data$rsids, pval=CXCL10.data$Pval)
retained_snps <- ld_clump(dat, clump_kb = 10000, clump_r2 = 0.001,
                            bfile='D:/LD_reference/EUR', plink_bin = 'D:/LD_reference/plink/plink.exe')
ivs <- retained_snps$rsid

CXCL10_CXCL9_CXCL11_to_LADA.data <- LADA.data[LADA.data$variant_id%in%ivs, ]
CXCL10.data <- CXCL10[CXCL10$rsids%in%ivs, ]
CXCL9.data <- CXCL9[CXCL9$rsids%in%ivs, ]
CXCL11.data <- CXCL11[CXCL11$rsids%in%ivs, ]
CXCL10.data <- format_data(CXCL10.data, snp_col = 'rsids', effect_allele_col = 'effectAllele',
                         other_allele_col = 'otherAllele', beta_col = 'Beta',
                         se_col = 'SE', pval_col = 'Pval')
CXCL9.data <- format_data(CXCL9.data, snp_col = 'rsids', effect_allele_col = 'effectAllele',
                        other_allele_col = 'otherAllele', beta_col = 'Beta',
                        se_col = 'SE', pval_col = 'Pval')
CXCL11.data <- format_data(CXCL11.data, snp_col = 'rsids', effect_allele_col = 'effectAllele',
                          other_allele_col = 'otherAllele', beta_col = 'Beta',
                          se_col = 'SE', pval_col = 'Pval')
CXCL10_CXCL9_CXCL11_to_LADA.data <- format_data(CXCL10_CXCL9_CXCL11_to_LADA.data, type='outcome', snp_col = 'variant_id', 
                                      effect_allele_col = 'reference_allele',
                                      other_allele_col = 'other_allele', eaf_col = 'eaf',
                                      beta_col = 'beta', se_col = 'OR_se', pval_col = 'P')

### change the effect allele, other allele, and beta values of CXCL9, CXCL11 and LADA according to the CXCL10 data
CXCL10.LADA <- harmonise_data(CXCL10.data, CXCL10_CXCL9_CXCL11_to_LADA.data)
CXCL10.LADA <- CXCL10.LADA[CXCL10.LADA$mr_keep==T, ]

CXCL9.data <- CXCL9.data[CXCL9.data$SNP%in%CXCL10.LADA$SNP, ]
for(rsid in CXCL9.data$SNP){
  CXCL9.EA <- CXCL9.data[CXCL9.data$SNP==rsid, 'effect_allele.exposure']
  CXCL10.EA <- CXCL10.LADA[CXCL10.LADA$SNP==rsid, 'effect_allele.exposure']
  if(CXCL9.EA == CXCL10.EA){
    next
  }else{
    CXCL9.data[CXCL9.data$SNP==rsid, 'effect_allele.exposure'] <- CXCL10.EA
    CXCL9.data[CXCL9.data$SNP==rsid, 'other_allele.exposure'] <- CXCL9.EA
    CXCL9.data[CXCL9.data$SNP==rsid, 'beta.exposure'] <- -1 * CXCL9.data[CXCL9.data$SNP==rsid, 'beta.exposure']
  }
}

CXCL11.data <- CXCL11.data[CXCL11.data$SNP%in%CXCL10.LADA$SNP, ]
for(rsid in CXCL11.data$SNP){
  CXCL11.EA <- CXCL11.data[CXCL11.data$SNP==rsid, 'effect_allele.exposure']
  CXCL10.EA <- CXCL10.LADA[CXCL10.LADA$SNP==rsid, 'effect_allele.exposure']
  if(CXCL11.EA == CXCL10.EA){
    next
  }else{
    CXCL11.data[CXCL11.data$SNP==rsid, 'effect_allele.exposure'] <- CXCL10.EA
    CXCL11.data[CXCL11.data$SNP==rsid, 'other_allele.exposure'] <- CXCL11.EA
    CXCL11.data[CXCL11.data$SNP==rsid, 'beta.exposure'] <- -1 * CXCL11.data[CXCL11.data$SNP==rsid, 'beta.exposure']
  }
}

CXCL10.LADA <- data.frame(SNP=CXCL10.LADA$SNP, EA=CXCL10.LADA$effect_allele.exposure, OA=CXCL10.LADA$other_allele.exposure,
                          EAF=CXCL10.LADA$eaf.outcome,beta.CXCL10=CXCL10.LADA$beta.exposure, beta.LADA=CXCL10.LADA$beta.outcome,
                        se.CXCL10=CXCL10.LADA$se.exposure, se.LADA=CXCL10.LADA$se.outcome, pval.CXCL10=CXCL10.LADA$pval.exposure,
                        pval.LADA=CXCL10.LADA$pval.outcome)
CXCL9.data <- data.frame(SNP=CXCL9.data$SNP, beta.CXCL9=CXCL9.data$beta.exposure, se.CXCL9=CXCL9.data$se.exposure, pval.CXCL9=CXCL9.data$pval.exposure)
CXCL11.data <- data.frame(SNP=CXCL11.data$SNP, beta.CXCL11=CXCL11.data$beta.exposure, se.CXCL11=CXCL11.data$se.exposure, pval.CXCL11=CXCL11.data$pval.exposure)

CXCL10.CXCL9.LADA <- merge(CXCL10.LADA, CXCL9.data, by='SNP', all.x=T)
CXCL11.CXCL10.CXCL9.LADA <- merge(CXCL10.CXCL9.LADA, CXCL11.data, by='SNP', all.x=T)
write.csv(CXCL11.CXCL10.CXCL9.LADA, 'CXCL11.CXCL10.CXCL9_to_LADA_IVs.csv', quote = F, row.names = F)
CXCL10.CXCL9.LADA.object <- mr_mvinput(bx=matrix(c(CXCL11.CXCL10.CXCL9.LADA$beta.CXCL10, CXCL11.CXCL10.CXCL9.LADA$beta.CXCL9,
                                                   CXCL11.CXCL10.CXCL9.LADA$beta.CXCL11), ncol=3),
                                    bxse=matrix(c(CXCL11.CXCL10.CXCL9.LADA$se.CXCL10, CXCL11.CXCL10.CXCL9.LADA$se.CXCL9,
                                                  CXCL11.CXCL10.CXCL9.LADA$se.CXCL11), ncol=3),
                                    by=CXCL11.CXCL10.CXCL9.LADA$beta.LADA, byse = CXCL11.CXCL10.CXCL9.LADA$se.LADA,
                                    exposure=c('CXCL10', 'CXCL9', 'CXCL11'), outcome='LADA')
generate_upper_lower_limit <- function(beta, lower_ci, upper_ci){
  return(paste0(round(exp(beta), 2), ' [', round(exp(lower_ci), 2), ',', round(exp(upper_ci), 2), ']'))
}

mr_mvivw(CXCL10.CXCL9.LADA.object, robust = T, nx=35559)
mr_mvmedian(CXCL10.CXCL9.LADA.object)
mr_mvlasso(CXCL10.CXCL9.LADA.object)
mr_mvegger(CXCL10.CXCL9.LADA.object)
