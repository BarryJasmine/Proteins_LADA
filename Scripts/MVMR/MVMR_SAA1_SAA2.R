library(MendelianRandomization)
library(MVMR)
library(TwoSampleMR)
library(readr)
library(ieugwasr)

### SAA1: chr11: 18266264-18269967
### read SAA1 data ####
SAA1 <- read.table(gzfile('E:/pQTL/LADA/15515_2_SAA1_SAA.txt.gz'), sep='\t', header=T, stringsAsFactors = F)
head(SAA1)

### SAA2:chr11: 18238236-18248668
### read SAA2 data ###
SAA2 <- read.table(gzfile('E:/pQTL/LADA/18832_65_SAA2_SAA2.txt.gz'), sep='\t', header=T, stringsAsFactors = F)
head(SAA2)

### read LADA summary data
LADA.data <- read_tsv("D:/Others/药物和LADA/data/LADACTRL_MA_filtered_2018_DLC.txt", col_types=cols(.default="c"))
LADA.data[,5:10] <- sapply(LADA.data[,5:10], as.numeric)
LADA.data$beta <- log(LADA.data$OR)

LTL <- read_tsv(gzfile("D:/Others/药物和LADA/data/UKB_telomere_gwas_summarystats.tsv.gz"), col_types=cols(.default="c"))
LTL <- LTL[,c('variant_id', 'chromosome', 'base_pair_location')]
LTL$rs_number <- paste0(LTL$chromosome, ':', LTL$base_pair_location)
LTL <- LTL[,c('rs_number', 'variant_id')]

LADA.data <- merge(LADA.data, LTL, by='rs_number', all.x = T)

### cis-pqtl of SAA1
SAA1.ivs <- SAA1[SAA1$Pval<5e-8, ]
upper <- 18269967
down <- 18266264
region <- 1000000
SAA1.ivs <- SAA1.ivs[SAA1.ivs$Chrom=='chr11',]
SAA1.ivs <- SAA1.ivs[((SAA1.ivs$Pos)<=upper+region)&((SAA1.ivs$Pos)>=down-region), ]
dat <- data.frame(rsid=SAA1.ivs$rsids, pval=SAA1.ivs$Pval)
remained_snps <- ld_clump(dat, clump_kb = 10000, clump_r2 = 0.1, bfile = 'D:/LD_reference/EUR',
                          plink_bin = 'D:/LD_reference/plink/plink.exe')
SAA1.ivs <- SAA1.ivs[SAA1.ivs$rsids%in%remained_snps$rsid, ]

### cis-pqtl of SAA2 
SAA2.ivs <- SAA2[SAA2$Pval<5e-8, ]
upper <- 18248668
down <- 18238236
region <- 1000000
SAA2.ivs <- SAA2.ivs[SAA2.ivs$Chrom=='chr11',]
SAA2.ivs <- SAA2.ivs[((SAA2.ivs$Pos)<=upper+region)&((SAA2.ivs$Pos)>=down-region), ]
dat <- data.frame(rsid=SAA2.ivs$rsids, pval=SAA2.ivs$Pval)
remained_snps <- ld_clump(dat, clump_kb = 10000, clump_r2 = 0.1, bfile = 'D:/LD_reference/EUR',
                          plink_bin = 'D:/LD_reference/plink/plink.exe')
SAA2.ivs <- SAA2.ivs[SAA2.ivs$rsids%in%remained_snps$rsid, ]

### both SAA1 and SAA2 significant IVs
ivs <- union(SAA1.ivs$rsids, SAA2.ivs$rsids)
ivs <- ivs[ivs%in%LADA.data$variant_id]
SAA1.data <- SAA1[(SAA1$rsids%in%ivs)&!is.na(SAA1$rsids), ]
dat <- data.frame(rsid=SAA1.data$rsids, pval=SAA1.data$Pval)
retained_SNPs <- ld_clump(dat, clump_kb = 10000, clump_r2 = 0.1, 
                          bfile='D:/LD_reference/EUR', plink_bin = 'D:/LD_reference/plink/plink.exe')
ivs <- retained_SNPs$rsid

SAA1_SAA2_to_LADA.data <- LADA.data[LADA.data$variant_id%in%ivs, ]
SAA1.data <- SAA1[SAA1$rsids%in%ivs, ]
SAA2.data <- SAA2[SAA2$rsids%in%ivs, ]
SAA1.data <- format_data(SAA1.data, snp_col = 'rsids', effect_allele_col = 'effectAllele',
                         other_allele_col = 'otherAllele', beta_col = 'Beta',
                         se_col = 'SE', pval_col = 'Pval')
SAA2.data <- format_data(SAA2.data, snp_col = 'rsids', effect_allele_col = 'effectAllele',
                        other_allele_col = 'otherAllele', beta_col = 'Beta',
                        se_col = 'SE', pval_col = 'Pval')
SAA1_SAA2_to_LADA.data <- format_data(SAA1_SAA2_to_LADA.data, type='outcome', snp_col = 'variant_id', 
                                      effect_allele_col = 'reference_allele',
                                      other_allele_col = 'other_allele', eaf_col = 'eaf',
                                      beta_col = 'beta', se_col = 'OR_se', pval_col = 'P')
### change the effect allele, other allele, and beta values of SAA2 and LADA according to the SAA1 data
SAA1.LADA <- harmonise_data(SAA1.data, SAA1_SAA2_to_LADA.data)
SAA1.LADA <- SAA1.LADA[SAA1.LADA$mr_keep==T, ]

SAA2.data <- SAA2.data[SAA2.data$SNP%in%SAA1.LADA$SNP, ]
for(rsid in SAA2.data$SNP){
  SAA2.EA <- SAA2.data[SAA2.data$SNP==rsid, 'effect_allele.exposure']
  SAA1.EA <- SAA1.LADA[SAA1.LADA$SNP==rsid, 'effect_allele.exposure']
  if(SAA2.EA == SAA1.EA){
    next
  }else{
    SAA2.data[SAA2.data$SNP==rsid, 'effect_allele.exposure'] <- SAA1.EA
    SAA2.data[SAA2.data$SNP==rsid, 'other_allele.exposure'] <- SAA2.EA
    SAA2.data[SAA2.data$SNP==rsid, 'beta.exposure'] <- -1 * SAA2.data[SAA2.data$SNP==rsid, 'beta.exposure']
  }
}
SAA1.LADA <- data.frame(SNP=SAA1.LADA$SNP, beta.SAA1=SAA1.LADA$beta.exposure, beta.LADA=SAA1.LADA$beta.outcome,
                        se.SAA1=SAA1.LADA$se.exposure, se.LADA=SAA1.LADA$se.outcome)
SAA2.data <- data.frame(SNP=SAA2.data$SNP, beta.SAA2=SAA2.data$beta.exposure, se.SAA2=SAA2.data$se.exposure)
SAA1.SAA2.LADA <- merge(SAA1.LADA, SAA2.data, by='SNP', all.x=T)
SAA1.SAA2.LADA.object <- mr_mvinput(bx=matrix(c(SAA1.SAA2.LADA$beta.SAA1, SAA1.SAA2.LADA$beta.SAA2), ncol=2),
                                    bxse=matrix(c(SAA1.SAA2.LADA$se.SAA1, SAA1.SAA2.LADA$se.SAA2), ncol=2),
                                    by=SAA1.SAA2.LADA$beta.LADA, byse = SAA1.SAA2.LADA$se.LADA,
                                    exposure=c('SAA1', 'SAA2'), outcome='LADA')
strength_mvmr(SAA1.SAA2.LADA.object)

mr_mvivw(SAA1.SAA2.LADA.object, nx=35559)
mr_mvegger(SAA1.SAA2.LADA.object)
mr_mvlasso(SAA1.SAA2.LADA.object)
mr_mvcML(SAA1.SAA2.LADA.object, n=8581, rho_mat = matrix(c(1,0.9,0.05,0.9,1,0.04,0.05,0.04,1), nrow=3,ncol=3))
