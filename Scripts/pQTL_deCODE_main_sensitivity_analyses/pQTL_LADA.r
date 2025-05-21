library(dplyr)
library(readr)
library(TwoSampleMR)
library(tidyverse)
library(ieugwasr)
library(MRPRESSO)
library(stringr)

LADA.data <- read_tsv("D:/Others/药物和LADA/data/LADACTRL_MA_filtered_2018_DLC.txt", col_types=cols(.default="c"))
LADA.data[,5:10] <- sapply(LADA.data[,5:10], as.numeric)
LADA.data$beta <- log(LADA.data$OR)

LTL <- read_tsv(gzfile("D:/Others/药物和LADA/data/UKB_telomere_gwas_summarystats.tsv.gz"), col_types=cols(.default="c"))
LTL <- LTL[,c('variant_id', 'chromosome', 'base_pair_location')]
LTL$rs_number <- paste0(LTL$chromosome, ':', LTL$base_pair_location)
LTL <- LTL[,c('rs_number', 'variant_id')]

LADA.data <- merge(LADA.data, LTL, by='rs_number', all.x = T)

#### calculate F-statistics ####

file <- read.csv('pQTLs_num.csv', sep=',', header = T)
colnames(file)[1] <- 'filename'
filenames <- file$filename[file$MHC=='No']

Decode_genes <- read.table('D:/Others/中介-蛋白组-MM/data/Decode_genes_with_position.tsv',
                           sep='\t', header = T, stringsAsFactors = F)

files <- NULL
name <- NULL
num_SNPs <- NULL
SNPs <- NULL
EA <- NULL
OA <- NULL
EAF <- NULL
beta.exposure <- NULL
se.exposure <- NULL
p.exposure <- NULL
beta.outcome <- NULL
se.outcome <- NULL
p.outcome <- NULL
F.stat <- NULL
cnt <- 0

for (filename in filenames){
  cnt <- cnt + 1
  #### cis-acting pQTL ####
  input_file <- paste0('D:/Others/中介-蛋白组-MM/data/pval_5e-08_cis-acting/', filename)
  seqid <- paste0(unlist(str_split(filename, '_'))[1], '_', unlist(str_split(filename, '_'))[2])
  gene_name <- Decode_genes$Gene[Decode_genes$SeqId==seqid]
  fp <- read.table(input_file, sep=' ', header = T, stringsAsFactors = F)
  protein <- format_data(fp, snp_col = 'rsids', effect_allele_col = 'effectAllele',
                         other_allele_col = 'otherAllele', beta_col = 'Beta',
                         se_col = 'SE', pval_col = 'Pval', eaf_col = 'effectAlleleFreq')
  LADA <- LADA.data[LADA.data$variant_id%in%protein$SNP, ]
  LADA <- format_data(LADA, type='outcome', snp_col = 'variant_id', effect_allele_col = 'reference_allele',
                      other_allele_col = 'other_allele', eaf_col = 'eaf',
                      beta_col = 'beta', se_col = 'OR_se', pval_col = 'P')
  protein.LADA <- harmonise_data(protein, LADA)
  protein.LADA$F.stat <- (protein.LADA$beta.exposure/protein.LADA$se.exposure)^2
  protein.LADA <- protein.LADA[protein.LADA$mr_keep==T, ]
  if (nrow(protein.LADA[protein.LADA$mr_keep==T, ])==0){
    cat('harmonize IVs个数不够\n')
  }else if (nrow(protein.LADA[protein.LADA$mr_keep==T, ])==1){
    res <- generate_odds_ratios(mr(protein.LADA, method_list = 'mr_wald_ratio'))
    files <- c(files, rep(filename, res$nsnp))
    name <- c(name, rep(gene_name, res$nsnp))
    num_SNPs <- c(num_SNPs, rep(res$nsnp, res$nsnp))
    SNPs <- c(SNPs, protein.LADA$SNP)
    EA <- c(EA, protein.LADA$effect_allele.exposure)
    OA <- c(OA, protein.LADA$other_allele.exposure)
    EAF <- c(EAF, protein.LADA$eaf.exposure)
    beta.exposure <- c(beta.exposure, protein.LADA$beta.exposure)
    se.exposure <- c(se.exposure, protein.LADA$se.exposure)
    p.exposure <- c(p.exposure, protein.LADA$pval.exposure)
    beta.outcome <- c(beta.outcome, protein.LADA$beta.outcome)
    se.outcome <- c(se.outcome, protein.LADA$se.outcome)
    p.outcome <- c(p.outcome, protein.LADA$pval.outcome)
    F.stat <- c(F.stat, protein.LADA$F.stat)
  }else{
    res <- generate_odds_ratios(mr(protein.LADA, method_list = 'mr_ivw'))
    files <- c(files, rep(filename, res$nsnp))
    name <- c(name, rep(gene_name, res$nsnp))
    num_SNPs <- c(num_SNPs, rep(res$nsnp, res$nsnp))
    SNPs <- c(SNPs, protein.LADA$SNP)
    EA <- c(EA, protein.LADA$effect_allele.exposure)
    OA <- c(OA, protein.LADA$other_allele.exposure)
    EAF <- c(EAF, protein.LADA$eaf.exposure)
    beta.exposure <- c(beta.exposure, protein.LADA$beta.exposure)
    se.exposure <- c(se.exposure, protein.LADA$se.exposure)
    p.exposure <- c(p.exposure, protein.LADA$pval.exposure)
    beta.outcome <- c(beta.outcome, protein.LADA$beta.outcome)
    se.outcome <- c(se.outcome, protein.LADA$se.outcome)
    p.outcome <- c(p.outcome, protein.LADA$pval.outcome)
    F.stat <- c(F.stat, protein.LADA$F.stat)
  }
  if (cnt %% 200 == 0){
    cat(paste0('\n----------------', cnt,'/1389-----------------\n'))
    results <- data.frame(filename=files, name=name, num_SNPs=num_SNPs, SNPs=SNPs, EA=EA, OA=OA, EAF=EAF,
                          beta.exposure=beta.exposure, se.exposure=se.exposure, p.exposure=p.exposure,
                          beta.outcome=beta.outcome, se.outcome=se.outcome, p.outcome=p.outcome, 
                          F.stat=F.stat)
    write.csv(results, 'pQTL_LADA_5e-8_all_SNPs_information.csv', quote = F, row.names = F)
  }
}

results <- data.frame(filename=files, name=name, SNPs=SNPs, EA=EA, OA=OA, EAF=EAF,
                      beta.exposure=beta.exposure, se.exposure=se.exposure, p.exposure=p.exposure,
                      beta.outcome=beta.outcome, se.outcome=se.outcome, p.outcome=p.outcome, 
                      F.stat=F.stat)
write.csv(results, 'pQTL_LADA_5e-8_all_SNPs_information.csv', quote = F, row.names = F)


#### Analysis ####
SNPs <- NULL
files <- NULL
name <- NULL
beta <- NULL
se <- NULL
OR <- NULL
CI_95 <- NULL
P <- NULL
beta_ivw_ <- NULL
se_ivw_ <- NULL
OR_ivw_ <- NULL
CI_95_ivw_ <- NULL
P_ivw_ <- NULL
cnt <- 0

Decode_genes <- read.table('D:/代做/中介-蛋白组-MM/data/Decode_genes_with_position.tsv',
                                                      sep='\t', header = T, stringsAsFactors = F)

#### cis-acting pQTL ####
filenames <- list.files('D:/代做/中介-蛋白组-MM/data/pval_5e-08_cis-acting/')

for (filename in filenames){
  cnt <- cnt + 1
  # #### cis-acting pQTL ####
  input_file <- paste0('D:/代做/中介-蛋白组-MM/data/pval_5e-08_cis-acting/', filename)
  seqid <- paste0(unlist(str_split(filename, '_'))[1], '_', unlist(str_split(filename, '_'))[2])
  gene_name <- Decode_genes$Gene[Decode_genes$SeqId==seqid]
  fp <- read.table(input_file, sep=' ', header = T, stringsAsFactors = F)
  protein <- format_data(fp, snp_col = 'rsids', effect_allele_col = 'effectAllele',
                         other_allele_col = 'otherAllele', beta_col = 'Beta',
                         se_col = 'SE', pval_col = 'Pval', eaf_col = 'effectAlleleFreq')
  LADA <- LADA.data[LADA.data$variant_id%in%protein$SNP, ]
  if (is.null(LADA)){
    cat('LADA IVs个数不够\n')
    next
  }
  if (nrow(LADA)<1){
    cat('LADA IVs个数不够\n')
    next
  }
  LADA <- format_data(LADA, type='outcome', snp_col = 'variant_id', effect_allele_col = 'reference_allele',
                      other_allele_col = 'other_allele', eaf_col = 'eaf',
                      beta_col = 'beta', se_col = 'OR_se', pval_col = 'P')
  protein.LADA <- harmonise_data(protein, LADA)
  if (nrow(protein.LADA[protein.LADA$mr_keep==T, ])==0){
    cat('harmonize IVs个数不够\n')
  }else if (nrow(protein.LADA[protein.LADA$mr_keep==T, ])==1){
    res <- generate_odds_ratios(mr(protein.LADA, method_list = 'mr_wald_ratio'))
    files <- c(files, filename)
    SNPs <- c(SNPs, res$nsnp)
    beta <- c(beta, round(res$b, 4))
    se <- c(se, round(res$se, 4))
    OR <- c(OR, round(res$or, 4))
    CI_95 <- c(CI_95, paste0(round(res$or_lci95, 4), '-', round(res$or_uci95, 4)))
    P <- c(P, res$pval)
    beta_ivw_ <- c(beta_ivw_, round(res$b, 4))
    se_ivw_ <- c(se_ivw_, round(res$se, 4))
    OR_ivw_ <- c(OR_ivw_, round(res$or, 4))
    CI_95_ivw_ <- c(CI_95_ivw_, paste0(round(res$or_lci95, 4), '-', round(res$or_uci95, 4)))
    P_ivw_ <- c(P_ivw_, res$pval)
  }else{
    res <- generate_odds_ratios(mr(protein.LADA, method_list = 'mr_ivw'))
    files <- c(files, filename)
    SNPs <- c(SNPs, res$nsnp)
    beta <- c(beta, round(res$b, 4))
    se <- c(se, round(res$se, 4))
    OR <- c(OR, round(res$or, 4))
    CI_95 <- c(CI_95, paste0(round(res$or_lci95, 4), '-', round(res$or_uci95, 4)))
    P <- c(P, res$pval)
    Q <- mr_heterogeneity(protein.LADA, method_list = 'mr_ivw')
    Q <- Q$Q_pval
    if (Q>0.05){
      res <- generate_odds_ratios(mr(protein.LADA, method_list = 'mr_ivw_fe'))
      beta_ivw_ <- c(beta_ivw_, round(res$b, 4))
      se_ivw_ <- c(se_ivw_, round(res$se, 4))
      OR_ivw_ <- c(OR_ivw_, round(res$or, 4))
      CI_95_ivw_ <- c(CI_95_ivw_, paste0(round(res$or_lci95, 4), '-', round(res$or_uci95, 4)))
      P_ivw_ <- c(P_ivw_, res$pval)
    }else{
      res <- generate_odds_ratios(mr(protein.LADA, method_list = 'mr_ivw'))
      beta_ivw_ <- c(beta_ivw_, round(res$b, 4))
      OR_ivw_ <- c(OR_ivw_, round(res$or, 4))
      se_ivw_ <- c(se_ivw_, round(res$se, 4))
      CI_95_ivw_ <- c(CI_95_ivw_, paste0(round(res$or_lci95, 4), '-', round(res$or_uci95, 4)))
      P_ivw_ <- c(P_ivw_, res$pval)
    }
  }
  if (cnt %% 200 == 0){
    cat(paste0('\n----------------', cnt,'/4920-----------------\n'))
    results <- data.frame(filename=files, SNPs=SNPs, Beta=beta, OR=OR, CI_95=CI_95, P_value=P,
                          Beta_ivw_=beta_ivw_, OR_ivw_=OR_ivw_, CI_95_ivw_=CI_95_ivw_, P_ivw_=P_ivw_)
    write.csv(results, 'pQTL_LADA_5e-8.csv', quote = F, row.names = F)
  }
}

results <- data.frame(filename=files, SNPs=SNPs, Beta=beta, OR=OR, SE=se, CI_95=CI_95, P_value=P,
                      Beta_ivw_=beta_ivw_, SE_ivw_=se_ivw_, OR_ivw_=OR_ivw_, CI_95_ivw_=CI_95_ivw_, P_ivw_=P_ivw_)
write.csv(results, 'pQTL_LADA_5e-8.csv', quote = F, row.names = F)

#### Multiple Locus Analysis ####

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

#### cis- and trans-acting pQTL ####
input_file <- paste0('D:/Others/中介-蛋白组-MM/data/pval_5e-08_clump_10000_0.001/', filename)

fp <- read.table(input_file, sep=' ', header = T, stringsAsFactors = F)
protein <- format_data(fp, snp_col = 'rsids', effect_allele_col = 'effectAllele',
                       other_allele_col = 'otherAllele', beta_col = 'Beta',
                       se_col = 'SE', pval_col = 'Pval', eaf_col = 'effectAlleleFreq', samplesize_col = 'N')
LADA <- LADA.data[LADA.data$variant_id%in%protein$SNP, ]

LADA <- format_data(LADA, type='outcome', snp_col = 'variant_id', effect_allele_col = 'reference_allele',
                    other_allele_col = 'other_allele', eaf_col = 'eaf',
                    beta_col = 'beta', se_col = 'OR_se', pval_col = 'P', samplesize_col = 'n_samples')
protein.LADA <- harmonise_data(protein, LADA)
protein.LADA <- protein.LADA[protein.LADA$mr_keep==T, ]
steiger_filtering(protein.LADA)
CXCL10.LADA.multi_locus <- data.frame(SNP=protein.LADA$SNP, EA=protein.LADA$effect_allele.exposure,
                                     OA=protein.LADA$other_allele.exposure, EAF=protein.LADA$eaf.exposure,
                                     beta.exposure=protein.LADA$beta.exposure, se.exposure=protein.LADA$se.exposure,
                                     pval.exposure=protein.LADA$pval.exposure, beta.outcome=protein.LADA$beta.outcome,
                                     se.outcome=protein.LADA$se.outcome, pval.outcome=protein.LADA$pval.outcome)
write.csv(CXCL10.LADA.multi_locus, 'CXCL10.LADA.multi_locus_data.csv',
          row.names = F, quote = F)

res <- generate_odds_ratios(mr(protein.LADA, method_list = c('mr_ivw', 'mr_egger_regression', 'mr_weighted_median', 'mr_raps')));res

#### Radial MR ####
radial_data <- RadialMR::format_radial(BXG=protein.LADA$beta.exposure, BYG = protein.LADA$beta.outcome,
                             seBXG=protein.LADA$se.exposure, seBYG = protein.LADA$se.outcome,
                             RSID=protein.LADA$SNP)
radial_ivw <- RadialMR::ivw_radial(radial_data, alpha=0.05, weights = 3, summary = T);radial_ivw

Q <- mr_heterogeneity(protein.LADA, method_list = 'mr_ivw');Q
I_sqr <- paste0(min(round((Q$Q-Q$Q_df)/Q$Q, 4),1)*100, '%');I_sqr

mr_pleiotropy_test(protein.LADA)

Protein <- c(Protein, rep('CXCL10', sum(protein.LADA$mr_keep)))
SNP <- c(SNP, protein.LADA[protein.LADA$mr_keep==T, ]$SNP)
EA <- c(EA, protein.LADA[protein.LADA$mr_keep==T, ]$effect_allele.exposure)
OA <- c(OA, protein.LADA[protein.LADA$mr_keep==T, ]$other_allele.exposure)
EAF <- c(EAF, protein.LADA[protein.LADA$mr_keep==T, ]$eaf.exposure)
Beta <- c(Beta, protein.LADA[protein.LADA$mr_keep==T, ]$beta.exposure)
SE <- c(SE, protein.LADA[protein.LADA$mr_keep==T, ]$se.exposure)
Pval <- c(Pval, protein.LADA[protein.LADA$mr_keep==T, ]$pval.exposure)

### SAA1 ###

filename <- '15515_2_SAA1_SAA.txt'
#### cis- and trans-acting pQTL ####
input_file <- paste0('D:/Others/中介-蛋白组-MM/data/pval_5e-08_clump_10000_0.001/', filename)

fp <- read.table(input_file, sep=' ', header = T, stringsAsFactors = F)
protein <- format_data(fp, snp_col = 'rsids', effect_allele_col = 'effectAllele',
                       other_allele_col = 'otherAllele', beta_col = 'Beta',
                       se_col = 'SE', pval_col = 'Pval', eaf_col = 'effectAlleleFreq', samplesize_col = 'N')
LADA <- LADA.data[LADA.data$variant_id%in%protein$SNP, ]

LADA <- format_data(LADA, type='outcome', snp_col = 'variant_id', effect_allele_col = 'reference_allele',
                    other_allele_col = 'other_allele', eaf_col = 'eaf',
                    beta_col = 'beta', se_col = 'OR_se', pval_col = 'P', samplesize_col = 'n_samples')
protein.LADA <- harmonise_data(protein, LADA)
protein.LADA <- protein.LADA[protein.LADA$mr_keep==T, ]
steiger_filtering(protein.LADA)
SAA1.LADA.multi_locus <- data.frame(SNP=protein.LADA$SNP, EA=protein.LADA$effect_allele.exposure,
                                      OA=protein.LADA$other_allele.exposure, EAF=protein.LADA$eaf.exposure,
                                      beta.exposure=protein.LADA$beta.exposure, se.exposure=protein.LADA$se.exposure,
                                      pval.exposure=protein.LADA$pval.exposure, beta.outcome=protein.LADA$beta.outcome,
                                      se.outcome=protein.LADA$se.outcome, pval.outcome=protein.LADA$pval.outcome)
write.csv(SAA1.LADA.multi_locus, 'SAA1.LADA.multi_locus_data.csv',
          row.names = F, quote = F)

res <- generate_odds_ratios(mr(protein.LADA, method_list = c('mr_ivw', 'mr_egger_regression', 'mr_weighted_median', 'mr_raps', 'mr_ivw_radial')));res
Q <- mr_heterogeneity(protein.LADA, method_list = 'mr_ivw');Q
I_sqr <- paste0(min(round((Q$Q-Q$Q_df)/Q$Q, 4),1)*100, '%');I_sqr

mr_pleiotropy_test(protein.LADA)

Protein <- c(Protein, rep('SAA1', sum(protein.LADA$mr_keep)))
SNP <- c(SNP, protein.LADA[protein.LADA$mr_keep==T, ]$SNP)
EA <- c(EA, protein.LADA[protein.LADA$mr_keep==T, ]$effect_allele.exposure)
OA <- c(OA, protein.LADA[protein.LADA$mr_keep==T, ]$other_allele.exposure)
EAF <- c(EAF, protein.LADA[protein.LADA$mr_keep==T, ]$eaf.exposure)
Beta <- c(Beta, protein.LADA[protein.LADA$mr_keep==T, ]$beta.exposure)
SE <- c(SE, protein.LADA[protein.LADA$mr_keep==T, ]$se.exposure)
Pval <- c(Pval, protein.LADA[protein.LADA$mr_keep==T, ]$pval.exposure)

### SAA2 ###

filename <- '18832_65_SAA2_SAA2.txt'
#### cis- and trans-acting pQTL ####
input_file <- paste0('D:/Others/中介-蛋白组-MM/data/pval_5e-08_clump_10000_0.001/', filename)

fp <- read.table(input_file, sep=' ', header = T, stringsAsFactors = F)
protein <- format_data(fp, snp_col = 'rsids', effect_allele_col = 'effectAllele',
                       other_allele_col = 'otherAllele', beta_col = 'Beta',
                       se_col = 'SE', pval_col = 'Pval', eaf_col = 'effectAlleleFreq', samplesize_col = 'N')
LADA <- LADA.data[LADA.data$variant_id%in%protein$SNP, ]

LADA <- format_data(LADA, type='outcome', snp_col = 'variant_id', effect_allele_col = 'reference_allele',
                    other_allele_col = 'other_allele', eaf_col = 'eaf',
                    beta_col = 'beta', se_col = 'OR_se', pval_col = 'P', samplesize_col = 'n_samples')
protein.LADA <- harmonise_data(protein, LADA)
protein.LADA <- protein.LADA[protein.LADA$mr_keep==T, ]
steiger_filtering(protein.LADA)
SAA2.LADA.multi_locus <- data.frame(SNP=protein.LADA$SNP, EA=protein.LADA$effect_allele.exposure,
                                      OA=protein.LADA$other_allele.exposure, EAF=protein.LADA$eaf.exposure,
                                      beta.exposure=protein.LADA$beta.exposure, se.exposure=protein.LADA$se.exposure,
                                      pval.exposure=protein.LADA$pval.exposure, beta.outcome=protein.LADA$beta.outcome,
                                      se.outcome=protein.LADA$se.outcome, pval.outcome=protein.LADA$pval.outcome)
write.csv(SAA2.LADA.multi_locus, 'SAA2.LADA.multi_locus_data.csv',
          row.names = F, quote = F)


res <- generate_odds_ratios(mr(protein.LADA, method_list = c('mr_ivw', 'mr_egger_regression', 'mr_weighted_median', 'mr_raps', 'mr_ivw_radial')));res
Q <- mr_heterogeneity(protein.LADA, method_list = 'mr_ivw');Q
I_sqr <- paste0(max(round((Q$Q-Q$Q_df)/Q$Q, 4),0)*100, '%');I_sqr

mr_pleiotropy_test(protein.LADA)

Protein <- c(Protein, rep('SAA2', sum(protein.LADA$mr_keep)))
SNP <- c(SNP, protein.LADA[protein.LADA$mr_keep==T, ]$SNP)
EA <- c(EA, protein.LADA[protein.LADA$mr_keep==T, ]$effect_allele.exposure)
OA <- c(OA, protein.LADA[protein.LADA$mr_keep==T, ]$other_allele.exposure)
EAF <- c(EAF, protein.LADA[protein.LADA$mr_keep==T, ]$eaf.exposure)
Beta <- c(Beta, protein.LADA[protein.LADA$mr_keep==T, ]$beta.exposure)
SE <- c(SE, protein.LADA[protein.LADA$mr_keep==T, ]$se.exposure)
Pval <- c(Pval, protein.LADA[protein.LADA$mr_keep==T, ]$pval.exposure)

multi.locus.IVs <- data.frame(Protein=Protein, SNP=SNP, EA=EA, OA=OA,
                              EAF=EAF, Beta=Beta, SE=SE, Pval=Pval)
write.csv(multi.locus.IVs, 'Multi_locus_analysis_IVs.csv', quote = F,
          row.names = F)

#### Analysis u####
SNPs <- NULL
files <- NULL
name <- NULL
beta <- NULL
se <- NULL
OR <- NULL
CI_95 <- NULL
P <- NULL
beta_ivw_ <- NULL
se_ivw_ <- NULL
OR_ivw_ <- NULL
CI_95_ivw_ <- NULL
P_ivw_ <- NULL
cnt <- 0

Decode_genes <- read.table('D:/代做/中介-蛋白组-MM/data/Decode_genes_with_position.tsv',
                           sep='\t', header = T, stringsAsFactors = F)

#### cis-acting pQTL ####
filenames <- list.files('D:/代做/中介-蛋白组-MM/data/pval_5e-08_cis-acting/')

for (filename in filenames){
  cnt <- cnt + 1
  #### cis-acting pQTL ####
  input_file <- paste0('D:/代做/中介-蛋白组-MM/data/pval_5e-08_cis-acting/', filename)
  seqid <- paste0(unlist(str_split(filename, '_'))[1], '_', unlist(str_split(filename, '_'))[2])
  gene_name <- Decode_genes$Gene[Decode_genes$SeqId==seqid]
  fp <- read.table(input_file, sep=' ', header = T, stringsAsFactors = F)
  protein <- format_data(fp, snp_col = 'rsids', effect_allele_col = 'effectAllele',
                         other_allele_col = 'otherAllele', beta_col = 'Beta',
                         se_col = 'SE', pval_col = 'Pval', eaf_col = 'effectAlleleFreq')
  LADA <- LADA.data[LADA.data$variant_id%in%protein$SNP, ]
  if (is.null(LADA)){
    cat('LADA IVs个数不够\n')
    next
  }
  if (nrow(LADA)<1){
    cat('LADA IVs个数不够\n')
    next
  }
  LADA <- format_data(LADA, type='outcome', snp_col = 'variant_id', effect_allele_col = 'reference_allele',
                      other_allele_col = 'other_allele', eaf_col = 'eaf',
                      beta_col = 'beta', se_col = 'OR_se', pval_col = 'P')
  protein.LADA <- harmonise_data(protein, LADA)
  if (nrow(protein.LADA[protein.LADA$mr_keep==T, ])==0){
    cat('harmonize IVs个数不够\n')
  }else if (nrow(protein.LADA[protein.LADA$mr_keep==T, ])==1){
    res <- generate_odds_ratios(mr(protein.LADA, method_list = 'mr_wald_ratio'))
    files <- c(files, filename)
    SNPs <- c(SNPs, res$nsnp)
    beta <- c(beta, round(res$b, 4))
    se <- c(se, round(res$se, 4))
    OR <- c(OR, round(res$or, 4))
    CI_95 <- c(CI_95, paste0(round(res$or_lci95, 4), '-', round(res$or_uci95, 4)))
    P <- c(P, res$pval)
    beta_ivw_ <- c(beta_ivw_, round(res$b, 4))
    se_ivw_ <- c(se_ivw_, round(res$se, 4))
    OR_ivw_ <- c(OR_ivw_, round(res$or, 4))
    CI_95_ivw_ <- c(CI_95_ivw_, paste0(round(res$or_lci95, 4), '-', round(res$or_uci95, 4)))
    P_ivw_ <- c(P_ivw_, res$pval)
  }else{
    res <- generate_odds_ratios(mr(protein.LADA, method_list = 'mr_ivw'))
    files <- c(files, filename)
    SNPs <- c(SNPs, res$nsnp)
    beta <- c(beta, round(res$b, 4))
    se <- c(se, round(res$se, 4))
    OR <- c(OR, round(res$or, 4))
    CI_95 <- c(CI_95, paste0(round(res$or_lci95, 4), '-', round(res$or_uci95, 4)))
    P <- c(P, res$pval)
    Q <- mr_heterogeneity(protein.LADA, method_list = 'mr_ivw')
    Q <- Q$Q_pval
    if (Q>0.05){
      res <- generate_odds_ratios(mr(protein.LADA, method_list = 'mr_ivw_fe'))
      beta_ivw_ <- c(beta_ivw_, round(res$b, 4))
      se_ivw_ <- c(se_ivw_, round(res$se, 4))
      OR_ivw_ <- c(OR_ivw_, round(res$or, 4))
      CI_95_ivw_ <- c(CI_95_ivw_, paste0(round(res$or_lci95, 4), '-', round(res$or_uci95, 4)))
      P_ivw_ <- c(P_ivw_, res$pval)
    }else{
      res <- generate_odds_ratios(mr(protein.LADA, method_list = 'mr_ivw'))
      beta_ivw_ <- c(beta_ivw_, round(res$b, 4))
      OR_ivw_ <- c(OR_ivw_, round(res$or, 4))
      se_ivw_ <- c(se_ivw_, round(res$se, 4))
      CI_95_ivw_ <- c(CI_95_ivw_, paste0(round(res$or_lci95, 4), '-', round(res$or_uci95, 4)))
      P_ivw_ <- c(P_ivw_, res$pval)
    }
  }
  if (cnt %% 200 == 0){
    cat(paste0('\n----------------', cnt,'/1788-----------------\n'))
    results <- data.frame(filename=files, SNPs=SNPs, Beta=beta, OR=OR, CI_95=CI_95, P_value=P,
                          Beta_ivw_=beta_ivw_, OR_ivw_=OR_ivw_, CI_95_ivw_=CI_95_ivw_, P_ivw_=P_ivw_)
    write.csv(results, 'pQTL_LADA_5e-8_cis.csv', quote = F, row.names = F)
  }
}

results <- data.frame(filename=files, SNPs=SNPs, Beta=beta, OR=OR, SE=se, CI_95=CI_95, P_value=P,
                      Beta_ivw_=beta_ivw_, SE_ivw_=se_ivw_, OR_ivw_=OR_ivw_, CI_95_ivw_=CI_95_ivw_, P_ivw_=P_ivw_)
write.csv(results, 'pQTL_LADA_5e-8_cis.csv', quote = F, row.names = F)


#### extract IVs ####

file <- readxl::read_excel("D:/代做/药物和LADA/蛋白组_LADA结果.xlsx", sheet='deCODE_pQTL_LADA_5e-8_cis')
filenames <- file$filename[1:7]

Decode_genes <- read.table('D:/代做/中介-蛋白组-MM/data/Decode_genes_with_position.tsv',
                           sep='\t', header = T, stringsAsFactors = F)
name <- NULL
SNP <- NULL
EA <- NULL
OA <- NULL
EAF <- NULL
Pval <- NULL
Beta <- NULL
SE <- NULL
F.stat <- NULL

for (filename in filenames){
  #### cis-acting pQTL ####
  input_file <- paste0('D:/代做/中介-蛋白组-MM/data/pval_5e-08_cis-acting/', filename)
  seqid <- paste0(unlist(str_split(filename, '_'))[1], '_', unlist(str_split(filename, '_'))[2])
  gene_name <- Decode_genes$Gene[Decode_genes$SeqId==seqid]
  fp <- read.table(input_file, sep=' ', header = T, stringsAsFactors = F)
  protein <- format_data(fp, snp_col = 'rsids', effect_allele_col = 'effectAllele',
                         other_allele_col = 'otherAllele', beta_col = 'Beta',
                         se_col = 'SE', pval_col = 'Pval', eaf_col = 'effectAlleleFreq')
  LADA <- LADA.data[LADA.data$variant_id%in%protein$SNP, ]
  if (is.null(LADA)){
    cat('LADA IVs个数不够\n')
    next
  }
  if (nrow(LADA)<1){
    cat('LADA IVs个数不够\n')
    next
  }
  LADA <- format_data(LADA, type='outcome', snp_col = 'variant_id', effect_allele_col = 'reference_allele',
                      other_allele_col = 'other_allele', eaf_col = 'eaf',
                      beta_col = 'beta', se_col = 'OR_se', pval_col = 'P')
  protein.LADA <- harmonise_data(protein, LADA)
  protein.LADA$F.stat <- ((protein.LADA$beta.exposure)/(protein.LADA$se.exposure))^2
  name <- c(name, rep(gene_name, nrow(protein.LADA)))
  SNP <- c(SNP, protein.LADA$SNP)
  EA <- c(EA, protein.LADA$effect_allele.exposure)
  OA <- c(OA, protein.LADA$other_allele.exposure)
  EAF <- c(EAF, protein.LADA$eaf.exposure)
  Beta <- c(Beta, protein.LADA$beta.exposure)
  SE <- c(SE, protein.LADA$se.exposure)
  Pval <- c(Pval, protein.LADA$pval.exposure)
  F.stat <- c(F.stat, protein.LADA$F.stat)
}

results <- data.frame(Gene=name, SNP=SNP, EA=EA, OA=OA, EAF=EAF, Beta=Beta, SE=SE, Pval=Pval)
write.csv(results, '../pQTL_deCODE/deCODE_genetic_IVs.csv', quote = F, row.names = F)


#### pQTL replicted ####

#### MIA chr19:41281065..41283395
up <- 41283395 + 1000000
down <- 41281446 - 1000000
MIA <- read.table("F:/pQTL/LADA/MIA.tsv", sep='\t', header = T, stringsAsFactors = F)
MIA <- MIA[MIA$p_value<5e-8, ]
MIA <- MIA[(MIA$effect_allele_frequency>0.01)&(MIA$effect_allele_frequency<0.99),]
MIA <- MIA[MIA$chromosome==19, ]
MIA <- MIA[(MIA$base_pair_location<up)&(MIA$base_pair_location>down), ]
dat <- data.frame(rsid=MIA$variant_id, pval=MIA$p_value)
retained_SNPs <- ieugwasr::ld_clump(dat = dat, clump_kb = 1000, clump_r2 = 0.001, plink_bin = 'D:/LD_reference/plink.exe',
                                    bfile = 'D:/LD_reference/EUR')
MIA <- MIA[MIA$variant_id%in%retained_SNPs$rsid, ]
MIA <- format_data(MIA, type='exposure', snp_col = 'variant_id', pval_col = 'p_value',
                   chr_col = 'chromosome', pos_col = 'base_pair_location', effect_allele_col = 'effect_allele',
                   other_allele_col = 'other_allele', eaf_col='effect_allele_frequency', beta_col = 'beta',
                   se_col = 'standard_error')
LADA <- LADA.data[LADA.data$variant_id%in%MIA$SNP, ]
LADA <- format_data(LADA, type='outcome', snp_col = 'variant_id', effect_allele_col = 'reference_allele',
                    other_allele_col = 'other_allele', eaf_col = 'eaf',
                    beta_col = 'beta', se_col = 'OR_se', pval_col = 'P')
MIA.LADA <- harmonise_data(MIA, LADA)
generate_odds_ratios(mr(MIA.LADA))

mr_forest_plot(mr_singlesnp(MIA.LADA), exponentiate = T)

#### SAA1 chr11:18287811..18291514

down <- 18259783 - 1000000
up <- 18270215 + 1000000
SAA1 <- read.csv('../pQTL_all/emilsson_2018_for_mr.csv', sep=',', stringsAsFactors = F, header = T)
SAA1 <- SAA1[SAA1$Protein=='SAA1', ]
colnames(SAA1)[1] <- 'SNP'
SAA1 <- format_data(SAA1, snp_col = 'SNP', effect_allele_col = 'effect_allele', other_allele_col = 'other_allele',
                    eaf_col = 'eaf', beta_col = 'Beta.coeff', se_col = 'Stderr', pval_col = 'P.value',
                    samplesize_col = 'N.subjects')
LADA <- LADA.data[LADA.data$variant_id%in%SAA1$SNP, ]
LADA <- format_data(LADA, type='outcome', snp_col = 'variant_id', effect_allele_col = 'reference_allele',
                    other_allele_col = 'other_allele', eaf_col = 'eaf',
                    beta_col = 'beta', se_col = 'OR_se', pval_col = 'P', samplesize_col = 'n_samples')
SAA1.LADA <- harmonise_data(SAA1, LADA)
steiger_filtering(SAA1.LADA)
generate_odds_ratios(mr(SAA1.LADA))

#### SAA2 chr11:18259783..18270215

down <- 18259783 - 1000000
up <- 18270215 + 1000000
SAA2 <- read.csv("gwas-association-downloaded_2024-01-08-accessionId_GCST90249417.tsv", 
                 sep='\t', stringsAsFactors = F, header = T)
SAA2 <- SAA2[,-c(1:11)]
SAA2 <- SAA2[c(2,4), ]
SAA2$beta <- SAA2$OR.or.BETA
SAA2$pval <- SAA2$P.VALUE
SAA2$se <- get_se(SAA2$beta, SAA2$pval)
SAA2$effect_allele <- c('A', 'A')
SAA2$other_allele <- c('C', 'G')
SAA2$N <- 10708
SAA2 <- format_data(SAA2, snp_col = 'SNPS', effect_allele_col = 'effect_allele', other_allele_col = 'other_allele',
                    eaf_col = 'RISK.ALLELE.FREQUENCY', beta_col = 'beta', se_col = 'se', pval_col = 'pval',
                    samplesize_col = 'N')
LADA <- LADA.data[LADA.data$variant_id%in%SAA2$SNP, ]
LADA <- format_data(LADA, type='outcome', snp_col = 'variant_id', effect_allele_col = 'reference_allele',
                    other_allele_col = 'other_allele', eaf_col = 'eaf',
                    beta_col = 'beta', se_col = 'OR_se', pval_col = 'P', samplesize_col = 'n_samples')
SAA2.LADA <- harmonise_data(SAA2, LADA)
steiger_filtering(SAA2.LADA)
generate_odds_ratios(mr(SAA2.LADA))

#### CXCL10 chr4:76942271..76944650

down <- 76942271 - 1000000
up <- 76944650 + 1000000
CXCL10 <- read.csv('../pQTL_all/emilsson_2018_for_mr.csv', sep=',', stringsAsFactors = F, header = T)
CXCL10 <- CXCL10[CXCL10$Protein=='CXCL10', ]
colnames(CXCL10)[1] <- 'SNP'
CXCL10 <- format_data(CXCL10, snp_col = 'SNP', effect_allele_col = 'effect_allele', other_allele_col = 'other_allele',
                      eaf_col = 'eaf', beta_col = 'Beta.coeff', se_col = 'Stderr', pval_col = 'P.value',
                      samplesize_col = 'N.subjects')
LADA <- LADA.data[LADA.data$variant_id%in%CXCL10$SNP, ]
LADA <- format_data(LADA, type='outcome', snp_col = 'variant_id', effect_allele_col = 'reference_allele',
                    other_allele_col = 'other_allele', eaf_col = 'eaf',
                    beta_col = 'beta', se_col = 'OR_se', pval_col = 'P', samplesize_col = 'n_samples')
CXCL10.LADA <- harmonise_data(CXCL10, LADA)
steiger_filtering(CXCL10.LADA)
generate_odds_ratios(mr(CXCL10.LADA))
