library(dplyr)
library(readr)
library(TwoSampleMR)
library(tidyverse)
library(ieugwasr)
library(MRPRESSO)

#### analysis with type 1 diabetes and type 2 diabetes ####

#### ebi-a-GCST90014023 T1D ####

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

file <- readxl::read_excel("D:/Others/药物和LADA/蛋白组_LADA结果_20240206.xlsx", sheet='comparison_LADA_other_diabetes')
filenames <- file$filename

Decode_genes <- read.table('D:/Others/中介-蛋白组-MM/data/Decode_genes_with_position.tsv',
                           sep='\t', header = T, stringsAsFactors = F)
input_file <- paste0('D:/Others/中介-蛋白组-MM/data/pval_5e-08_cis-acting/', filenames[1])
fp <- read.table(input_file, sep=' ', header = T, stringsAsFactors = F)

for (filename in filenames){
  #### cis-acting pQTL ####
  input_file <- paste0('D:/Others/中介-蛋白组-MM/data/pval_5e-08_cis-acting/', filename)
  seqid <- paste0(unlist(str_split(filename, '_'))[1], '_', unlist(str_split(filename, '_'))[2])
  gene_name <- Decode_genes$Gene[Decode_genes$SeqId==seqid]
  fp_add <- read.table(input_file, sep=' ', header = T, stringsAsFactors = F)
  fp <- rbind(fp, fp_add)
}

fp <- fp[!duplicated(fp$rsids), ]

T1D <- extract_outcome_data(fp$rsids, 'ebi-a-GCST90014023')
write.csv(T1D, 'T1D_summary_statistics_GCST90014023.csv', quote = F, row.names = F)

T1D.data <- read.csv('T1D_summary_statistics_GCST90014023.csv', sep = ',', stringsAsFactors = F, header = T)

for (filename in filenames){
  #### cis-acting pQTL ####
  input_file <- paste0('D:/Others/中介-蛋白组-MM/data/pval_5e-08_cis-acting/', filename)
  seqid <- paste0(unlist(str_split(filename, '_'))[1], '_', unlist(str_split(filename, '_'))[2])
  gene_name <- Decode_genes$Gene[Decode_genes$SeqId==seqid]
  fp <- read.table(input_file, sep=' ', header = T, stringsAsFactors = F)
  fp <- fp[fp$effectAlleleFreq>0.01, ]
  protein <- format_data(fp, snp_col = 'rsids', effect_allele_col = 'effectAllele',
                         other_allele_col = 'otherAllele', beta_col = 'Beta',
                         se_col = 'SE', pval_col = 'Pval', eaf_col = 'effectAlleleFreq')
  T1D <- T1D.data[T1D.data$SNP%in%protein$SNP, ]
  if (is.null(T1D)){
    cat('T1D IVs个数不够\n')
    next
  }
  if (nrow(T1D)<1){
    cat('T1D IVs个数不够\n')
    next
  }
  protein.T1D <- harmonise_data(protein, T1D)
  if (nrow(protein.T1D[protein.T1D$mr_keep==T, ])==0){
    cat('harmonize IVs个数不够\n')
  }else if (nrow(protein.T1D[protein.T1D$mr_keep==T, ])==1){
    res <- generate_odds_ratios(mr(protein.T1D, method_list = 'mr_wald_ratio'))
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
    res <- generate_odds_ratios(mr(protein.T1D, method_list = 'mr_ivw'))
    files <- c(files, filename)
    SNPs <- c(SNPs, res$nsnp)
    beta <- c(beta, round(res$b, 4))
    se <- c(se, round(res$se, 4))
    OR <- c(OR, round(res$or, 4))
    CI_95 <- c(CI_95, paste0(round(res$or_lci95, 4), '-', round(res$or_uci95, 4)))
    P <- c(P, res$pval)
    Q <- mr_heterogeneity(protein.T1D, method_list = 'mr_ivw')
    Q <- Q$Q_pval
    if (Q>0.05){
      res <- generate_odds_ratios(mr(protein.T1D, method_list = 'mr_ivw_fe'))
      beta_ivw_ <- c(beta_ivw_, round(res$b, 4))
      se_ivw_ <- c(se_ivw_, round(res$se, 4))
      OR_ivw_ <- c(OR_ivw_, round(res$or, 4))
      CI_95_ivw_ <- c(CI_95_ivw_, paste0(round(res$or_lci95, 4), '-', round(res$or_uci95, 4)))
      P_ivw_ <- c(P_ivw_, res$pval)
    }else{
      res <- generate_odds_ratios(mr(protein.T1D, method_list = 'mr_ivw'))
      beta_ivw_ <- c(beta_ivw_, round(res$b, 4))
      OR_ivw_ <- c(OR_ivw_, round(res$or, 4))
      se_ivw_ <- c(se_ivw_, round(res$se, 4))
      CI_95_ivw_ <- c(CI_95_ivw_, paste0(round(res$or_lci95, 4), '-', round(res$or_uci95, 4)))
      P_ivw_ <- c(P_ivw_, res$pval)
    }
  }
}

results <- data.frame(filename=files, SNPs=SNPs, Beta=beta, OR=OR, SE=se, CI_95=CI_95, P_value=P,
                      Beta_ivw_=beta_ivw_, SE_ivw_=se_ivw_, OR_ivw_=OR_ivw_, CI_95_ivw_=CI_95_ivw_, P_ivw_=P_ivw_)
write.csv(results, 'pQTL_T1D_5e-8_cis_GCST90014023.csv', quote = F, row.names = F)

#### 2018 T2D ####

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

T2D.data <- read.table('Mahajan.NatGenet2018b.T2D.European.txt', sep='\t', header = T, stringsAsFactors = F)

file <- readxl::read_excel("D:/Others/药物和LADA/蛋白组_LADA结果_20240206.xlsx", sheet='comparison_LADA_other_diabetes')
filenames <- file$filename

Decode_genes <- read.table('D:/Others/中介-蛋白组-MM/data/Decode_genes_with_position.tsv',
                           sep='\t', header = T, stringsAsFactors = F)
input_file <- paste0('D:/Others/中介-蛋白组-MM/data/pval_5e-08_cis-acting/', filenames[1])
fp <- read.table(input_file, sep=' ', header = T, stringsAsFactors = F)


for (filename in filenames[2:100]){
  #### cis-acting pQTL ####
  input_file <- paste0('D:/Others/中介-蛋白组-MM/data/pval_5e-08_cis-acting/', filename)
  seqid <- paste0(unlist(str_split(filename, '_'))[1], '_', unlist(str_split(filename, '_'))[2])
  gene_name <- Decode_genes$Gene[Decode_genes$SeqId==seqid]
  fp_add <- read.table(input_file, sep=' ', header = T, stringsAsFactors = F)
  fp <- rbind(fp, fp_add)
}

fp <- fp[!duplicated(fp$rsids), ]

write.csv(fp, 'LADA-associated_proteins_IVs.csv', row.names = F, quote = F)

fp <- read.csv('LADA-associated_proteins_IVs.csv', header = T, stringsAsFactors = F)
fp <- fp[fp$effectAlleleFreq>0.01, ]
T2D.data$variant <- paste0(T2D.data$Chr, ':', T2D.data$Pos)
fp$variant <- paste0(fp$Chr,':',fp$Pos)
T2D.data <- T2D.data[T2D.data$variant%in%fp$variant, ]
rsids <- data.frame(variant=fp$variant, rsids=fp$rsids)
T2D.data <- merge(T2D.data, rsids, by='variant', all.x=T)

for (filename in filenames){
  #### cis-acting pQTL ####
  input_file <- paste0('D:/Others/中介-蛋白组-MM/data/pval_5e-08_cis-acting/', filename)
  seqid <- paste0(unlist(str_split(filename, '_'))[1], '_', unlist(str_split(filename, '_'))[2])
  gene_name <- Decode_genes$Gene[Decode_genes$SeqId==seqid]
  fp <- read.table(input_file, sep=' ', header = T, stringsAsFactors = F)
  protein <- format_data(fp, snp_col = 'rsids', effect_allele_col = 'effectAllele',
                         other_allele_col = 'otherAllele', beta_col = 'Beta',
                         se_col = 'SE', pval_col = 'Pval', eaf_col = 'effectAlleleFreq')
  T2D <- T2D.data[T2D.data$rsids%in%protein$SNP, ]
  if (is.null(T2D)){
    cat('T2D IVs个数不够\n')
    next
  }
  if (nrow(T2D)<1){
    cat('T2D IVs个数不够\n')
    next
  }
  T2D <- format_data(T2D, type='outcome', snp_col='rsids', effect_allele_col='EA',
                     other_allele_col='NEA', eaf_col='EAF',
                     se_col='SE', beta_col = 'Beta', pval_col = 'Pvalue')
  protein.T2D <- harmonise_data(protein, T2D)
  if (nrow(protein.T2D[protein.T2D$mr_keep==T, ])==0){
    cat('harmonize IVs个数不够\n')
  }else if (nrow(protein.T2D[protein.T2D$mr_keep==T, ])==1){
    res <- generate_odds_ratios(mr(protein.T2D, method_list = 'mr_wald_ratio'))
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
    res <- generate_odds_ratios(mr(protein.T2D, method_list = 'mr_ivw'))
    files <- c(files, filename)
    SNPs <- c(SNPs, res$nsnp)
    beta <- c(beta, round(res$b, 4))
    se <- c(se, round(res$se, 4))
    OR <- c(OR, round(res$or, 4))
    CI_95 <- c(CI_95, paste0(round(res$or_lci95, 4), '-', round(res$or_uci95, 4)))
    P <- c(P, res$pval)
    Q <- mr_heterogeneity(protein.T2D, method_list = 'mr_ivw')
    Q <- Q$Q_pval
    if (Q>0.05){
      res <- generate_odds_ratios(mr(protein.T2D, method_list = 'mr_ivw_fe'))
      beta_ivw_ <- c(beta_ivw_, round(res$b, 4))
      se_ivw_ <- c(se_ivw_, round(res$se, 4))
      OR_ivw_ <- c(OR_ivw_, round(res$or, 4))
      CI_95_ivw_ <- c(CI_95_ivw_, paste0(round(res$or_lci95, 4), '-', round(res$or_uci95, 4)))
      P_ivw_ <- c(P_ivw_, res$pval)
    }else{
      res <- generate_odds_ratios(mr(protein.T2D, method_list = 'mr_ivw'))
      beta_ivw_ <- c(beta_ivw_, round(res$b, 4))
      OR_ivw_ <- c(OR_ivw_, round(res$or, 4))
      se_ivw_ <- c(se_ivw_, round(res$se, 4))
      CI_95_ivw_ <- c(CI_95_ivw_, paste0(round(res$or_lci95, 4), '-', round(res$or_uci95, 4)))
      P_ivw_ <- c(P_ivw_, res$pval)
    }
  }
}

results <- data.frame(filename=files, SNPs=SNPs, Beta=beta, OR=OR, SE=se, CI_95=CI_95, P_value=P,
                      Beta_ivw_=beta_ivw_, SE_ivw_=se_ivw_, OR_ivw_=OR_ivw_, CI_95_ivw_=CI_95_ivw_, P_ivw_=P_ivw_)
write.csv(results, 'pQTL_T2D_5e-8_cis_2018.csv', quote = F, row.names = F)