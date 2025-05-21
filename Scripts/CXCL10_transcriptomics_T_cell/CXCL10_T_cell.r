##call necessary libraries
library(TwoSampleMR)
library(ieugwasr)
library(readxl)

#### outcomes ####
outcome_dat <- read.table('Transcriptomics_T_cell_Genes.txt', sep='\t', header = T, stringsAsFactors = F)
outcome_dat$pval <- 10^(-1*outcome_dat$LOG10P)
proteins <- outcome_dat$Protein_name
protein_names <- NULL
Beta <- NULL
SE <- NULL
CI_95 <- NULL
Pval <- NULL

#### CXCL10 ####
CXCL10 <- read.table('D:/Others/中介-蛋白组-MM/data/pval_5e-08_cis-acting/4141_79_CXCL10_IP_10.txt', 
                     sep=' ', header=T, stringsAsFactors = F)
CXCL10 <- format_data(CXCL10, snp_col = 'rsids', effect_allele_col = 'effectAllele',
                      other_allele_col = 'otherAllele', beta_col = 'Beta',
                      se_col = 'SE', pval_col = 'Pval', eaf_col = 'effectAlleleFreq')
CXCL10 <- CXCL10[CXCL10$eaf.exposure>0.01, ]

for(protein_name in proteins){
  protein <- outcome_dat[outcome_dat$Protein_name==protein_name, ]
  protein <- format_data(protein, type='outcome', snp_col='SNP', effect_allele_col = 'ALLELE1',
                         other_allele_col = 'ALLELE0', eaf_col = 'A1FREQ', beta_col = 'BETA',
                         se_col = 'SE', pval_col = 'pval')
  CXCL10.protein <- harmonise_data(CXCL10, protein)
  res <- generate_odds_ratios(mr(CXCL10.protein))
  protein_names <- c(protein_names, protein_name)
  Beta <- c(Beta, round(res$b, 4))
  SE <- c(SE, round(res$se, 4))
  CI_95 <- c(CI_95, paste0(round(res$lo_ci, 4), ' - ', round(res$up_ci, 4)))
  Pval <- c(Pval, res$pval)
}
result <- data.frame(protein_names=protein_names, Beta=Beta, CI_95=CI_95, Pval=Pval, SE=SE)
write.table(result, 'CXCL10_Transcriptomics_T_cells_one_IV.txt', quote=F, row.names=F, sep='\t')

results <- read.csv('CXCL10_Transcriptomics_T_cells_one_IV.txt', sep='\t',
                    header = T, stringsAsFactors = F)
p.adjust(results$Pval, method = 'fdr')