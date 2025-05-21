##call necessary libraries
library(TwoSampleMR)
library(ieugwasr)
library("readxl") #plesae install the package

Sys.setenv(OPENGWAS_JWT="eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiJ3enFfaW5kdXN0cmlvdXNAMTYzLmNvbSIsImlhdCI6MTc0MTY1OTUwMCwiZXhwIjoxNzQyODY5MTAwfQ.m4XjtoJQm6QaaeU24yR9BVB5OPT2EgsYMQd_kdF20p5Ay0VsZaoySYm5mNLSFGOBrafN6nYzfbuhKRlveNVU5wBgykijCJUU67FFPOV3flYWdaL37vnsPzz-nFfElEhxAURpB9p_-g-dvJq9K6yJm4u2B0RZOQy0t5yIK4Zui4gtKYKNUK-uSCny8osuCVTIE_TZ_vR22YhUBXyFNaka4piNz4iJQ4Ei757fIs4RKFWrX505DS8SosiX9WIx-PphK6eqx2K9DC6FX3XE-o-e9TWuyjLdOkg4GXG1EQgWYBl0GHy78vwJ6h_HVeOhWv_09asQXgUjKAXOzeIvgOmAaQ")
outcome_id <- read_excel('OpenGWAS_ids.xlsx')
ids<-as.character(outcome_id$GWAS_id)

#### CXCL10 ####
CXCL10 <- read.table('D:/Others/中介-蛋白组-MM/data/pval_5e-08_cis-acting/4141_79_CXCL10_IP_10.txt', 
                     sep=' ', header=T, stringsAsFactors = F)
CXCL10 <- format_data(CXCL10, snp_col = 'rsids', effect_allele_col = 'effectAllele',
                      other_allele_col = 'otherAllele', beta_col = 'Beta',
                      se_col = 'SE', pval_col = 'Pval', eaf_col = 'effectAlleleFreq')
CXCL10 <- CXCL10[CXCL10$eaf.exposure>0.01, ]
outcome_dat<- extract_outcome_data(snps = CXCL10$SNP, outcomes = ids)
dat <- NULL
dat <- harmonise_data(
  exposure_dat = CXCL10, 
  outcome_dat = outcome_dat
)
mr_results <- NULL
try(mr_results <- mr(dat, method_list=c("mr_wald_ratio", "mr_ivw")))  # main MR analysis
result_file <- "Phewas.CXCL10.mr_1024.txt"
if (exists("mr_results")==TRUE){ 
  write.table(mr_results, file=result_file, sep="\t", col.names=T, row.names=F, quote=F)
}
