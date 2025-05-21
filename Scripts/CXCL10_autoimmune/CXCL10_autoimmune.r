##call necessary libraries
library(TwoSampleMR)
library(ieugwasr)
library("readxl") #plesae install the package

Sys.setenv(OPENGWAS_JWT="eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiJ3enFfaW5kdXN0cmlvdXNAMTYzLmNvbSIsImlhdCI6MTc0MTg3NDM0NSwiZXhwIjoxNzQzMDgzOTQ1fQ.v0g32-1WYowaakmDaf6R3OCZBA5Vjg-x_VLNcvMp9tbSI1enYifiUUL8ya8Cecf5d1MccI2UdQuA_Zqwx5sPQewzxwTz7GvsBJygaz7VcoZCihh1mmUQXOMZEQBAdaPCRgd45sPEj1pVT5bLf7rziItptopxTkNH6__u7cWBqBJ-erUqkSdiJMcIKlSyleJN4SlD1UwjKbNOrOK1Xmmsaf0zsbfy56lhVeAMw8v3u3V6mz5QYmvyymyMt4E6Cq-XHpv6EwjLz4GOk6CQ4ROqW5i5gSZUM6pAdqJY1sS5AauHX9WJyGo-8CN3xVxKUGC8a0RD3W0NHC4iALypS68t2A")
#### CXCL10 ####
CXCL10 <- read.table('D:/Others/中介-蛋白组-MM/data/pval_5e-08_cis-acting/4141_79_CXCL10_IP_10.txt', 
                     sep=' ', header=T, stringsAsFactors = F)
CXCL10 <- CXCL10[CXCL10$effectAlleleFreq>0.01,]
CXCL10 <- format_data(CXCL10, snp_col = 'rsids', effect_allele_col = 'effectAllele',
                      other_allele_col = 'otherAllele', beta_col = 'Beta',
                      se_col = 'SE', pval_col = 'Pval', eaf_col = 'effectAlleleFreq')

### primary biliary cirrhosis ####
outcome_dat<- extract_outcome_data(snps = CXCL10$SNP, outcomes = 'ebi-a-GCST90061440')
write.table(outcome_dat, 'D:/Others/药物和LADA/data/autoimmune_diseases/Primary_Biliary_Cirrhosis.txt',
            sep='\t', quote = F, row.names = F)

PBC_dat <- read.table('D:/Others/药物和LADA/data/autoimmune_diseases/Primary_Biliary_Cirrhosis.txt',
                          sep='\t', header = T, stringsAsFactors = F)
CXCL10.PBC <- harmonise_data(CXCL10, PBC_dat)
CXCL10.PBC$mr_keep <- T
generate_odds_ratios(mr(CXCL10.PBC))

### multiple sclerosis ####
outcome_dat<- extract_outcome_data(snps = CXCL10$SNP, outcomes = 'ieu-b-18')
write.table(outcome_dat, 'D:/Others/药物和LADA/data/autoimmune_diseases/Multiple_Sclerosis.txt',
            sep='\t', quote = F, row.names = F)

MS_dat <- read.table('D:/Others/药物和LADA/data/autoimmune_diseases/Multiple_Sclerosis.txt',
                     sep='\t', header = T, stringsAsFactors = F)
CXCL10.MS <- harmonise_data(CXCL10, MS_dat)
generate_odds_ratios(mr(CXCL10.MS))

### Crohn's disease ####
outcome_dat<- extract_outcome_data(snps = CXCL10$SNP, outcomes = 'ebi-a-GCST004132')
write.csv(outcome_dat, 'D:/Others/药物和LADA/data/autoimmune_diseases/Crohn_disease.csv',
          quote = F, row.names = F)

CD_dat <- read.csv('D:/Others/药物和LADA/data/autoimmune_diseases/Crohn_disease.csv',
                     sep=',', header = T, stringsAsFactors = F)
CXCL10.CD <- harmonise_data(CXCL10, CD_dat)
generate_odds_ratios(mr(CXCL10.CD))

### Ulcerative Colitis ####
outcome_dat<- extract_outcome_data(snps = CXCL10$SNP, outcomes = 'ebi-a-GCST90038684')
write.csv(outcome_dat, 'D:/Others/药物和LADA/data/autoimmune_diseases/Ulcerative_colites.csv',
          quote = F, row.names = F)

UC_dat <- read.csv('D:/Others/药物和LADA/data/autoimmune_diseases/Ulcerative_colites.csv',
                   sep=',', header = T, stringsAsFactors = F)
CXCL10.UC <- harmonise_data(CXCL10, UC_dat)
generate_odds_ratios(mr(CXCL10.UC))


### SLE ####
SLE <- read.table('D:/Others/药物和LADA/data/autoimmune_diseases/SLE/Spain_Results/Spain_Results.txt',
                  header = T, stringsAsFactors = F, sep='\t')
SLE <- SLE[SLE$SNP%in%CXCL10$SNP, ]
SLE$SE <- (log(SLE$OR_upper) - log(SLE$OR))/qnorm(0.975)
SLE$BETA <- log(SLE$OR)
SLE <- format_data(SLE, type='outcome', snp_col = 'SNP', effect_allele_col = 'A1leleA',
                   other_allele_col = 'AlleleB', beta_col = 'BETA', se_col = 'SE',
                   pval_col = 'P')
CXCL10.SLE <- harmonise_data(CXCL10, SLE)
generate_odds_ratios(mr(CXCL10.SLE))

### GRAVE's disease ####
finngen_cols <- c('#chrom', 'pos', 'ref', 'alt', 'rsids', 'nearest_genes',
                  'pval', 'mlogp', 'beta', 'sebeta', 'af_alt', 'af_alt_cases',
                  'af_alt_controls')
Grave <- read.table(gzfile('D:/Others/药物和LADA/data/autoimmune_diseases/summary_stats_release_finngen_R12_E4_GRAVES_STRICT.gz'),
                    header=F, stringsAsFactors = F, sep='\t')
colnames(Grave) <- finngen_cols
Grave <- Grave[Grave$rsids%in%CXCL10$SNP, ]
Grave <- format_data(Grave, type='outcome', snp_col = 'rsids', effect_allele_col = 'alt',
                     other_allele_col = 'ref', beta_col = 'beta', se_col = 'sebeta',
                     pval_col = 'pval', eaf_col = 'af_alt')
CXCL10.Grave <- harmonise_data(CXCL10, Grave)
generate_odds_ratios(mr(CXCL10.Grave))

### celiac's disease ####
celiac <- read.table(gzfile('D:/Others/药物和LADA/data/autoimmune_diseases/summary_stats_release_finngen_R12_K11_COELIAC.gz'),
                     header=F, stringsAsFactors = F, sep='\t')
colnames(celiac) <- finngen_cols
celiac <- celiac[celiac$rsids%in%CXCL10$SNP, ]
celiac <- format_data(celiac, type='outcome', snp_col = 'rsids', effect_allele_col = 'alt',
                      other_allele_col = 'ref', beta_col = 'beta', se_col = 'sebeta',
                      pval_col = 'pval', eaf_col = 'af_alt')
CXCL10.celiac <- harmonise_data(CXCL10, celiac)
generate_odds_ratios(mr(CXCL10.celiac))

### psoriasis disease ####
psoriasis <- read.table(gzfile('D:/Others/药物和LADA/data/autoimmune_diseases/summary_stats_release_finngen_R12_L12_PSORIASIS.gz'),
                        header=F, stringsAsFactors = F, sep='\t', col.names = finngen_cols)
psoriasis <- psoriasis[psoriasis$rsids%in%CXCL10$SNP, ]
psoriasis <- format_data(psoriasis, type='outcome', snp_col = 'rsids', effect_allele_col = 'alt',
                         other_allele_col = 'ref', beta_col = 'beta', se_col = 'sebeta',
                         pval_col = 'pval', eaf_col = 'af_alt')
CXCL10.psoriasis <- harmonise_data(CXCL10, psoriasis)
generate_odds_ratios(mr(CXCL10.psoriasis))

### rheuma_arth disease ####
rheuma_arth <- read.table(gzfile('D:/Others/药物和LADA/data/autoimmune_diseases/summary_stats_release_finngen_R12_M13_RHEUMA.gz'),
                          header=F, stringsAsFactors = F, sep='\t', col.names = finngen_cols)
rheuma_arth <- rheuma_arth[rheuma_arth$rsids%in%CXCL10$SNP, ]
rheuma_arth <- format_data(rheuma_arth, type='outcome', snp_col = 'rsids', effect_allele_col = 'alt',
                           other_allele_col = 'ref', beta_col = 'beta', se_col = 'sebeta',
                           pval_col = 'pval', eaf_col = 'af_alt')
CXCL10.rheuma_arth <- harmonise_data(CXCL10, rheuma_arth)
generate_odds_ratios(mr(CXCL10.rheuma_arth))

### sjogren disease ####
sjogren <- read.table(gzfile('D:/Others/药物和LADA/data/autoimmune_diseases/summary_stats_release_finngen_R12_M13_SJOGREN.gz'),
                      header=F, stringsAsFactors = F, sep='\t', col.names = finngen_cols)
sjogren <- sjogren[sjogren$rsids%in%CXCL10$SNP, ]
sjogren <- format_data(sjogren, type='outcome', snp_col = 'rsids', effect_allele_col = 'alt',
                       other_allele_col = 'ref', beta_col = 'beta', se_col = 'sebeta',
                       pval_col = 'pval', eaf_col = 'af_alt')
CXCL10.sjogren <- harmonise_data(CXCL10, sjogren)
generate_odds_ratios(mr(CXCL10.sjogren))

p.adjust(c(0.005, 0.098, 0.006, 0.707, 0.568, 0.827, 0.726, 0.077, 0.507, 0.024), method = 'fdr')

### IBD disease ####
IBD <- read.table(gzfile('D:/Others/药物和LADA/data/autoimmune_diseases/summary_stats_release_finngen_R12_K11_IBD_STRICT.gz'),
                  header=F, stringsAsFactors = F, sep='\t', col.names = finngen_cols)
IBD <- IBD[IBD$rsids%in%CXCL10$SNP, ]
IBD <- format_data(IBD, type='outcome', snp_col = 'rsids', effect_allele_col = 'alt',
                   other_allele_col = 'ref', beta_col = 'beta', se_col = 'sebeta',
                   pval_col = 'pval', eaf_col = 'af_alt')
CXCL10.IBD <- harmonise_data(CXCL10, IBD)
generate_odds_ratios(mr(CXCL10.IBD))
