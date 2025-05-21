library(susieR)
library(readr)
library(stringr)
library(tidyverse)
library(ieugwasr)
library(NMOF)

#### rs6532086 ####
### grch 37 ###
### up and down 250 kb ####
Chr <- 4
grch37_pos <- 76937278
region <- 10000

grch37_up <- grch37_pos + region
grch37_down <- grch37_pos - region

### up and down 250 kb ####
Chr <- 4
grch38_pos <- 76016125
region <- 10000

grch38_up <- grch38_pos + region
grch38_down <- grch38_pos - region

#### CXCL10 ####
#### 250 kb ####
CXCL10 <- read_tsv(gzfile('E:/pQTL/LADA/4141_79_CXCL10_IP_10.txt.gz'), col_types=cols(.default="c"))
CXCL10 <- CXCL10[CXCL10$Chrom=='chr4', ]
CXCL10$SE <- as.numeric(CXCL10$SE)
CXCL10$Beta <- as.numeric(CXCL10$Beta)
CXCL10$Pos <- as.numeric(CXCL10$Pos)
CXCL10.data <- CXCL10[(CXCL10$Pos<grch38_up)&(CXCL10$Pos>grch38_down), ]
CXCL10.data <- CXCL10.data[!is.na(CXCL10.data$rsids), ]
CXCL10.data <- CXCL10.data[!duplicated(CXCL10.data$rsids), ]
### estimate LD matrix ####
LD_matrix <- ld_matrix(CXCL10.data$rsids, bfile = 'D:/LD_reference/EUR', plink_bin = 'D:/LD_reference/plink/plink', with_alleles = T)

LD_matrix.sqr <- LD_matrix^2
logic_value <- LD_matrix.sqr[grepl('rs6532086', colnames(LD_matrix.sqr)),]>=0.1
extracted_cols <- colnames(LD_matrix.sqr)[logic_value==T]
LD_matrix <- LD_matrix[row.names(LD_matrix)%in%extracted_cols, colnames(LD_matrix)%in%extracted_cols]

rsids_in_matrix <- sapply(row.names(LD_matrix), function(x) strsplit(x,'_')[[1]][1])
CXCL10.data <- CXCL10.data[CXCL10.data$rsids%in%rsids_in_matrix, ]
CXCL10.data <- CXCL10.data[order(match(CXCL10.data$rsids, rsids_in_matrix)), ]

r <- LD_matrix[grepl('rs6532086', row.names(LD_matrix)), ];r
r <- data.frame(rsids=names(r), r=unname(r))
r$rsids <- sapply(r$rsids, function(x) strsplit(x, '_')[[1]][1])

# ### flip ###
# ld_alleles <- data.frame(rsids=rsids_in_matrix, ld_allele1=sapply(row.names(LD_matrix), function(x) strsplit(x, '_')[[1]][2]),
#                          ld_allele2=sapply(row.names(LD_matrix), function(x) strsplit(x, '_')[[1]][3]))
# merged <- merge(CXCL10.data[,c('rsids', 'effectAllele', 'otherAllele')], ld_alleles, by = 'rsids')
# merged <- merged[order(match(merged$rsids, rsids_in_matrix)),]
# merged$flip <- with(merged, (
#   (effectAllele == ld_allele2 & otherAllele == ld_allele1)
# ))
# sum(merged$flip)

estimate_s_rss(CXCL10.data$Beta/CXCL10.data$SE, R=LD_matrix, n=35559)

fitted_rss1 <- susie_rss(bhat = CXCL10.data$Beta, shat = CXCL10.data$SE, R=LD_matrix, n=35559, L=5)
sets_95_CI <- print(fitted_rss1$sets)
fitted_rss1$pip[sets_95_CI$cs$L1]
fitted_rss1$pip[grepl('rs6532086', names(fitted_rss1$pip))]

pip <- fitted_rss1$pip
pip <- data.frame(rsids=names(pip), pip=unname(pip))
pip$rsids <- sapply(pip$rsids, function(x) strsplit(x, '_')[[1]][1])

susie_res <- merge(CXCL10.data, r, by='rsids')
susie_res <- merge(susie_res, pip, by='rsids')
write.csv(susie_res, 'SuSiE_fine_mapping_results.csv', quote = F, row.names = F)