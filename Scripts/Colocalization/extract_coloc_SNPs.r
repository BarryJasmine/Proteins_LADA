library(readr)

#### MIA ####
#### hg38: chr19:40775160..40777490
MIA <- read_tsv(gzfile('E:/pQTL/LADA/2687_2_MIA_MIA.txt.gz'), col_types=cols(.default="c"))
#### top cis-pqtl 500 kb ####
top_pqtl <- 'rs2279011'
CHR <- MIA$Chrom[(MIA$rsids==top_pqtl)&(!is.na(MIA$rsids))]
MIA_coloc <- MIA[MIA$Chrom==CHR, ]
MIA_coloc$Pos <- as.numeric(MIA_coloc$Pos)

POS <- MIA_coloc$Pos[(MIA_coloc$rsids==top_pqtl)&(!is.na(MIA_coloc$rsids))]
region <- 500000
MIA_coloc <- MIA_coloc[(MIA_coloc$Pos<POS+region)&(MIA_coloc$Pos>POS-region),]
write.table(MIA_coloc, 'colocalization_top_pqtl_500kb/MIA_coloc_SNPs.txt', sep='\t', quote = F, row.names = F)

#### CXCL10 ####
#### hg38: chr4:76021118..76023497

#### top cis-pqtl 500 kb ####

top_pqtl <- 'rs6532086'
CHR <- CXCL10$Chrom[(CXCL10$rsids==top_pqtl)&(!is.na(CXCL10$rsids))]
CXCL10_coloc <- CXCL10[CXCL10$Chrom==CHR, ]
CXCL10_coloc$Pos <- as.numeric(CXCL10_coloc$Pos)

POS <- CXCL10_coloc$Pos[(CXCL10_coloc$rsids==top_pqtl)&(!is.na(CXCL10_coloc$rsids))]
region <- 500000
CXCL10_coloc <- CXCL10_coloc[(CXCL10_coloc$Pos<POS+region)&(CXCL10_coloc$Pos>POS-region),]
write.table(CXCL10_coloc, 'colocalization_top_pqtl_500kb/CXCL10_coloc_SNPs.txt', sep='\t', quote = F, row.names = F)

#### SAA1 ####
#### hg38: chr11:18266264..18269967
#### top cis-pqtl 500 kb ####
top_pqtl <- 'rs11024589'
CHR <- SAA1$Chrom[(SAA1$rsids==top_pqtl)&(!is.na(SAA1$rsids))]
SAA1_coloc <- SAA1[SAA1$Chrom==CHR, ]
SAA1_coloc$Pos <- as.numeric(SAA1_coloc$Pos)

POS <- SAA1_coloc$Pos[(SAA1_coloc$rsids==top_pqtl)&(!is.na(SAA1_coloc$rsids))]
region <- 500000
SAA1_coloc <- SAA1_coloc[(SAA1_coloc$Pos<POS+region)&(SAA1_coloc$Pos>POS-region),]
write.table(SAA1_coloc, 'colocalization_top_pqtl_500kb/SAA1_coloc_SNPs.txt', sep='\t', quote = F, row.names = F)
#### top cis-pqtl 250 kb ####
top_pqtl <- 'rs11024589'
CHR <- SAA1$Chrom[(SAA1$rsids==top_pqtl)&(!is.na(SAA1$rsids))]
SAA1_coloc <- SAA1[SAA1$Chrom==CHR, ]
SAA1_coloc$Pos <- as.numeric(SAA1_coloc$Pos)

POS <- SAA1_coloc$Pos[(SAA1_coloc$rsids==top_pqtl)&(!is.na(SAA1_coloc$rsids))]
region <- 250000
SAA1_coloc <- SAA1_coloc[(SAA1_coloc$Pos<POS+region)&(SAA1_coloc$Pos>POS-region),]
write.table(SAA1_coloc, 'colocalization_top_pqtl_250kb/SAA1_coloc_SNPs.txt', sep='\t', quote = F, row.names = F)

#### SAA2 ####
#### hg38: chr11:18238236..18248668
SAA2 <- read_tsv(gzfile('E:/pQTL/LADA/18832_65_SAA2_SAA2.txt.gz'), col_types=cols(.default="c"))

#### top cis-pqtl 500 kb ####
top_pqtl <- 'rs11024589'
CHR <- SAA2$Chrom[(SAA2$rsids==top_pqtl)&(!is.na(SAA2$rsids))]
SAA2_coloc <- SAA2[SAA2$Chrom==CHR, ]
SAA2_coloc$Pos <- as.numeric(SAA2_coloc$Pos)

POS <- SAA2_coloc$Pos[(SAA2_coloc$rsids==top_pqtl)&(!is.na(SAA2_coloc$rsids))]
region <- 500000
SAA2_coloc <- SAA2_coloc[(SAA2_coloc$Pos<POS+region)&(SAA2_coloc$Pos>POS-region),]
write.table(SAA2_coloc, 'colocalization_top_pqtl_500kb/SAA2_coloc_SNPs.txt', sep='\t', quote = F, row.names = F)