library(coloc)
library(snpStats)
library(locuscomparer)
library(readr)
library(TwoSampleMR)
library(stringr)
library(locusplotr)
library(tidyverse)

####### read summary genetic data of LADA ##### 
LADA.data <- read_tsv("D:/Others/药物和LADA/data/LADACTRL_MA_filtered_2018_DLC.txt", col_types=cols(.default="c"))
LADA.data[,5:10] <- sapply(LADA.data[,5:10], as.numeric)
LADA.data$beta <- log(LADA.data$OR)

LTL <- read_tsv(gzfile("D:/Others/药物和LADA/data/UKB_telomere_gwas_summarystats.tsv.gz"), col_types=cols(.default="c"))
LTL <- LTL[,c('variant_id', 'chromosome', 'base_pair_location')]
LTL$rs_number <- paste0(LTL$chromosome, ':', LTL$base_pair_location)
LTL <- LTL[,c('rs_number', 'variant_id', 'chromosome', 'base_pair_location')]

LADA.data <- merge(LADA.data, LTL, by='rs_number', all.x = T)


filenames.list <- list.files('colocalization_data_100kb/')

LADA_cc_ratio <- 2634/8581

Decode_genes <- read.table('D:/Others/中介-蛋白组-MM/data/Decode_genes_with_position.tsv',
                           sep='\t', header = T, stringsAsFactors = F)

### top pqtl 500kb ####

filenames.list <- list.files('colocalization_top_pqtl_500kb/')

LADA_cc_ratio <- 2634/8581

Decode_genes <- read.table('D:/Others/中介-蛋白组-MM/data/Decode_genes_with_position.tsv',
                           sep='\t', header = T, stringsAsFactors = F)

# set priors
p1_comb <- c(1e-4, 1e-5, 1e-6)
p2_comb <- c(1e-4, 1e-5, 1e-6)
p12_comb <- c(5e-5, 1e-5, 5e-6, 1e-6)

coloc.res <- data.frame(filename=0, Protein=0, p1=0, p2=0, p12=0, snps=0, H0=0, H1=0, H2=0, H3=0, H4=0)
for(p1 in p1_comb){
  for(p2 in p2_comb){
    for(p12 in p12_comb){
      for (filename in filenames.list){
        input_file <- paste0('colocalization_top_pqtl_500kb/', filename)
        protein <- unlist(str_split(filename, '_'))[1]
        fp <- read.table(input_file, sep='\t', header = T, stringsAsFactors = F)
        protein.data <- format_data(fp, snp_col = 'rsids', effect_allele_col = 'effectAllele',
                                    other_allele_col = 'otherAllele', beta_col = 'Beta',
                                    se_col = 'SE', pval_col = 'Pval',
                                    samplesize_col = 'N')
        LADA <- LADA.data[LADA.data$variant_id%in%protein.data$SNP, ]
        LADA <- format_data(LADA, type='outcome', pval_col = 'P', beta_col = 'beta',
                            snp_col ='variant_id', effect_allele_col = 'reference_allele', other_allele_col = 'other_allele',
                            se_col = 'OR_se', eaf_col = 'eaf', samplesize_col = 'n_samples')
        protein.LADA <- harmonise_data(protein.data, LADA)
        protein.LADA <- protein.LADA[protein.LADA$mr_keep==T, ]
        protein.LADA$var.exposure <- protein.LADA$se.exposure^2
        protein.LADA$var.outcome <- protein.LADA$se.outcome^2
        result <- coloc.abf(dataset1 = list( beta = protein.LADA$beta.exposure, # beta
                                             varbeta = protein.LADA$var.exposure,  # standard error squared
                                             N = protein.LADA$samplesize.exposure,  # sample size # might work without it
                                             type = "quant",
                                             sdY = 1,
                                             #MAF = protein.LADA$eaf.exposure,
                                             snp = protein.LADA$SNP),
                            dataset2 = list( beta = protein.LADA$beta.outcome,
                                             varbeta = protein.LADA$var.outcome,
                                             type = "cc", # case control
                                             N = protein.LADA$samplesize.outcome,
                                             s = LADA_cc_ratio,
                                             MAF = protein.LADA$eaf.outcome,
                                             snp = protein.LADA$SNP), p1 = p1, p2 = p2, p12 = p12)
        summary_coloc <- as.numeric(result$summary)
        sub_res <- c(filename, protein, p1, p2, p12, summary_coloc[1], paste0(round(summary_coloc[2:6]*100, 4), '%'))
        coloc.res <- rbind(coloc.res, sub_res)
      }
    }
  }
}

write.csv(coloc.res, 'pQTL_LADA_coloc_top_pqtl_500kb.csv', quote = F, row.names = F)

#### Locus plot ####

#### CXCL10 ####
filename <- 'CXCL10_coloc_SNPs.txt'
input_file <- paste0('colocalization_top_pqtl_500kb/', filename)
fp <- read.table(input_file, sep='\t', header = T, stringsAsFactors = F)
fp$chr <- 4
protein.data <- format_data(fp, snp_col = 'rsids', effect_allele_col = 'effectAllele',
                            other_allele_col = 'otherAllele', beta_col = 'Beta',
                            se_col = 'SE', pval_col = 'Pval',
                            samplesize_col = 'N')
LADA <- LADA.data[LADA.data$variant_id%in%protein.data$SNP, ]
LADA <- format_data(LADA, type='outcome', pval_col = 'P', beta_col = 'beta',
                    snp_col ='variant_id', effect_allele_col = 'reference_allele', other_allele_col = 'other_allele',
                    se_col = 'OR_se', eaf_col = 'eaf', samplesize_col = 'n_samples', chr_col = 'chromosome',
                    pos_col = 'base_pair_location')
protein.LADA <- harmonise_data(protein.data, LADA)
protein.LADA <- protein.LADA[protein.LADA$mr_keep==T, ]
protein.LADA$chr.outcome <- as.numeric(protein.LADA$chr.outcome)
protein.LADA$pos.outcome <- as.numeric(protein.LADA$pos.outcome)
protein.LADA$pval.exposure <- as.numeric(protein.LADA$pval.exposure)
protein <- data.frame(chromosome=protein.LADA$chr.outcome, position=protein.LADA$pos.outcome, 
                    rsid=protein.LADA$SNP, effect_allele=protein.LADA$effect_allele.exposure, 
                    other_allele=protein.LADA$other_allele.exposure,
                    effect=protein.LADA$beta.exposure, std_err=protein.LADA$se.exposure, 
                    pval=protein.LADA$pval.exposure)
p1 <- gg_locusplot(df=protein, lead_snp = 'rs6532086', chrom=chromosome, pos=position,
                   ref=effect_allele, alt=other_allele, p_value = pval, genome_build = 'GRCh37', population='EUR',
                   plot_genes = F, plot_pvalue_threshold = 1);p1
p1 <- p1 + 
  labs(x = "Chromosome 4 (Mb)", fill=parse(text = "LD  (r^2)")) +
  theme(panel.grid=element_blank(), legend.background = element_rect(fill = 'white', colour = 'black'),
        legend.position = c(0.14, 0.99));p1
Ipaper::write_fig(p1,
                  file='regional_plot/cis-pQTL for CXCL10.jpg',
                  width=10,
                  height=4,
                  devices='jpg',
                  res=300,
                  show=T)
### LADA.data
LADA.data.CXCL10 <- data.frame(chromosome=protein.LADA$chr.outcome, position=protein.LADA$pos.outcome, 
                              rsid=protein.LADA$SNP, effect_allele=protein.LADA$effect_allele.outcome, 
                              other_allele=protein.LADA$other_allele.outcome,
                              effect=protein.LADA$beta.outcome, std_err=protein.LADA$se.outcome, 
                              pval=protein.LADA$pval.outcome)
source('D:/Others/ACLY与cancer/共定位/gg_locusplot_self.r')
p2 <- gg_locusplot_self(df=LADA.data.CXCL10, lead_snp = 'rs6532086', chrom=chromosome, pos=position,
                        ref=effect_allele, alt=other_allele, p_value = pval, genome_build = 'GRCh37', population='EUR',
                        plot_genes = T, plot_pvalue_threshold = 1);p2
Ipaper::write_fig(p2,
                  file='regional_plot/CXCL10_LADA.jpg',
                  width=10,
                  height=4,
                  devices='jpg',
                  res=300,
                  show=T)

#### SAA1 ####
filename <- 'SAA1_coloc_SNPs.txt'
input_file <- paste0('colocalization_top_pqtl_500kb/', filename)
fp <- read.table(input_file, sep='\t', header = T, stringsAsFactors = F)
protein.data <- format_data(fp, snp_col = 'rsids', effect_allele_col = 'effectAllele',
                            other_allele_col = 'otherAllele', beta_col = 'Beta',
                            se_col = 'SE', pval_col = 'Pval',
                            samplesize_col = 'N')
LADA <- LADA.data[LADA.data$variant_id%in%protein.data$SNP, ]
LADA <- format_data(LADA, type='outcome', pval_col = 'P', beta_col = 'beta',
                    snp_col ='variant_id', effect_allele_col = 'reference_allele', other_allele_col = 'other_allele',
                    se_col = 'OR_se', eaf_col = 'eaf', samplesize_col = 'n_samples', chr_col = 'chromosome',
                    pos_col = 'base_pair_location')
protein.LADA <- harmonise_data(protein.data, LADA)
protein.LADA <- protein.LADA[protein.LADA$mr_keep==T, ]
protein.LADA$chr.outcome <- as.numeric(protein.LADA$chr.outcome)
protein.LADA$pos.outcome <- as.numeric(protein.LADA$pos.outcome)
protein.LADA$pval.exposure <- as.numeric(protein.LADA$pval.exposure)
protein <- data.frame(chromosome=protein.LADA$chr.outcome, position=protein.LADA$pos.outcome, 
                      rsid=protein.LADA$SNP, effect_allele=protein.LADA$effect_allele.exposure, 
                      other_allele=protein.LADA$other_allele.exposure,
                      effect=protein.LADA$beta.exposure, std_err=protein.LADA$se.exposure, 
                      pval=protein.LADA$pval.exposure)
p1 <- gg_locusplot(df=protein, lead_snp = 'rs11024589', chrom=chromosome, pos=position,
                   ref=effect_allele, alt=other_allele, p_value = pval, genome_build = 'GRCh37', population='EUR',
                   plot_genes = F, plot_pvalue_threshold = 1);p1
p1 <- p1 + 
  labs(x = "Chromosome 11 (Mb)", fill=parse(text = "LD  (r^2)")) +
  theme(panel.grid=element_blank(), legend.background = element_rect(fill = 'white', colour = 'black'),
        legend.position = c(0.14, 0.99));p1
Ipaper::write_fig(p1,
                  file='regional_plot/cis-pQTL for SAA1.jpg',
                  width=10,
                  height=4,
                  devices='jpg',
                  res=300,
                  show=T)
### LADA.data
LADA.data.SAA1 <- data.frame(chromosome=protein.LADA$chr.outcome, position=protein.LADA$pos.outcome, 
                             rsid=protein.LADA$SNP, effect_allele=protein.LADA$effect_allele.outcome, 
                             other_allele=protein.LADA$other_allele.outcome,
                             effect=protein.LADA$beta.outcome, std_err=protein.LADA$se.outcome, 
                             pval=protein.LADA$pval.outcome)
source('D:/Others/ACLY与cancer/共定位/gg_locusplot_self.r')
p2 <- gg_locusplot_self(df=LADA.data.SAA1, lead_snp = 'rs11024589', chrom=chromosome, pos=position,
                        ref=effect_allele, alt=other_allele, p_value = pval, genome_build = 'GRCh37', population='EUR',
                        plot_genes = T, plot_pvalue_threshold = 1);p2
Ipaper::write_fig(p2,
                  file='regional_plot/SAA1_LADA.jpg',
                  width=10,
                  height=4,
                  devices='jpg',
                  res=300,
                  show=T)

#### SAA2 ####
filename <- 'SAA2_coloc_SNPs.txt'
input_file <- paste0('colocalization_top_pqtl_500kb/', filename)
fp <- read.table(input_file, sep='\t', header = T, stringsAsFactors = F)
protein.data <- format_data(fp, snp_col = 'rsids', effect_allele_col = 'effectAllele',
                            other_allele_col = 'otherAllele', beta_col = 'Beta',
                            se_col = 'SE', pval_col = 'Pval',
                            samplesize_col = 'N')
LADA <- LADA.data[LADA.data$variant_id%in%protein.data$SNP, ]
LADA <- format_data(LADA, type='outcome', pval_col = 'P', beta_col = 'beta',
                    snp_col ='variant_id', effect_allele_col = 'reference_allele', other_allele_col = 'other_allele',
                    se_col = 'OR_se', eaf_col = 'eaf', samplesize_col = 'n_samples', chr_col = 'chromosome',
                    pos_col = 'base_pair_location')
protein.LADA <- harmonise_data(protein.data, LADA)
protein.LADA <- protein.LADA[protein.LADA$mr_keep==T, ]
protein.LADA$chr.outcome <- as.numeric(protein.LADA$chr.outcome)
protein.LADA$pos.outcome <- as.numeric(protein.LADA$pos.outcome)
protein.LADA$pval.exposure <- as.numeric(protein.LADA$pval.exposure)
protein <- data.frame(chromosome=protein.LADA$chr.outcome, position=protein.LADA$pos.outcome, 
                      rsid=protein.LADA$SNP, effect_allele=protein.LADA$effect_allele.exposure, 
                      other_allele=protein.LADA$other_allele.exposure,
                      effect=protein.LADA$beta.exposure, std_err=protein.LADA$se.exposure, 
                      pval=protein.LADA$pval.exposure)
p1 <- gg_locusplot(df=protein, lead_snp = 'rs11024589', chrom=chromosome, pos=position,
                   ref=effect_allele, alt=other_allele, p_value = pval, genome_build = 'GRCh37', population='EUR',
                   plot_genes = F, plot_pvalue_threshold = 1);p1
p1 <- p1 + 
  labs(x = "Chromosome 11 (Mb)", fill=parse(text = "LD  (r^2)")) +
  theme(panel.grid=element_blank(), legend.background = element_rect(fill = 'white', colour = 'black'),
        legend.position = c(0.14, 0.99));p1
Ipaper::write_fig(p1,
                  file='regional_plot/cis-pQTL for SAA2.jpg',
                  width=10,
                  height=4,
                  devices='jpg',
                  res=300,
                  show=T)
### LADA.data
LADA.data.SAA2 <- data.frame(chromosome=protein.LADA$chr.outcome, position=protein.LADA$pos.outcome, 
                             rsid=protein.LADA$SNP, effect_allele=protein.LADA$effect_allele.outcome, 
                             other_allele=protein.LADA$other_allele.outcome,
                             effect=protein.LADA$beta.outcome, std_err=protein.LADA$se.outcome, 
                             pval=protein.LADA$pval.outcome)
source('D:/Others/ACLY与cancer/共定位/gg_locusplot_self.r')
p2 <- gg_locusplot_self(df=LADA.data.SAA2, lead_snp = 'rs11024589', chrom=chromosome, pos=position,
                        ref=effect_allele, alt=other_allele, p_value = pval, genome_build = 'GRCh37', population='EUR',
                        plot_genes = T, plot_pvalue_threshold = 1);p2
Ipaper::write_fig(p2,
                  file='regional_plot/SAA2_LADA.jpg',
                  width=10,
                  height=4,
                  devices='jpg',
                  res=300,
                  show=T)

#### MIA ####
filename <- 'MIA_coloc_SNPs.txt'
input_file <- paste0('colocalization_top_pqtl_500kb/', filename)
fp <- read.table(input_file, sep='\t', header = T, stringsAsFactors = F)
protein.data <- format_data(fp, snp_col = 'rsids', effect_allele_col = 'effectAllele',
                            other_allele_col = 'otherAllele', beta_col = 'Beta',
                            se_col = 'SE', pval_col = 'Pval',
                            samplesize_col = 'N')
LADA <- LADA.data[LADA.data$variant_id%in%protein.data$SNP, ]
LADA <- format_data(LADA, type='outcome', pval_col = 'P', beta_col = 'beta',
                    snp_col ='variant_id', effect_allele_col = 'reference_allele', other_allele_col = 'other_allele',
                    se_col = 'OR_se', eaf_col = 'eaf', samplesize_col = 'n_samples', chr_col = 'chromosome',
                    pos_col = 'base_pair_location')
protein.LADA <- harmonise_data(protein.data, LADA)
protein.LADA <- protein.LADA[protein.LADA$mr_keep==T, ]
protein.LADA$chr.outcome <- as.numeric(protein.LADA$chr.outcome)
protein.LADA$pos.outcome <- as.numeric(protein.LADA$pos.outcome)
protein.LADA$pval.exposure <- as.numeric(protein.LADA$pval.exposure)
protein <- data.frame(chromosome=protein.LADA$chr.outcome, position=protein.LADA$pos.outcome, 
                      rsid=protein.LADA$SNP, effect_allele=protein.LADA$effect_allele.exposure, 
                      other_allele=protein.LADA$other_allele.exposure,
                      effect=protein.LADA$beta.exposure, std_err=protein.LADA$se.exposure, 
                      pval=protein.LADA$pval.exposure)
p1 <- gg_locusplot(df=protein, lead_snp = 'rs2279011', chrom=chromosome, pos=position,
                   ref=effect_allele, alt=other_allele, p_value = pval, genome_build = 'GRCh37', population='EUR',
                   plot_genes = F, plot_pvalue_threshold = 1);p1
p1 <- p1 + 
  labs(x = "Chromosome 19 (Mb)", fill=parse(text = "LD  (r^2)")) +
  theme(panel.grid=element_blank(), legend.background = element_rect(fill = 'white', colour = 'black'),
        legend.position = c(0.14, 0.99));p1
Ipaper::write_fig(p1,
                  file='regional_plot/cis-pQTL for MIA.jpg',
                  width=10,
                  height=4,
                  devices='jpg',
                  res=300,
                  show=T)
### LADA.data
LADA.data.MIA <- data.frame(chromosome=protein.LADA$chr.outcome, position=protein.LADA$pos.outcome, 
                            rsid=protein.LADA$SNP, effect_allele=protein.LADA$effect_allele.outcome, 
                            other_allele=protein.LADA$other_allele.outcome,
                            effect=protein.LADA$beta.outcome, std_err=protein.LADA$se.outcome, 
                            pval=protein.LADA$pval.outcome)
source('D:/Others/ACLY与cancer/共定位/gg_locusplot_self.r')
p2 <- gg_locusplot_self(df=LADA.data.MIA, lead_snp = 'rs2279011', chrom=chromosome, pos=position,
                        ref=effect_allele, alt=other_allele, p_value = pval, genome_build = 'GRCh37', population='EUR',
                        plot_genes = T, plot_pvalue_threshold = 1);p2
Ipaper::write_fig(p2,
                  file='regional_plot/MIA_LADA.jpg',
                  width=10,
                  height=4,
                  devices='jpg',
                  res=300,
                  show=T)