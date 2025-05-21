BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")  # For human gene annotations (change as needed)
library(clusterProfiler)
library(org.Hs.eg.db)  # Human gene annotation database
library(ggplot2)

# risk factors gene list of LADA
risk_factors <- c('CXCL10')

# Background gene list (the 1,389 proteins)
background_proteins <- read.csv('background_gene_list.csv', sep=',', stringsAsFactors = F, header = T)
background_proteins <- background_proteins$Genes

# Map protein names to Entrez IDs
risk_factors_entrez <- bitr(risk_factors, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
background_proteins_entrez <- bitr(background_proteins, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
background_proteins_unmapped_genes <- background_proteins[!background_proteins%in%background_proteins_entrez$SYMBOL]
cat(background_proteins_unmapped_genes)

# Extract the Entrez IDs
risk_factors_entrez_ids <- risk_factors_entrez$ENTREZID
background_proteins_entrez_ids <- background_proteins_entrez$ENTREZID

# Perform GO enrichment analysis
# Biological Process
go_enrichment_BP <- enrichGO(
  gene = risk_factors_entrez_ids,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  universe = background_proteins_entrez_ids,
  ont = "BP",  # Biological Process (you can choose "MF" for Molecular Function or "CC" for Cellular Component)
  pAdjustMethod = "BH",  # Benjamini-Hochberg correction
  qvalueCutoff = 0.05,  # Set significance threshold
  readable = TRUE
)
summary(go_enrichment_BP)

# Molecular Function
go_enrichment_MF <- enrichGO(
  gene = risk_factors_entrez_ids,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  universe = background_proteins_entrez_ids,
  ont = "MF",  # Biological Process (you can choose "MF" for Molecular Function or "CC" for Cellular Component)
  pAdjustMethod = "BH",  # Benjamini-Hochberg correction
  qvalueCutoff = 0.05,  # Set significance threshold
  readable = TRUE
)
summary(go_enrichment_MF)

# Cellular Component
go_enrichment_CC <- enrichGO(
  gene = risk_factors_entrez_ids,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  universe = background_proteins_entrez_ids,
  ont = "CC",  # Biological Process (you can choose "MF" for Molecular Function or "CC" for Cellular Component)
  qvalueCutoff = 0.05,  # Set significance threshold
  readable = TRUE
)
summary(go_enrichment_CC)

# KEGG
kegg_enrichment <- enrichKEGG(
  gene = risk_factors_entrez_ids,
  organism = "hsa",
  universe = background_proteins_entrez_ids,
  pAdjustMethod = "fdr",  # Benjamini-Hochberg correction
  qvalueCutoff = 0.05
)
summary(kegg_enrichment)

write.table(kegg_enrichment, 'KEGG_enrichment.txt', sep='\t', row.names = F, quote = F)

# Convert the enrichment result to a data frame
kegg_df <- as.data.frame(kegg_enrichment)
# Create a customized bar plot
ggplot(kegg_df, aes(x = reorder(Description, -p.adjust), y = Count, fill = p.adjust)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(x = "KEGG Pathway", y = "Gene Count", title = "Top 10 Enriched KEGG Pathways") +
  scale_fill_gradient(low = "blue", high = "red") +
  theme_minimal()
