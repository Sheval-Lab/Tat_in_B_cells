## Load libraries --------------------------------------------------------------
library(tidyverse)
library(clusterProfiler)

input_dir <- "output/tables/01_DGEA/deseq"
output_dir <- "output/tables/03_GE/KEGG"
fig_dir <- "output/figures/03_GE/KEGG"


## Load DGEA results -----------------------------------------------------------
tat_vs_lcl <- read_tsv(str_c(input_dir, "Tat_vs_LCL.tsv", sep = "/"))
tat_vs_lcl_LFC <- read_tsv(str_c(input_dir, "Tat_vs_LCL.LFC.tsv", sep = "/"))


## Prepare input data for GSEA -------------------------------------------------
### Gene IDs conversion
geneIDs_table <- bitr(tat_vs_lcl$ensID, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db") %>% 
  filter(!duplicated(ENTREZID), !duplicated(ENSEMBL)) %>% 
  dplyr::rename(ensID = ENSEMBL, entrezID = ENTREZID)

tat <- tat_vs_lcl %>% 
  dplyr::select(ensID, gene_name, gene_type, baseMean, stat, padj, log2FC) %>% 
  ## Add shrunken log2FC values
  left_join(., dplyr::select(tat_vs_lcl_LFC, ensID, log2FC_LFC = log2FC), by = "ensID") %>% 
  ## Add entrezID
  inner_join(., geneIDs_table, by = "ensID")


### Prepare pre-ranked gene list
gene_list <- tat$stat
names(gene_list) = as.character(tat$entrezID)
gene_list = sort(gene_list, decreasing = TRUE)

### Prepare gene list of log2FC values
fc_list <- tat %>% filter(padj < 0.05) %>% pull(log2FC_LFC)
names(fc_list) <- tat %>% filter(padj < 0.05) %>% pull(entrezID) %>% as.character()
fc_list <- sort(fc_list, decreasing = TRUE)

fc_list_flt <- tat %>% filter(baseMean >= 100, padj < 0.05) %>% pull(log2FC_LFC)
names(fc_list_flt) <- tat %>% filter(baseMean >= 100, padj < 0.05) %>% pull(entrezID) %>% as.character()
fc_list_flt <- sort(fc_list_flt, decreasing = TRUE)

### Prepare gene lists of UP and DOWN regulated DEGs ---------------------------
genes_up <- tat %>% filter(padj < 0.05, log2FC >= 1) %>% pull(entrezID)
genes_dn <- tat %>% filter(padj < 0.05, log2FC <= (-1)) %>% pull(entrezID)


## Run ORA against KEGG database ----------------------------------------------
### ORA of UP-regulated DEGs
kegg_ora_up <- enrichKEGG(
  gene = genes_up,
  organism = "hsa",
  pvalueCutoff = 0.05,
  universe = tat$entrezID)

head(kegg_ora_up)

### ORA of DOWN-regulated DEGs
kegg_ora_dn <- enrichKEGG(
  gene = genes_dn,
  organism = "hsa",
  pvalueCutoff = 0.05,
  universe = tat$entrezID)

head(kegg_ora_dn)

### Combine 2 results to plot them together
kegg_ora <- kegg_ora_up
kegg_ora@gene <- c(kegg_ora_up@gene, kegg_ora_dn@gene)
kegg_ora@result <- rbind(head(kegg_ora_up), head(kegg_ora_dn))
# rownames(kegg_ora@result) <- kegg_ora@result$ID
head(kegg_ora)


### Make cnetplot
cnetplot(setReadable(kegg_ora, 'org.Hs.eg.db', 'ENTREZID'), 
         foldChange = fc_list, #node_label = "none",
         cex_category = 0.7,
         cex_label_gene = 0.5,
         cex_label_category = 0.5) +
  scale_colour_gradient2(name = "log2FoldChange", low = "blue", mid = "white", high = "red") + 
  guides(color = guide_colourbar(barwidth = 0.7, barheight = 4))

ggsave(str_c(fig_dir, "KEGG_ORA_cnet.png", sep = "/"), units = "cm", width = 16, height = 8)
ggsave(str_c(fig_dir, "KEGG_ORA_cnet.svg", sep = "/"), units = "cm", width = 16, height = 8)
  

## Run GSEA against KEGG database ----------------------------------------------
set.seed(42) ## use seed for reproducibility
kegg_gsea <- gseKEGG(
  geneList = gene_list, 
  organism = "hsa", 
  minGSSize = 110, 
  pvalueCutoff = 0.05, 
  verbose = FALSE,
  seed = TRUE)

head(kegg_gsea)

### Make dotplot of all enriched categories
dotplot(kegg_gsea, showCategory = 11, split=".sign") + 
  facet_grid(.~.sign) +
  theme(
    strip.background = element_rect(fill = "white", colour = "black", size = 0.5, linetype = "solid"),
    strip.text = element_text(face = "bold", size = 12))

ggsave(str_c(fig_dir, "KEGG_GSEA_dot.png", sep = "/"), units = "cm", width = 6, height = 4, scale = 3)
ggsave(str_c(fig_dir, "KEGG_GSEA_dot.svg", sep = "/"), units = "cm", width = 6, height = 4, scale = 3)


### Modify gseaResult object: 
### keep only DEGs in core_enrichment,
### re-calculate size of category nodes (for plotting)
kegg_gsea_mod <- kegg_gsea
kegg_gsea_mod@result$core_enrichment <- 
  kegg_gsea_mod@result$core_enrichment %>% 
  str_split("/") %>% 
  map(function(x) x[x %in% names(fc_list_flt)]) %>% 
  map_chr(str_c, collapse = "/") 

kegg_gsea_mod <- setReadable(kegg_gsea_mod, "org.Hs.eg.db", "ENTREZID")

### Make cnetplot of all enriched categories
cnetplot(kegg_gsea_mod, showCategory = 11, foldChange = fc_list_flt,
         cex_category = 0.7,
         cex_label_gene = 0.5,
         cex_label_category = 0.5) +
  scale_colour_gradient2(name = "log2FoldChange", low = "blue", mid = "white", high = "red") + 
  guides(color = guide_colourbar(barwidth = 0.7, barheight = 4))

ggsave(str_c(fig_dir, "KEGG_GSEA_cnet.png", sep = "/"), units = "cm", width = 16, height = 10, scale = 1.5)
ggsave(str_c(fig_dir, "KEGG_GSEA_cnet.svg", sep = "/"), units = "cm", width = 16, height = 10, scale = 1.5)


### Make cnetplot of pre-filtered enriched categories
cat2keep <- c("hsa04060", "hsa05168", "hsa04621", "hsa05205", "hsa05165", "hsa04510", "hsa04120")
kegg_gsea_mod_flt <- kegg_gsea_mod
kegg_gsea_mod_flt@result <- kegg_gsea_mod_flt@result %>% 
  filter(ID %in% cat2keep)

cnetplot(kegg_gsea_mod_flt, showCategory = 11, foldChange = fc_list_flt,
         cex_category = 0.7,
         cex_label_gene = 0.5,
         cex_label_category = 0.5) +
  scale_colour_gradient2(name = "log2FoldChange", low = "blue", mid = "white", high = "red") + 
  guides(color = guide_colourbar(barwidth = 0.7, barheight = 4))

ggsave(str_c(fig_dir, "KEGG_GSEA_cnet_flt.png", sep = "/"), units = "cm", width = 16, height = 10, scale = 1.5)
ggsave(str_c(fig_dir, "KEGG_GSEA_cnet_flt.svg", sep = "/"), units = "cm", width = 16, height = 10, scale = 1.5)

