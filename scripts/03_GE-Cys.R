## Load libraries --------------------------------------------------------------
library(tidyverse)
library(scales)
library(clusterProfiler)
library(ggraph)
library(patchwork)
source("scripts/functions/exclude_nested_clusters.R")

input_dir <- "output/tables/01_DGEA/deseq"
output_dir <- "output/tables/03_GE/Cys"
fig_dir <- "output/figures/03_GE/Cys"


## Specify log2FoldChange threshold for DEGs -----------------------------------
log2fc_threshold <- log2(1.5)


## Load DGEA results -----------------------------------------------------------
tat_vs_cys <- read_tsv(str_c(input_dir, "Tat_vs_Cys.tsv", sep = "/"))
tat_vs_cys_LFC <- read_tsv(str_c(input_dir, "Tat_vs_Cys.LFC.tsv", sep = "/"))


## Prepare input data for GSEA -------------------------------------------------
### Gene IDs conversion
geneIDs_table <- bitr(tat_vs_cys$ensID, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db") %>% 
  filter(!duplicated(ENTREZID), !duplicated(ENSEMBL)) %>% 
  dplyr::rename(ensID = ENSEMBL, entrezID = ENTREZID)

cys <- tat_vs_cys %>% 
  dplyr::select(ensID, gene_name, gene_type, baseMean, stat, padj, log2FC) %>% 
  ## Add shrunken log2FC values
  left_join(., dplyr::select(tat_vs_cys_LFC, ensID, log2FC_LFC = log2FC), by = "ensID") %>% 
  ## Add entrezID
  inner_join(., geneIDs_table, by = "ensID")


### Prepare gene list of log2FC values
fc_list <- cys %>% filter(padj < 0.05) %>% mutate(log2FC_LFC = -log2FC_LFC) %>% pull(log2FC_LFC)
names(fc_list) <- cys %>% filter(padj < 0.05) %>% pull(entrezID) %>% as.character()
fc_list <- sort(fc_list, decreasing = TRUE)

### Prepare gene lists of UP and DOWN regulated DEGs ---------------------------
### Fold change >= 1.5
### Take a reverse comparison: Cys vs Tat
genes_up <- cys %>% filter(padj < 0.05, log2FC <= (-log2fc_threshold)) %>% pull(entrezID)
genes_dn <- cys %>% filter(padj < 0.05, log2FC >= log2fc_threshold) %>% pull(entrezID)


## Run ORA against KEGG database ----------------------------------------------
### ORA of UP-regulated DEGs
kegg_ora_up <- enrichKEGG(
  gene = genes_up,
  organism = "hsa",
  pvalueCutoff = 0.05,
  universe = cys$entrezID)

head(kegg_ora_up)

### ORA of DOWN-regulated DEGs
kegg_ora_dn <- enrichKEGG(
  gene = genes_dn,
  organism = "hsa",
  pvalueCutoff = 0.05,
  universe = cys$entrezID)

head(kegg_ora_dn)

### Combine 2 results to plot them together
kegg_ora <- kegg_ora_up
kegg_ora@gene <- c(kegg_ora_up@gene, kegg_ora_dn@gene)
kegg_ora@result <- rbind(head(kegg_ora_up), head(kegg_ora_dn))
rownames(kegg_ora@result) <- kegg_ora@result$ID

#### Add UP/DOWN labels to enriched categories
kegg_ora %<>% 
  mutate(
    sign = if_else(ID %in% head(kegg_ora_up)$ID, "activated", "suppressed"),
    FoldEnrichment = DOSE::parse_ratio(GeneRatio) / DOSE::parse_ratio(BgRatio)) %>% 
  arrange(desc(sign), FoldEnrichment)

head(kegg_ora)


### Make dotplot of all enriched categories
dotpl_ora <- dotplot(kegg_ora, x = "FoldEnrichment", showCategory = nrow(kegg_ora), split = "sign") + 
  facet_grid(.~sign) +
  theme(
    strip.background = element_rect(fill = "white", colour = "black", size = 0.5, linetype = "solid"),
    strip.text = element_text(face = "bold", size = 12)) +
  scale_y_discrete(limits = kegg_ora@result$Description)


dotpl_ora + ggh4x::force_panelsizes(cols = unit(3.5, "cm")) + xlim(4,23)
ggsave(str_c(fig_dir, "KEGG_ORA_dot.png", sep = "/"), units = "cm", width = 7, height = 3, scale = 3)
ggsave(str_c(fig_dir, "KEGG_ORA_dot.pdf", sep = "/"), units = "cm", width = 7, height = 3, scale = 3)


### Make cnetplot
cnetplot(setReadable(kegg_ora, 'org.Hs.eg.db', 'ENTREZID'), 
         showCategory = nrow(kegg_ora),
         foldChange = fc_list, 
         # cex_category = 0.7,
         cex_label_gene = 0.8,
         # cex_label_category = 0.5,
         node_label = "gene") +
  scale_colour_gradient(name = "log2FoldChange", low = "royalblue4", high = "skyblue2") + 
  guides(color = guide_colourbar(barwidth = 0.7, barheight = 4))

ggsave(str_c(fig_dir, "KEGG_ORA_cnet.png", sep = "/"), units = "cm", width = 5, height = 3, scale = 2.8)
ggsave(str_c(fig_dir, "KEGG_ORA_cnet.pdf", sep = "/"), units = "cm", width = 5, height = 3, scale = 2.8)

