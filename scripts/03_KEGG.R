## Load libraries --------------------------------------------------------------
library(tidyverse)
library(scales)
library(clusterProfiler)
library(ggraph)
library(patchwork)
source("scripts/functions/exclude_nested_clusters.R")

input_dir <- "output/tables/01_DGEA/deseq"
output_dir <- "output/tables/03_GE/KEGG"
fig_dir <- "output/figures/03_GE/KEGG"


## Specify log2FoldChange threshold for DEGs -----------------------------------
log2fc_threshold <- log2(1.5)


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
### Fold change >= 1.5
genes_up <- tat %>% filter(padj < 0.05, log2FC >= log2fc_threshold) %>% pull(entrezID)
genes_dn <- tat %>% filter(padj < 0.05, log2FC <= (-log2fc_threshold)) %>% pull(entrezID)


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
rownames(kegg_ora@result) <- kegg_ora@result$ID

#### Add UP/DOWN labels to enriched categories
kegg_ora %<>% 
  mutate(
    sign = if_else(ID %in% head(kegg_ora_up)$ID, "activated", "suppressed"),
    FoldEnrichment = DOSE::parse_ratio(GeneRatio) / DOSE::parse_ratio(BgRatio)) %>% 
  arrange(desc(sign), FoldEnrichment)

head(kegg_ora)

### Save as RDS file
saveRDS(kegg_ora, str_c(output_dir, "Tat_kegg_ora.rds", sep = "/"))

kegg_ora <- readRDS(str_c(output_dir, "Tat_kegg_ora.rds", sep = "/"))


### Make dotplot of all enriched categories
dotpl_ora <- dotplot(kegg_ora, x = "FoldEnrichment", showCategory = nrow(kegg_ora), split = "sign") + 
  facet_grid(.~sign) +
  theme(
    strip.background = element_rect(fill = "white", colour = "black", size = 0.5, linetype = "solid"),
    strip.text = element_text(face = "bold", size = 12)) +
  scale_y_discrete(limits = kegg_ora@result$Description)
  

dotpl_ora + ggh4x::force_panelsizes(cols = unit(3.5, "cm"))
ggsave(str_c(fig_dir, "KEGG_ORA_dot.png", sep = "/"), units = "cm", width = 7, height = 4, scale = 3.4)
ggsave(str_c(fig_dir, "KEGG_ORA_dot.svg", sep = "/"), units = "cm", width = 7, height = 4, scale = 3.4)


### Make cnetplot
cnetplot(setReadable(kegg_ora, 'org.Hs.eg.db', 'ENTREZID'), 
         showCategory = nrow(kegg_ora),
         foldChange = fc_list, #node_label = "none",
         cex_category = 0.7,
         cex_label_gene = 0.5,
         cex_label_category = 0.5) +
  scale_colour_gradient2(name = "log2FoldChange", low = "blue", mid = "white", high = "red") + 
  guides(color = guide_colourbar(barwidth = 0.7, barheight = 4))

ggsave(str_c(fig_dir, "KEGG_ORA_cnet.png", sep = "/"), units = "cm", width = 16, height = 10, scale = 1.5)
ggsave(str_c(fig_dir, "KEGG_ORA_cnet.svg", sep = "/"), units = "cm", width = 16, height = 10, scale = 1.5)


### Exclude nested clusters
kegg_ora_flt_1 <- exclude_nested_clusters(kegg_ora, 0)

#### Make cnetplot on filtered results
cnetplot(setReadable(kegg_ora_flt_1, 'org.Hs.eg.db', 'ENTREZID'), 
         showCategory = nrow(kegg_ora),
         foldChange = fc_list, #node_label = "none",
         cex_category = 0.7,
         cex_label_gene = 0.5,
         cex_label_category = 0.5) +
  scale_colour_gradient2(name = "log2FoldChange", low = "blue", mid = "white", high = "red") + 
  guides(color = guide_colourbar(barwidth = 0.7, barheight = 4))

ggsave(str_c(fig_dir, "KEGG_ORA_cnet_flt_1.png", sep = "/"), units = "cm", width = 16, height = 9, scale = 1.5)
ggsave(str_c(fig_dir, "KEGG_ORA_cnet_flt_1.svg", sep = "/"), units = "cm", width = 16, height = 9, scale = 1.5)


#### Filter even more
clusters2exclude <- c("hsa05323", "hsa04640", "hsa04061", "hsa05164")

kegg_ora_flt_2 <- exclude_nested_clusters(kegg_ora, 2)
## or manually:
# kegg_ora_flt_2 <- kegg_ora %>% filter(!(ID %in% clusters2exclude))

cnet_ora <- cnetplot(setReadable(kegg_ora_flt_2, 'org.Hs.eg.db', 'ENTREZID'), 
         showCategory = nrow(kegg_ora),
         foldChange = fc_list, #node_label = "none",
         cex_category = 0.7,
         cex_label_gene = 0.5,
         cex_label_category = 0.5) +
  scale_colour_gradient2(name = "log2FoldChange", low = "blue", mid = "white", high = "red") + 
  guides(color = guide_colourbar(barwidth = 0.7, barheight = 4))

ggsave(str_c(fig_dir, "KEGG_ORA_cnet_flt.png", sep = "/"), cnet_ora, units = "cm", width = 16, height = 8, scale = 1.5)
ggsave(str_c(fig_dir, "KEGG_ORA_cnet_flt.svg", sep = "/"), cnet_ora, units = "cm", width = 16, height = 8, scale = 1.5)



### Make plots for publishing --------------------------------------------------
#### Dotplots ------------------------------------------------------------------
dotpl_up <- dotplot(filter(kegg_ora, sign == "activated"), x = "FoldEnrichment", showCategory = nrow(kegg_ora), split = "sign") + 
  facet_grid(.~sign) +
  theme(
    strip.background = element_rect(fill = "white", colour = "black", size = 0.5, linetype = "solid"),
    strip.text = element_text(face = "bold", size = 12)) +
  scale_x_continuous(limits = c(2, 4.2), breaks = c(2, 3, 4)) +
  scale_y_discrete(limits = filter(kegg_ora, sign == "activated")@result$Description) +
  ggh4x::force_panelsizes(cols = unit(3, "cm"), rows = unit(5, "cm")) +
  guides(
    size = guide_legend(order = 1),
    color = guide_colourbar(order = 2, reverse = TRUE))


dotpl_dn <- dotplot(filter(kegg_ora, sign == "suppressed"), x = "FoldEnrichment", showCategory = nrow(kegg_ora), split = "sign") + 
  facet_grid(.~sign) +
  theme(
    strip.background = element_rect(fill = "white", colour = "black", size = 0.5, linetype = "solid"),
    strip.text = element_text(face = "bold", size = 12)) +
  scale_x_continuous(limits = c(2.9, 6.9), breaks = c(3, 5, 7)) +
  scale_y_discrete(limits = filter(kegg_ora, sign == "suppressed")@result$Description) +
  ggh4x::force_panelsizes(cols = unit(3, "cm"), rows = unit(5, "cm")) +
  guides(
    size = guide_legend(order = 1),
    color = guide_colourbar(order = 2, reverse = TRUE))

#### Edit legends
dt_size_limits <- c(min(c(dotpl_up$data$Count, dotpl_dn$data$Count)), max(c(dotpl_up$data$Count, dotpl_dn$data$Count)))
dt_fill_limits <- c(min(c(dotpl_up$data$p.adjust, dotpl_dn$data$p.adjust)), max(c(dotpl_up$data$p.adjust, dotpl_dn$data$p.adjust)))

dotpl_up <- dotpl_up +
  scale_size_continuous(limits = dt_size_limits, breaks = c(10, 15, 20)) +
  scale_color_continuous(limits = dt_fill_limits, breaks = c(0.005, 0.015, 0.025),
                         low = "red", high = "blue") 
ggsave(str_c(fig_dir, "KEGG_ORA_dot_up.png", sep = "/"), dotpl_up, units = "cm", width = 8, height = 6, scale = 2)
ggsave(str_c(fig_dir, "KEGG_ORA_dot_up.pdf", sep = "/"), dotpl_up, units = "cm", width = 8, height = 6, scale = 2)


dotpl_dn <- dotpl_dn +
  scale_size_continuous(limits = dt_size_limits, breaks = c(10, 15, 20)) +
  scale_color_continuous(limits = dt_fill_limits, breaks = c(0.005, 0.015, 0.025), 
                         low = "red", high = "blue") 
ggsave(str_c(fig_dir, "KEGG_ORA_dot_dn.png", sep = "/"), dotpl_dn, units = "cm", width = 10, height = 6, scale = 2)
ggsave(str_c(fig_dir, "KEGG_ORA_dot_dn.pdf", sep = "/"), dotpl_dn, units = "cm", width = 10, height = 6, scale = 2)


#### Cnetplots -----------------------------------------------------------------
cnet_up <- cnetplot(filter(setReadable(kegg_ora_flt_2, 'org.Hs.eg.db', 'ENTREZID'), sign == "activated"), 
                    showCategory = nrow(kegg_ora),
                    foldChange = fc_list, 
                    # cex_category = 1,
                    cex_label_gene = 0.8,
                    # cex_label_category = 0.9, 
                    node_label = "gene") +
  scale_colour_gradient(name = "log2FoldChange", low = "lightcoral", high = "red4") + 
  guides(color = guide_colourbar(barwidth = 0.7, barheight = 10)) 

cnet_dn <- cnetplot(filter(setReadable(kegg_ora_flt_2, 'org.Hs.eg.db', 'ENTREZID'), sign == "suppressed"), 
                    showCategory = nrow(kegg_ora),
                    foldChange = fc_list, 
                    # cex_category = 1,
                    cex_label_gene = 0.8,
                    # cex_label_category = 0.9, 
                    node_label = "gene") +
  scale_colour_gradient(name = "log2FoldChange", low = "royalblue4", high = "skyblue2") + 
  guides(
    color = guide_colourbar(barwidth = 0.7, barheight = 4, order = 1))


ggsave(str_c(fig_dir, "KEGG_ORA_cnet_up.png", sep = "/"), cnet_up, units = "cm", width = 5, height = 5, scale = 3.5)
ggsave(str_c(fig_dir, "KEGG_ORA_cnet_up.pdf", sep = "/"), cnet_up, units = "cm", width = 5, height = 5, scale = 3.5)

ggsave(str_c(fig_dir, "KEGG_ORA_cnet_dn.png", sep = "/"), cnet_dn, units = "cm", width = 5, height = 5, scale = 3.5)
ggsave(str_c(fig_dir, "KEGG_ORA_cnet_dn.pdf", sep = "/"), cnet_dn, units = "cm", width = 5, height = 5, scale = 3.5)




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
dotplot(kegg_gsea, showCategory = nrow(kegg_gsea), split = ".sign") + 
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
cnetplot(kegg_gsea_mod, showCategory = nrow(kegg_gsea_mod), foldChange = fc_list_flt,
         cex_category = 0.7,
         cex_label_gene = 0.5,
         cex_label_category = 0.5) +
  scale_colour_gradient2(name = "log2FoldChange", low = "blue", mid = "white", high = "red") + 
  guides(color = guide_colourbar(barwidth = 0.7, barheight = 4))

ggsave(str_c(fig_dir, "KEGG_GSEA_cnet.png", sep = "/"), units = "cm", width = 16, height = 10, scale = 1.5)
ggsave(str_c(fig_dir, "KEGG_GSEA_cnet.svg", sep = "/"), units = "cm", width = 16, height = 10, scale = 1.5)


### Make cnetplot of pre-filtered enriched categories
clusters2exclude <- c("hsa05161", "hsa05160", "hsa05164", "hsa05171")

kegg_gsea_mod_flt <- kegg_gsea_mod
kegg_gsea_mod_flt@result <- kegg_gsea_mod_flt@result %>% 
  filter(!(ID %in% clusters2exclude))

cnetplot(kegg_gsea_mod_flt, showCategory = nrow(kegg_gsea_mod_flt), foldChange = fc_list_flt,
         cex_category = 0.7,
         cex_label_gene = 0.5,
         cex_label_category = 0.5) +
  scale_colour_gradient2(name = "log2FoldChange", low = "blue", mid = "white", high = "red") + 
  guides(color = guide_colourbar(barwidth = 0.7, barheight = 4))

ggsave(str_c(fig_dir, "KEGG_GSEA_cnet_flt.png", sep = "/"), units = "cm", width = 16, height = 10, scale = 1.5)
ggsave(str_c(fig_dir, "KEGG_GSEA_cnet_flt.svg", sep = "/"), units = "cm", width = 16, height = 10, scale = 1.5)

