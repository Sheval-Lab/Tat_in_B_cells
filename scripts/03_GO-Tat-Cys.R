## Load libraries --------------------------------------------------------------
library(tidyverse)
library(clusterProfiler)

input_dir <- "output/tables/01_DGEA/deseq"
output_dir <- "output/tables/03_GE/Tat-Cys"
fig_dir <- "output/figures/03_GE/Tat-Cys"


## Specify log2FoldChange threshold for DEGs -----------------------------------
log2fc_threshold <- log2(1.5)


## Load DGEA results -----------------------------------------------------------
tat_vs_lcl <- read_tsv(str_c(input_dir, "Tat_vs_LCL.LFC.DE.tsv", sep = "/"))
cys_vs_lcl <- read_tsv(str_c(input_dir, "Cys_vs_LCL.LFC.DE.tsv", sep = "/"))
all_counts <- read_tsv("output/tables/01_DGEA/processed_counts/counts_raw.tsv")


## Genes -----------------------------------------------------------------------
### All genes
tat_cys <- list(
  tat = filter(tat_vs_lcl, abs(log2FC) >= log2fc_threshold) %>%  pull(gene_name), 
  cys = filter(cys_vs_lcl, abs(log2FC) >= log2fc_threshold) %>%  pull(gene_name),
  tat_up = filter(tat_vs_lcl, log2FC >= log2fc_threshold) %>%  pull(gene_name), 
  cys_up = filter(cys_vs_lcl, log2FC >= log2fc_threshold) %>%  pull(gene_name), 
  tat_dn = filter(tat_vs_lcl, log2FC <= (-log2fc_threshold)) %>%  pull(gene_name),
  cys_dn = filter(cys_vs_lcl, log2FC <= (-log2fc_threshold)) %>%  pull(gene_name))

tat_cys_groups <- list(
  # UP + DOWN genes
  tat_cys = intersect(tat_cys$tat, tat_cys$cys),
  tat_only = setdiff(tat_cys$tat, tat_cys$cys),
  cys_only = setdiff(tat_cys$cys, tat_cys$tat),
  # UP genes
  UP_tat_cys = intersect(tat_cys$tat_up, tat_cys$cys_up),
  UP_tat_only = setdiff(tat_cys$tat_up, tat_cys$cys_up),
  UP_cys_only = setdiff(tat_cys$cys_up, tat_cys$tat_up),
  # DOWN genes
  DN_tat_cys = intersect(tat_cys$tat_dn, tat_cys$cys_dn),
  DN_tat_only = setdiff(tat_cys$tat_dn, tat_cys$cys_dn),
  DN_cys_only = setdiff(tat_cys$cys_dn, tat_cys$tat_dn))


tibble(
  UP_tat_cys = c(tat_cys_groups$UP_tat_cys, rep(NA_character_, 350 - length(tat_cys_groups$UP_tat_cys))),
  UP_tat_only = c(tat_cys_groups$UP_tat_only, rep(NA_character_, 350 - length(tat_cys_groups$UP_tat_only))),
  UP_cys_only = c(tat_cys_groups$UP_cys_only, rep(NA_character_, 350 - length(tat_cys_groups$UP_cys_only))),
  DN_tat_cys = c(tat_cys_groups$DN_tat_cys, rep(NA_character_, 350 - length(tat_cys_groups$DN_tat_cys))),
  DN_tat_only = c(tat_cys_groups$DN_tat_only, rep(NA_character_, 350 - length(tat_cys_groups$DN_tat_only))),
  DN_cys_only = c(tat_cys_groups$DN_cys_only, rep(NA_character_, 350 - length(tat_cys_groups$DN_cys_only)))) %>% 
  write_tsv(paste(output_dir, "genes.tsv", sep = "/"))


## Prepare list of UP and DOWN-regulated DEGs ----------------------------------
### Gene IDs conversion
geneIDs_table <- bitr(str_sub(all_counts$geneID, 1, 15), 
                      fromType="ENSEMBL", toType="ENTREZID", 
                      OrgDb="org.Hs.eg.db") %>% 
  dplyr::filter(!duplicated(ENTREZID), !duplicated(ENSEMBL)) %>% 
  dplyr::rename(ensID = ENSEMBL, entrezID = ENTREZID)

tat_vs_lcl_entrez <- inner_join(tat_vs_lcl, geneIDs_table, by = "ensID")
cys_vs_lcl_entrez <- inner_join(cys_vs_lcl, geneIDs_table, by = "ensID")

### All genes
tat_cys <- list(
  tat = filter(tat_vs_lcl_entrez, abs(log2FC) >= log2fc_threshold) %>%  pull(entrezID), 
  cys = filter(cys_vs_lcl_entrez, abs(log2FC) >= log2fc_threshold) %>%  pull(entrezID),
  tat_up = filter(tat_vs_lcl_entrez, log2FC >= log2fc_threshold) %>%  pull(entrezID), 
  cys_up = filter(cys_vs_lcl_entrez, log2FC >= log2fc_threshold) %>%  pull(entrezID), 
  tat_dn = filter(tat_vs_lcl_entrez, log2FC <= (-log2fc_threshold)) %>%  pull(entrezID),
  cys_dn = filter(cys_vs_lcl_entrez, log2FC <= (-log2fc_threshold)) %>%  pull(entrezID))

tat_cys_groups <- list(
  # UP + DOWN genes
  tat_cys = intersect(tat_cys$tat, tat_cys$cys),
  tat_only = setdiff(tat_cys$tat, tat_cys$cys),
  cys_only = setdiff(tat_cys$cys, tat_cys$tat),
  # UP genes
  UP_tat_cys = intersect(tat_cys$tat_up, tat_cys$cys_up),
  UP_tat_only = setdiff(tat_cys$tat_up, tat_cys$cys_up),
  UP_cys_only = setdiff(tat_cys$cys_up, tat_cys$tat_up),
  # DOWN genes
  DN_tat_cys = intersect(tat_cys$tat_dn, tat_cys$cys_dn),
  DN_tat_only = setdiff(tat_cys$tat_dn, tat_cys$cys_dn),
  DN_cys_only = setdiff(tat_cys$cys_dn, tat_cys$tat_dn))


## Run ORA against GO-BP database ----------------------------------------------
### ORA of UP/DOWN-regulated DEGs ----------------------------------------------
GO_tat_cys <- enrichGO(
  gene = tat_cys_groups$tat_cys,
  universe = geneIDs_table$entrezID,
  OrgDb = "org.Hs.eg.db",
  ont = "BP",
  pvalueCutoff = 0.05, 
  readable = TRUE)@result %>% 
  dplyr::filter(p.adjust < 0.05)

GO_tat_only <- enrichGO(
  gene = tat_cys_groups$tat_only,
  universe = geneIDs_table$entrezID,
  OrgDb = "org.Hs.eg.db",
  ont = "BP",
  pvalueCutoff = 0.05, 
  readable = TRUE)@result %>% 
  dplyr::filter(p.adjust < 0.05)

GO_cys_only <- enrichGO(
  gene = tat_cys_groups$cys_only,
  universe = geneIDs_table$entrezID,
  OrgDb = "org.Hs.eg.db",
  ont = "BP",
  pvalueCutoff = 0.05, 
  readable = TRUE)@result %>% 
  dplyr::filter(p.adjust < 0.05)

### ORA of UP-regulated DEGs ---------------------------------------------------
GO_UP_tat_cys <- enrichGO(
  gene = tat_cys_groups$UP_tat_cys,
  universe = geneIDs_table$entrezID,
  OrgDb = "org.Hs.eg.db",
  ont = "BP",
  pvalueCutoff = 0.05, 
  readable = TRUE)@result %>% 
  dplyr::filter(p.adjust < 0.05)

GO_UP_tat_only <- enrichGO(
  gene = tat_cys_groups$UP_tat_only,
  universe = geneIDs_table$entrezID,
  OrgDb = "org.Hs.eg.db",
  ont = "BP",
  pvalueCutoff = 0.05, 
  readable = TRUE)@result %>% 
  dplyr::filter(p.adjust < 0.05)

GO_UP_cys_only <- enrichGO(
  gene = tat_cys_groups$UP_cys_only,
  universe = geneIDs_table$entrezID,
  OrgDb = "org.Hs.eg.db",
  ont = "BP",
  pvalueCutoff = 0.05, 
  readable = TRUE)@result %>% 
  dplyr::filter(p.adjust < 0.05)

### ORA of DOWN-regulated DEGs -------------------------------------------------
GO_DN_tat_cys <- enrichGO(
  gene = tat_cys_groups$DN_tat_cys,
  universe = geneIDs_table$entrezID,
  OrgDb = "org.Hs.eg.db",
  ont = "BP",
  pvalueCutoff = 0.05, 
  readable = TRUE)@result %>% 
  dplyr::filter(p.adjust < 0.05)

GO_DN_tat_only <- enrichGO(
  gene = tat_cys_groups$DN_tat_only,
  universe = geneIDs_table$entrezID,
  OrgDb = "org.Hs.eg.db",
  ont = "BP",
  pvalueCutoff = 0.05, 
  readable = TRUE)@result %>% 
  dplyr::filter(p.adjust < 0.05)

GO_DN_cys_only <- enrichGO(
  gene = tat_cys_groups$DN_cys_only,
  universe = geneIDs_table$entrezID,
  OrgDb = "org.Hs.eg.db",
  ont = "BP",
  pvalueCutoff = 0.05, 
  readable = TRUE)@result %>% 
  dplyr::filter(p.adjust < 0.05)


## Save results ----------------------------------------------------------------
write_tsv(GO_tat_cys, paste(output_dir, "GO_tat_cys.tsv", sep = "/"))
write_tsv(GO_tat_only, paste(output_dir, "GO_tat_only.tsv", sep = "/"))
write_tsv(GO_cys_only, paste(output_dir, "GO_cys_only.tsv", sep = "/"))
write_tsv(GO_UP_tat_cys, paste(output_dir, "GO_UP_tat_cys.tsv", sep = "/"))
write_tsv(GO_UP_tat_only, paste(output_dir, "GO_UP_tat_only.tsv", sep = "/"))
write_tsv(GO_UP_cys_only, paste(output_dir, "GO_UP_cys_only.tsv", sep = "/"))
write_tsv(GO_DN_tat_cys, paste(output_dir, "GO_DN_tat_cys.tsv", sep = "/"))
write_tsv(GO_DN_tat_only, paste(output_dir, "GO_DN_tat_only.tsv", sep = "/"))
write_tsv(GO_DN_cys_only, paste(output_dir, "GO_DN_cys_only.tsv", sep = "/"))


## Plots -----------------------------------------------------------------------
### Dotplots -------------------------------------------------------------------
cat2show <- 10 # number of categories to show in the plot

GO_UP_tat_cys <- read_tsv(paste(output_dir, "GO_UP_tat_cys.tsv", sep = "/"))
GO_UP_cys_only <- read_tsv(paste(output_dir, "GO_UP_cys_only.tsv", sep = "/"))

GO_DN_tat_cys <- read_tsv(paste(output_dir, "GO_DN_tat_cys.tsv", sep = "/"))
GO_DN_tat_only <- read_tsv(paste(output_dir, "GO_DN_tat_only.tsv", sep = "/"))

dotpl <- function(df, x_limits, x_limits_labs){
  df %>% 
  dplyr::select(Description, p.adjust, Count, GeneRatio) %>% 
    mutate(
      GeneRatio_num = as.numeric(str_extract(GeneRatio, "^\\d+")),
      GeneRatio_den = as.numeric(str_extract(GeneRatio, "\\d+$")),
      GeneRatio = GeneRatio_num/GeneRatio_den,
      Description = as_factor(Description),
      Description = fct_reorder(Description, GeneRatio)) %>% 
    ggplot(aes(x = GeneRatio, y = Description, color = p.adjust, size = Count)) +
    geom_point() +
    scale_x_continuous(
      limits = x_limits,
      breaks = x_limits_labs) +
    scale_color_gradient(
      low = "red", high = "blue", 
      limits = c(0.00, 0.05),
      breaks = c(0.01, 0.02, 0.03, 0.04)) +
    scale_size(limits = c(4, 20), breaks = c(5, 10, 15)) +
    guides(
      color = guide_colorbar(order = 2, reverse = TRUE),
      size = guide_legend(order = 1)) +
    labs(y = "") +
    theme_bw() +
    theme(
      axis.text = element_text(color = "black")) +
    ggh4x::force_panelsizes(cols = unit(3, "cm"), rows = unit(4, "cm"))
}

dotpl_GO_UP_tat_cys <- GO_UP_tat_cys %>% 
  dotpl(x_limits = c(0.035, 0.115), x_limits_labs = c(0.04, 0.07, 0.1))

dotpl_GO_UP_cys_only <- GO_UP_cys_only %>%  
  slice_head(n = cat2show) %>% 
  dotpl(x_limits = c(0.035, 0.115), x_limits_labs = c(0.04, 0.07, 0.1))

dotpl_GO_DN_tat_cys <- GO_DN_tat_cys %>%  
  slice_head(n = cat2show) %>% 
  dotpl(x_limits = c(0.015, 0.08), x_limits_labs = c(0.02, 0.05, 0.08))

dotpl_GO_DN_tat_only <- GO_DN_tat_only %>%  
  slice_head(n = cat2show) %>% 
  dotpl(x_limits = c(0.02, 0.085), x_limits_labs = c(0.02, 0.05, 0.08))

### Save plots
#### UP Tat&Cys
ggsave(str_c(fig_dir, "GOBP_ORA_UP_dot_TatCys.png", sep = "/"), dotpl_GO_UP_tat_cys, units = "cm", width = 8, height = 5, scale = 3)
ggsave(str_c(fig_dir, "GOBP_ORA_UP_dot_TatCys.pdf", sep = "/"), dotpl_GO_UP_tat_cys, units = "cm", width = 8, height = 4, scale = 3)

#### UP Cys only
ggsave(str_c(fig_dir, "GOBP_ORA_UP_dot_Cysonly.png", sep = "/"), dotpl_GO_UP_cys_only, units = "cm", width = 8, height = 5, scale = 3)
ggsave(str_c(fig_dir, "GOBP_ORA_UP_dot_Cysonly.pdf", sep = "/"), dotpl_GO_UP_cys_only, units = "cm", width = 8, height = 4, scale = 3)

#### DN Tat&Cys
ggsave(str_c(fig_dir, "GOBP_ORA_DN_dot_TatCys.png", sep = "/"), dotpl_GO_DN_tat_cys, units = "cm", width = 8, height = 5, scale = 3)
ggsave(str_c(fig_dir, "GOBP_ORA_DN_dot_TatCys.pdf", sep = "/"), dotpl_GO_DN_tat_cys, units = "cm", width = 8, height = 4, scale = 3)

#### DN Cys only
ggsave(str_c(fig_dir, "GOBP_ORA_DN_dot_Tatonly.png", sep = "/"), dotpl_GO_DN_tat_only, units = "cm", width = 8, height = 5, scale = 3)
ggsave(str_c(fig_dir, "GOBP_ORA_DN_dot_Tatonly.pdf", sep = "/"), dotpl_GO_DN_tat_only, units = "cm", width = 8, height = 4, scale = 3)


### Heatmap --------------------------------------------------------------------
go_terms <- c("GO:0051607", "GO:0009615", "GO:0060337", "GO:0071357", "GO:0034340", "GO:0009615")

GO_UP_tat_cys <- read_tsv(paste(output_dir, "GO_UP_tat_cys.tsv", sep = "/"))
GO_UP_cys_only <- read_tsv(paste(output_dir, "GO_UP_cys_only.tsv", sep = "/"))

GO_UP_tat_cys_hm <- GO_UP_tat_cys %>% 
  dplyr::filter(ID %in% go_terms) %>% 
  dplyr::select(ID, Description, geneID) %>% 
  separate_rows(geneID) %>% 
  mutate(group = "tat_cys")

GO_UP_cys_only_hm <- GO_UP_cys_only %>% 
  dplyr::filter(ID %in% go_terms) %>% 
  dplyr::select(ID, Description, geneID) %>% 
  separate_rows(geneID) %>% 
  mutate(group = "cys_only")

GO_UP <- bind_rows(GO_UP_tat_cys_hm, GO_UP_cys_only_hm)

GO_UP %>% 
  ggplot(aes(x = geneID, y = Description, fill = group)) +
  geom_tile(color = "black") +
  scale_fill_manual(
    values = c("tat_cys" = "#bfafba", "cys_only" ="#F39B7F"),
    labels = c("Tat & Cys", "Cys")) +
  labs(x = "", y = "", fill = "") +
  theme_classic() +
  theme(
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 12),
    legend.position = "top"
  )

ggsave(paste(fig_dir, "GO_UP_Tat_Cys.png", sep = "/"), width = 10, height = 3, dpi = 300)


tat_vs_lcl_all <- read_tsv(str_c(input_dir, "Tat_vs_LCL.LFC.tsv", sep = "/"))
cys_vs_lcl_all <- read_tsv(str_c(input_dir, "Cys_vs_LCL.LFC.tsv", sep = "/"))

tat_vs_lcl_all %>% filter(gene_name == "CXCL10")
