## Load libraries --------------------------------------------------------------
library(tidyverse)

input_dir <- "output/tables/01_DGEA/deseq"
output_dir <- "output/tables/02_DGEA-plots"
fig_dir <- "output/figures/02_DGEA-plots"


## Specify log2FoldChange threshold for DEGs -----------------------------------
log2FC_threshold <- log2(1.5)


## Load DGEA results -----------------------------------------------------------
tat_vs_LCL_LFC <- read_tsv(str_c(input_dir, "Tat_vs_LCL.LFC.tsv", sep = "/"))
cys_vs_LCL_LFC <- read_tsv(str_c(input_dir, "Cys_vs_LCL.LFC.tsv", sep = "/"))
GFP_vs_LCL_LFC <- read_tsv(str_c(input_dir, "GFP_vs_LCL.LFC.tsv", sep = "/"))


gene_list <- c("CD19", "CD38", "FCER2", "HLA-DRA", "MME", "MS4A1", "PTPRC")
gene_list_conv <- c("CD19", "CD38", "CD23", "HLA-DR", "CD10", "CD20", "CD45")

tat_vs_LCL_genes <- tat_vs_LCL_LFC %>% 
  filter(gene_name %in% gene_list, ensID != "ENSG00000262418") %>% 
  arrange(gene_name) %>% 
  mutate(gene_name_conv = gene_list_conv)

cys_vs_LCL_genes <- cys_vs_LCL_LFC %>% 
  filter(gene_name %in% gene_list, ensID != "ENSG00000262418") %>% 
  arrange(gene_name) %>% 
  mutate(gene_name_conv = gene_list_conv)

GFP_vs_LCL_genes <- GFP_vs_LCL_LFC %>% 
  filter(gene_name %in% gene_list, ensID != "ENSG00000262418") %>% 
  arrange(gene_name) %>% 
  mutate(gene_name_conv = gene_list_conv)

## Make volcano plot -----------------------------------------------------------
x_axis <- c(-10, -5, 0, 5, 10)

### Tat vs LCL
ggplot(tat_vs_LCL_LFC, aes(x = log2FC, y = -log10(padj))) +
  geom_point(alpha = 0.5, color = "grey") +
  geom_point(
    tat_vs_LCL_genes, 
    mapping = aes(x = log2FC, y = -log10(padj)),
    color = "red") +
  ggrepel::geom_text_repel(
    tat_vs_LCL_genes[tat_vs_LCL_genes$padj < 0.05,],
    mapping = aes(label = gene_name_conv), 
    size = 4, 
    min.segment.length = 0, 
    ylim  = c(-log10(0.05), NA),
    box.padding = 0.55) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_vline(xintercept = log2FC_threshold, linetype = "dashed", color = "black") +
  geom_vline(xintercept = -log2FC_threshold, linetype = "dashed", color = "black") +
  ## Edit axis names
  ylab(expression(-log[10]~padj)) +
  xlab(expression(log[2]~FoldChange)) +
  ## Edit X axis
  scale_x_continuous(breaks = x_axis, limits = c(-max(abs(tat_vs_LCL_LFC$log2FC)), max(abs(tat_vs_LCL_LFC$log2FC)))) +
  theme_bw() +
  theme(
    aspect.ratio = 1, 
    text = element_text(size = 12, color = "black"), 
    panel.background = element_rect(fill = "white", colour = "black", size = 0.5, linetype = "solid"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank())


ggsave(str_c(fig_dir, "Tat_vs_LCL_markers_volcano.png", sep = "/"), units = "cm", width = 16)
ggsave(str_c(fig_dir, "Tat_vs_LCL_markers_volcano.pdf", sep = "/"), units = "cm", width = 16)


### Cys vs LCL
ggplot(cys_vs_LCL_LFC, aes(x = log2FC, y = -log10(padj))) +
  geom_point(alpha = 0.5, color = "grey") +
  geom_point(
    cys_vs_LCL_genes, 
    mapping = aes(x = log2FC, y = -log10(padj)),
    color = "red") +
  ggrepel::geom_text_repel(
    cys_vs_LCL_genes[cys_vs_LCL_genes$padj < 0.05,],
    mapping = aes(label = gene_name_conv), 
    size = 4, 
    min.segment.length = 0, 
    ylim  = c(-log10(0.05), NA),
    box.padding = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_vline(xintercept = log2FC_threshold, linetype = "dashed", color = "black") +
  geom_vline(xintercept = -log2FC_threshold, linetype = "dashed", color = "black") +
  ## Edit axis names
  ylab(expression(-log[10]~padj)) +
  xlab(expression(log[2]~FoldChange)) +
  ## Edit X axis
  scale_x_continuous(breaks = x_axis, limits = c(-max(abs(cys_vs_LCL_LFC$log2FC)), max(abs(cys_vs_LCL_LFC$log2FC)))) +
  theme_bw() +
  theme(
    aspect.ratio = 1, 
    text = element_text(size = 12, color = "black"), 
    panel.background = element_rect(fill = "white", colour = "black", size = 0.5, linetype = "solid"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank())


ggsave(str_c(fig_dir, "Cys_vs_LCL_markers_volcano.png", sep = "/"), units = "cm", width = 16)
ggsave(str_c(fig_dir, "Cys_vs_LCL_markers_volcano.pdf", sep = "/"), units = "cm", width = 16)

### Tat vs LCL
ggplot(GFP_vs_LCL_LFC, aes(x = log2FC, y = -log10(padj))) +
  geom_point(alpha = 0.5, color = "grey") +
  geom_point(
    GFP_vs_LCL_genes, 
    mapping = aes(x = log2FC, y = -log10(padj)),
    color = "red") +
  ggrepel::geom_text_repel(
    GFP_vs_LCL_genes[GFP_vs_LCL_genes$padj < 0.05,],
    mapping = aes(label = gene_name_conv), 
    size = 4, 
    min.segment.length = 0, 
    ylim  = c(-log10(0.05), NA),
    box.padding = 0.55) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_vline(xintercept = log2FC_threshold, linetype = "dashed", color = "black") +
  geom_vline(xintercept = -log2FC_threshold, linetype = "dashed", color = "black") +
  ## Edit axis names
  ylab(expression(-log[10]~padj)) +
  xlab(expression(log[2]~FoldChange)) +
  ## Edit X axis
  scale_x_continuous(breaks = x_axis, limits = c(-max(abs(GFP_vs_LCL_LFC$log2FC)), max(abs(GFP_vs_LCL_LFC$log2FC)))) +
  theme_bw() +
  theme(
    aspect.ratio = 1, 
    text = element_text(size = 12, color = "black"), 
    panel.background = element_rect(fill = "white", colour = "black", size = 0.5, linetype = "solid"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank())


ggsave(str_c(fig_dir, "GFP_vs_LCL_markers_volcano.png", sep = "/"), units = "cm", width = 16)
ggsave(str_c(fig_dir, "GFP_vs_LCL_markers_volcano.pdf", sep = "/"), units = "cm", width = 16)


## Make table ---------------------------------------------
gene_table <- bind_cols(
  select(tat_vs_LCL_genes, gene_name, gene_name_conv, Tat_vs_LCL.padj = padj, Tat_vs_LCL.FC = FC),
  select(cys_vs_LCL_genes, Cys_vs_LCL.padj = padj, Cys_vs_LCL.FC = FC),
  select(GFP_vs_LCL_genes, GFP_vs_LCL.padj = padj, GFP_vs_LCL.FC = FC))

write_tsv(gene_table, str_c(output_dir, "markers_expression.tsv", sep = "/"))


