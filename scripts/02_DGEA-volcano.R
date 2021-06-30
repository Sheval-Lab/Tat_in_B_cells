## Load libraries --------------------------------------------------------------
library(tidyverse)

input_dir <- "output/tables/01_DGEA/deseq"
output_dir <- "output/tables/02_DGEA-plots"
fig_dir <- "output/figures/02_DGEA-plots"


## Specify log2FoldChange threshold for DEGs -----------------------------------
log2fc_threshold <- log2(1.5)


## Load DGEA results -----------------------------------------------------------
tat_vs_cys <- read_tsv(str_c(input_dir, "Tat_vs_Cys.tsv", sep = "/"))


## Prepare dataframe for plotting volcano plot ---------------------------------
vln_df <- tat_vs_cys %>% 
  transmute(
    geneID, gene_name, padj, log2FC,
    gene_type = ifelse(gene_type == "protein_coding", "protein_coding", "other"),
    group = case_when(
      # ((padj < 0.05) & (log2FC >= log2fc_threshold) & (gene_type == "protein_coding")) ~ "up_fc2_pc",
      # ((padj < 0.05) & (log2FC <= (-log2fc_threshold)) & (gene_type == "protein_coding")) ~ "dn_fc2_pc",
      ((padj < 0.05) & (log2FC >= log2fc_threshold)) ~ "up_fc2",
      ((padj < 0.05) & (log2FC <= (-log2fc_threshold))) ~ "dn_fc2",
      ((padj < 0.05) & (log2FC >= 0)) ~ "up",
      ((padj < 0.05) & (log2FC <= 0)) ~ "dn",
      TRUE ~ "not_significant"),
    group = factor(group, levels = c("up_fc2", "up", "dn_fc2",  "dn", "not_significant")))

write_tsv(vln_df, str_c(output_dir, "Tat_vs_Cys.vln_df.tsv", sep = "/"))

## Make volcano plot -----------------------------------------------------------
cols = c(
  "up_fc2" = "#db3a34", 
  "up" = "#f7a399",
  "dn_fc2" = "#03045e", 
  "dn" = "#00b4d8",
  "not_significant" = "grey")

group_names <- c(
  "p.adj < 0.05\nlog2FC \u2265 1",
  "p.adj < 0.05\nlog2FC \u2265 0",
  "p.adj < 0.05\nlog2FC \u2264 -1", 
  "p.adj < 0.05\nlog2FC \u2264 0",
  "Not significant")

x_axis <- c(-3, -2, -1, 0, 1, 2, 3)

ggplot(vln_df, aes(x = log2FC, y = -log10(padj), color = group)) +
  geom_point(alpha = 0.9) +
  # geom_text_repel(filter(res_def_df, gene_expr %in% c("up_FC2", "down_FC2")),
  #                 mapping = aes(label = gene_name), size = 3) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_vline(xintercept = log2fc_threshold, linetype = "dashed", color = "black") +
  geom_vline(xintercept = -log2fc_threshold, linetype = "dashed", color = "black") +
  ## Edit axis names
  ylab(expression(-log[10]~padj)) +
  xlab(expression(log[2]~FoldChange)) +
  ## Edit color scheme
  scale_color_manual(name = "", values = cols, labels = group_names) +
  ## Edit X axis
  scale_x_continuous(breaks = x_axis, limits = c(-max(abs(vln_df$log2FC)), max(abs(vln_df$log2FC)))) +
  theme_bw() +
  theme(
    aspect.ratio = 1, 
    text = element_text(size = 12, color = "black"), 
    # axis.text.x = element_text(color = "black"),
    # axis.text.y = element_text(color = "black"),
    panel.background = element_rect(fill = "white", colour = "black", size = 0.5, linetype = "solid"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank())


ggsave(str_c(fig_dir, "DEGs_volcano.png", sep = "/"), units = "cm", width = 16)
ggsave(str_c(fig_dir, "DEGs_volcano.svg", sep = "/"), units = "cm", width = 16)
