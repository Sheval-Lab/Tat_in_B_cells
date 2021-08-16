## Load libraries --------------------------------------------------------------
library(tidyverse)
library(DESeq2)
library(gplots)

input_dir <- "output/tables/01_DGEA"
output_dir <- "output/tables/01_DGEA"
fig_dir <- "output/figures/01_DGEA"


## Load DESeq object -----------------------------------------------------------
dds <- readRDS(str_c(input_dir, "dds.rds", sep = "/"))


### Specify sample names
sample_names_format <- c(
  bquote("RPMI"^{"Tat"}), 
  bquote("RPMI"^{"Cys"}),  
  bquote("RPMI"^{"EGFP"}),
  bquote("RPMI"))

sample_names_format_all <- c(
  expression("RPMI"^{"Cys"}~1), expression("RPMI"^{"Cys"}~2), expression("RPMI"^{"Cys"}~3),  
  expression("RPMI"^{"EGFP"}~1), expression("RPMI"^{"EGFP"}~2), expression("RPMI"^{"EGFP"}~3),
  expression("RPMI"~1), expression("RPMI"~2), expression("RPMI"~3),
  expression("RPMI"^{"Tat"}~1), expression("RPMI"^{"Tat"}~2), expression("RPMI"^{"Tat"}~3))

sample_names <- c(
  "Tat.1", "Tat.2", "Tat.3",
  "Cys.1", "Cys.2", "Cys.3",
  "GFP.1", 'GFP.2', "GFP.3",
  "LCL.1", "LCL.2", "LCL.3")

## Visualize expression data with PCA ------------------------------------------
cols4 <- c("#CC79A7", "#D55E00", "#009E73", "#0072B2")

pca <- plotPCA(rlog(dds), intgroup = "sample", ntop = 500, returnData = TRUE)
pca_percent_var <- round(100 * attr(pca, "percentVar"))


ggplot(pca, aes(PC1, PC2, color = sample)) + 
  geom_point(size = 4) +
  xlab(paste0("PC1: ", pca_percent_var[1], "% variance")) + 
  ylab(paste0("PC2: ", pca_percent_var[2], "% variance")) +
  scale_color_manual(
    name = "",
    values = c("Tat" = cols4[1], "Cys" = cols4[2], "GFP" = cols4[3], "LCL" = cols4[4]),
    labels = sample_names_format) +
  theme(text = element_text(size = 12, color = "black"), 
        panel.background = element_rect(fill = "white", colour = "black", size = 0.5, linetype = "solid"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.key = element_rect(fill = "white"),
        aspect.ratio = 1)

ggsave(str_c(fig_dir, "PCA.png", sep = "/"), units = "cm", width = 5, height = 5, scaling = 0.5)
ggsave(str_c(fig_dir, "PCA.pdf", sep = "/"), units = "cm", width = 5, height = 5, scale = 2)


## Visualize cross-sample correlation with heatmap -----------------------------
### Calculate pairwise correlation
sample_cor <- round(cor(counts(dds, normalized = TRUE),method = "spearman" ), 3)


### Plot with heatmap.2
lhei <-  c(1, 4)

#### PNG
ragg::agg_png(str_c(fig_dir, "heatmap.png", sep = "/"), units = "cm", width = 15, height = 11, res = 300, scaling = 3/4)
heatmap.2(
  sample_cor, 
  col = rev(heat.colors(1000)), 
  main = "Spearman Correlation", 
  cellnote = sample_cor, 
  notecol = "black", 
  notecex = 0.8, 
  trace = "none",
  density.info = "none",
  dendrogram = "row", 
  key.xlab = "",
  cexCol = 1, 
  cexRow = 0.9, 
  margins = c(7,7), 
  lhei = lhei,
  labCol = sample_names_format_all, 
  labRow = sample_names_format_all)
dev.off()

#### PDF
pdf(str_c(fig_dir, "heatmap.pdf", sep = "/"), width = 8, height = 6)
heatmap.2(
  sample_cor, 
  col = rev(heat.colors(1000)), 
  main = "Spearman Correlation", 
  cellnote = sample_cor, 
  notecol = "black", 
  notecex = 0.8, 
  trace = "none",
  density.info = "none",
  dendrogram = "row", 
  key.xlab = "",
  cexCol = 1, 
  cexRow = 0.9, 
  margins = c(7,7), 
  lhei = lhei,
  labCol = sample_names_format_all, 
  labRow = sample_names_format_all)
dev.off()


## Visualization of normalization effect on count distribution -----------------
### Import raw and normalized counts
counts_flt <- read_tsv(str_c(input_dir, "processed_counts", "counts_flt.tsv", sep = "/"))
counts_norm <- read_tsv(str_c(input_dir, "processed_counts", "counts_norm.tsv", sep = "/"))

counts_flt_long <- counts_flt %>% 
  pivot_longer(cols = -c("geneID", "ensID"), names_to = "sample", values_to = "count") %>% 
  select(-c("geneID", "ensID")) %>% 
  mutate(
    count = log2(count + 1),
    state = "Before normalization")

counts_norm_long <- counts_norm %>% 
  pivot_longer(cols = -c("geneID"), names_to = "sample", values_to = "count") %>% 
  select(-c("geneID")) %>% 
  mutate(
    count = log2(count + 1),
    state = "After normalization")

counts_long <- rbind(counts_flt_long, counts_norm_long) %>% 
  mutate(
    sample = factor(sample, levels = sample_names),
    state = factor(state, levels = c("Before normalization", "After normalization")))


ggplot(counts_long, aes(x = sample, y = count)) + 
  geom_boxplot() + 
  facet_wrap(~state) +
  labs(x = "", y = expression("log"[10]*"(count+1)")) +
  scale_x_discrete(labels = rep(sample_names_format_all[c(10:12, 1:9)], times = 2)) +
  ## Edit theme
  theme(text = element_text(size = 12, color = "black"), 
        axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
        axis.text.y = element_text(color = "black"),
        panel.background = element_rect(fill = "white", colour = "black", size = 0.5, linetype = "solid"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill="white", colour = "black", size = 0.5, linetype = "solid"),
        strip.text = element_text(face="bold", size = 10), 
        legend.position = "top",
        aspect.ratio = 1)
  
ggsave(str_c(fig_dir, "DESeq_normalization.png", sep = "/"), units = "cm", width = 10, scaling = 2/3)
ggsave(str_c(fig_dir, "DESeq_normalization.pdf", sep = "/"), units = "cm", width = 10, scale = 3/2)

