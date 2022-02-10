## Load libraries --------------------------------------------------------------
library(tidyverse)
library(magrittr)
library(DESeq2)
source("scripts/functions/extract_results.R")

input_dir <- "data/counts"
output_dir <- "output/tables/04_DGEA-EBV"
fig_dir <- "output/figures/04_DGEA-EBV"

log2FC_threshold = log2(1.5)

## Load gene annotation --------------------------------------------------------
gene_annotation_file <- "data/metadata/GRCh38.p10_ALL.annotation.IDs.txt"
gene_annotation <- read_tsv(gene_annotation_file,
                            col_names = c("geneID", "gene_name", "gene_type"))

## Load counts data ------------------------------------------------------------
ebv_counts_files <- list.files(input_dir, full.names = TRUE, pattern = "r.ebv_e.tab$")

### Extract sample names from file names
sample_names <- ebv_counts_files %>% 
  str_split("/") %>% 
  map_chr(last) %>% 
  str_remove(".r.ebv_e.tab")

### Read multiple htseq count files into one dataframe
ebv_counts <- map_dfc(ebv_counts_files, read_tsv, col_names = FALSE) %>% 
  select(1, where(is.numeric)) %>% 
  set_colnames(c("geneID", sample_names)) %>% 
  ## Filter out summary rows
  filter(!str_detect(geneID, "^__"))

### Pre-filter non-expressed genes
ebv_counts_flt <- ebv_counts %>% 
  filter(if_any(where(is.numeric), ~ . > 0))

#### Save pre-filtered counts data
write_tsv(ebv_counts_flt, str_c(output_dir, "processed_counts", "EBV_counts_flt.tsv", sep = "/"))

## Prepare data for DGEA -------------------------------------------------------
### Prepare matrix of counts
tat_counts_flt <- read_tsv("output/tables/01_DGEA/processed_counts/counts_flt.tsv")

ebv_counts_flt <- read_tsv(str_c(output_dir, "processed_counts", "EBV_counts_flt.tsv", sep = "/"))
ebv_counts_mtx <- ebv_counts_flt %>% 
  bind_rows(tat_counts_flt[,-2]) %>% 
  select_if(is.numeric) %>% 
  as.matrix()
rownames(ebv_counts_mtx) <- c(ebv_counts_flt$geneID, tat_counts_flt$geneID)

### Prepare design dataframe
design <- data.frame(
  sample = str_remove(colnames(ebv_counts_mtx), ".\\d"),
  row.names = colnames(ebv_counts_mtx))

### Create DESeq object
ebv_dds <- DESeqDataSetFromMatrix(
  countData = ebv_counts_mtx, 
  colData = design, 
  design = ~sample)

### Use size factors
### Load precomputed size factors
# sf <- dget("output/tables/01_DGEA/processed_counts/size_factors.tsv")
# sizeFactors(ebv_dds) <- sf


## Run DGEA --------------------------------------------------------------------
ebv_dds <- DESeq(ebv_dds)


## Extract DGEA results --------------------------------------------------------
### Create contrasts list
contrasts_list <- tibble(
  numerator = design$sample,
  denominator = design$sample) %>% 
  tidyr::expand(numerator, denominator)

contrasts_list <- contrasts_list[c(3,7,15,13),]

### Extract results
#### Tat vs Cys
walk2(contrasts_list$numerator[4], contrasts_list$denominator[4], 
      extract_results, dds = ebv_dds, gene_annotation = gene_annotation, log2FC_threshold = log2FC_threshold)

#### Tat vs LCL, Cys vs LCL, GFP vs LCL
ebv_dds$sample <- relevel(ebv_dds$sample, ref = "LCL")
ebv_dds <- DESeq(ebv_dds)

walk2(contrasts_list$numerator[1:3], contrasts_list$denominator[1:3], 
      extract_results, dds = ebv_dds, gene_annotation = gene_annotation, log2FC_threshold = log2FC_threshold)


## Visualize results -----------------------------------------------------------
### 'Tat vs LCL' ---------------------------------------------------------------
### Read in Deseq2 results 'Tat vs LCL'
tat_ebv <- read_tsv(str_c(output_dir, "deseq", "Tat_vs_LCL.LFC.tsv", sep = "/"))
tat_ebv_out <- tat_ebv %>% filter(!str_detect(geneID, "ENSG"))

write_tsv(tat_ebv_out, str_c(output_dir, "deseq", "Tat_vs_LCL.LFC.EBV_only.tsv", sep = "/"))

### Volcano plot 'Tat vs LCL'
tat_ebv_vln <- tat_ebv %>% 
  select(geneID, padj, log2FC) %>% 
  mutate(organism = if_else(str_detect(geneID, "ENSG"), "HS", "EBV"))

x_axis <- c(-9, -6, -3, 0, 3, 6, 9)

### Highlight EBV genes
ggplot() +
  geom_point(
    filter(tat_ebv_vln, organism == "HS"), 
    mapping = aes(x = log2FC, y = -log10(padj)), 
    color = "grey", alpha = 0.6) +
  geom_point(
    filter(tat_ebv_vln, organism == "EBV"), 
    mapping = aes(x = log2FC, y = -log10(padj)), 
    color = "red", alpha = 0.9) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_vline(xintercept = log2FC_threshold, linetype = "dashed", color = "black") +
  geom_vline(xintercept = -log2FC_threshold, linetype = "dashed", color = "black") +
  ## Edit axis names
  ylab(expression(-log[10]~padj)) +
  xlab(expression(log[2]~FoldChange)) +
  ## Edit color scheme
  # scale_color_manual(name = "", values = cols, labels = group_names) +
  ## Edit X axis
  scale_x_continuous(breaks = x_axis, limits = c(-max(abs(tat_ebv_vln$log2FC)), max(abs(tat_ebv_vln$log2FC)))) +
  theme_bw() +
  theme(
    aspect.ratio = 1, 
    text = element_text(size = 12, color = "black"), 
    # axis.text.x = element_text(color = "black"),
    # axis.text.y = element_text(color = "black"),
    panel.background = element_rect(fill = "white", colour = "black", size = 0.5, linetype = "solid"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank())


ggsave(str_c(fig_dir, "EBV_volcano.png", sep = "/"), units = "cm", width = 16)
ggsave(str_c(fig_dir, "EBV_volcano.pdf", sep = "/"), units = "cm", width = 16)

### 'Cys vs LCL' ---------------------------------------------------------------
### Read in Deseq2 results 'Cys vs LCL'
cys_ebv <- read_tsv(str_c(output_dir, "deseq", "Cys_vs_LCL.LFC.tsv", sep = "/"))
cys_ebv_out <- cys_ebv %>% filter(!str_detect(geneID, "ENSG"))

write_tsv(cys_ebv_out, str_c(output_dir, "deseq", "Cys_vs_LCL.LFC.EBV_only.tsv", sep = "/"))

### Volcano plot 'Tat vs LCL'
cys_ebv_vln <- cys_ebv %>% 
  select(geneID, padj, log2FC) %>% 
  mutate(organism = if_else(str_detect(geneID, "ENSG"), "HS", "EBV"))

x_axis <- c(-9, -6, -3, 0, 3, 6, 9)

### Highlight EBV genes
ggplot() +
  geom_point(
    filter(cys_ebv_vln, organism == "HS"), 
    mapping = aes(x = log2FC, y = -log10(padj)), 
    color = "grey", alpha = 0.6) +
  geom_point(
    filter(cys_ebv_vln, organism == "EBV"), 
    mapping = aes(x = log2FC, y = -log10(padj)), 
    color = "red", alpha = 0.9) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_vline(xintercept = log2FC_threshold, linetype = "dashed", color = "black") +
  geom_vline(xintercept = -log2FC_threshold, linetype = "dashed", color = "black") +
  ## Edit axis names
  ylab(expression(-log[10]~padj)) +
  xlab(expression(log[2]~FoldChange)) +
  ## Edit color scheme
  # scale_color_manual(name = "", values = cols, labels = group_names) +
  ## Edit X axis
  scale_x_continuous(breaks = x_axis, limits = c(-max(abs(cys_ebv_vln$log2FC)), max(abs(cys_ebv_vln$log2FC)))) +
  theme_bw() +
  theme(
    aspect.ratio = 1, 
    text = element_text(size = 12, color = "black"), 
    # axis.text.x = element_text(color = "black"),
    # axis.text.y = element_text(color = "black"),
    panel.background = element_rect(fill = "white", colour = "black", size = 0.5, linetype = "solid"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank())


ggsave(str_c(fig_dir, "EBV_Cys_volcano.png", sep = "/"), units = "cm", width = 16)
ggsave(str_c(fig_dir, "EBV_Cys_volcano.pdf", sep = "/"), units = "cm", width = 16)
