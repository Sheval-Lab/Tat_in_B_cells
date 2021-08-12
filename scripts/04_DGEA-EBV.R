## Load libraries --------------------------------------------------------------
library(tidyverse)
library(magrittr)
library(DESeq2)
source("scripts/functions/extract_results.R")

input_dir <- "data/counts"
output_dir <- "output/tables/04_DGEA-EBV"
fig_dir <- "output/figures/04_DGEA-EBV"

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
### Load precomputed size factors
sf <- dget("output/tables/01_DGEA/processed_counts/size_factors.tsv")

### Prepare matrix of counts
ebv_counts_flt <- read_tsv(str_c(output_dir, "processed_counts", "EBV_counts_flt.tsv", sep = "/"))
ebv_counts_mtx <- ebv_counts_flt %>% 
  select_if(is.numeric) %>% 
  as.matrix()
rownames(ebv_counts_mtx) <- ebv_counts_flt$geneID

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
sizeFactors(ebv_dds) <- sf


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
log2FC_threshold = log2(1.5)
#### Tat vs Cys
walk2(contrasts_list$numerator[4], contrasts_list$denominator[4], 
      extract_results, dds = ebv_dds, gene_annotation = gene_annotation, log2FC_threshold = log2FC_threshold)

#### Tat vs LCL, Cys vs LCL, GFP vs LCL
ebv_dds$sample <- relevel(ebv_dds$sample, ref = "LCL")
ebv_dds <- DESeq(ebv_dds)

walk2(contrasts_list$numerator[1:3], contrasts_list$denominator[1:3], 
      extract_results, dds = ebv_dds, gene_annotation = gene_annotation, log2FC_threshold = log2FC_threshold)


