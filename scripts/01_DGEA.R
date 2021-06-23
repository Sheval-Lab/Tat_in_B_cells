## Load libraries --------------------------------------------------------------
library(tidyverse)
library(magrittr)
library(DESeq2)
source("scripts/functions/extract_results.R")

input_dir <- "data/counts"
output_dir <- "output/tables/01_DGEA"
fig_dir <- "output/figures/01_DGEA"

## Load gene annotation --------------------------------------------------------
gene_annotation_file <- "data/metadata/GRCh38.p10_ALL.annotation.IDs.txt"
gene_annotation <- read_tsv(gene_annotation_file, 
                            col_names = c("geneID", "gene_name", "gene_type"))

## Load counts data ------------------------------------------------------------
counts_files <- list.files(input_dir, full.names = TRUE)

### Extract sample names from file names
sample_names <- counts_files %>% 
  str_split("/") %>% 
  map_chr(last) %>% 
  str_remove(".r.tab")

### Read multiple htseq count files into one dataframe
counts <- map_dfc(counts_files, read_tsv, col_names = FALSE) %>% 
  select(1, where(is.numeric)) %>% 
  set_colnames(c("geneID", sample_names)) %>% 
  ## Filter out summary rows
  filter(!str_detect(geneID, "^__"))

### Pre-filter non-expressed genes
counts_flt <- counts %>% 
  filter(if_any(where(is.numeric), ~ . > 0))

### Pre-filter highly expressed ribosomal protein genes
#### Save list of top 100 highly expressed genes
counts_flt %>% 
  mutate(
    ensID = str_remove(geneID, ".\\d+$"),
    rowsum = rowSums(select_if(., is.numeric), na.rm = TRUE)) %>% 
  arrange(desc(rowsum)) %>% 
  select(ensID) %>% 
  head(100) %>% 
  write_tsv(str_c(output_dir, "processed_counts", "top_100_genes.tsv", sep = "/"))

#### After annotating top 100 genes with DAVID, read in list with highly expressed ribosomal protein genes
ribo_genes <- read_tsv(str_c(output_dir, "processed_counts", "ribo_genes.tsv", sep = "/"))

#### Exclude ribosomal protein genes from counts data
counts_flt %<>%
  mutate(ensID = str_remove(geneID, ".\\d+$")) %>% 
  filter(!(ensID %in% ribo_genes$ID)) %>% 
  select(geneID, ensID, everything())

#### Save pre-filtered counts data
write_tsv(counts_flt, str_c(output_dir, "processed_counts", "counts_flt.tsv", sep = "/"))


## Prepare data for DGEA -------------------------------------------------------
### Prepare matrix of counts
counts_flt <- read_tsv(str_c(output_dir, "processed_counts", "counts_flt.tsv", sep = "/"))
counts_mtx <- counts_flt %>% 
  select_if(is.numeric) %>% 
  as.matrix()
rownames(counts_mtx) <- counts_flt$geneID

### Prepare design dataframe
design <- data.frame(
  sample = str_remove(colnames(counts_mtx), ".\\d"),
  row.names = colnames(counts_mtx))

### Create DESeq object
dds <- DESeqDataSetFromMatrix(
  countData = counts_mtx, 
  colData = design, 
  design = ~sample)


### Calculate size factors and save normalized counts
dds <- estimateSizeFactors(dds)
sizeFactors(dds)

counts_norm <- counts(dds, normalized=TRUE) %>% 
  as.data.frame() %>% 
  mutate(geneID = row.names(.)) %>% 
  select(geneID, everything())
  
write_tsv(counts_norm, str_c(output_dir, "processed_counts", "counts_norm.tsv", sep = "/"))


## Run DGEA --------------------------------------------------------------------
dds <- DESeq(dds)

saveRDS(dds, str_c(output_dir, "dds.rds", sep = "/"))


## Extract DGEA results --------------------------------------------------------
### Create contrasts list
contrasts_list <- tibble(
  numerator = design$sample,
  denominator = design$sample) %>% 
  tidyr::expand(numerator, denominator)

contrasts_list <- contrasts_list[c(3,7,15,13),]


### Extract results
#### Tat vs Cys
walk2(contrasts_list$numerator[4], contrasts_list$denominator[4], extract_results, dds = dds, gene_annotation = gene_annotation)

#### Tat vs LCL, Cys vs LCL, GFP vs LCL
dds$sample <- relevel(dds$sample, ref = "LCL")
dds <- DESeq(dds)

walk2(contrasts_list$numerator[1:3], contrasts_list$denominator[1:3], extract_results, dds = dds, gene_annotation = gene_annotation)

