## Load libraries --------------------------------------------------------------
library(tidyverse)
library(clusterProfiler)


input_dir <- "output/tables/01_DGEA/deseq"
output_dir <- "output/tables/03_GE/TF"
fig_dir <- "output/figures/03_GE/TF"


## Specify log2FoldChange threshold for DEGs -----------------------------------
log2fc_threshold <- log2(1.5)


## Load DGEA results -----------------------------------------------------------
tat_vs_lcl <- read_tsv(str_c(input_dir, "Tat_vs_LCL.tsv", sep = "/"))
tat_vs_lcl_LFC <- read_tsv(str_c(input_dir, "Tat_vs_LCL.LFC.tsv", sep = "/"))


## Prepare input data for TFs search -------------------------------------------
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


### Prepare gene lists of UP and DOWN regulated DEGs ---------------------------
### Fold change >= 1.5
genes_up <- tat %>% filter(padj < 0.05, log2FC >= log2fc_threshold) %>% pull(entrezID)
genes_dn <- tat %>% filter(padj < 0.05, log2FC <= (-log2fc_threshold)) %>% pull(entrezID)


## Load GO (MF) annotation
go_annot <- select(org.Hs.eg.db, keys = keys(org.Hs.eg.db), columns = "GO")


### Extract genes with annotated GO Molecular function GO:0003700
### 'DNA-binding transcription factor activity'
go_annot_tf <- go_annot %>% filter(GO == "GO:0003700")

### Extract upregulated TFs
go_annot_tf_up <- go_annot_tf %>% filter(ENTREZID %in% genes_up)

go_annot_tf_up <- bitr(go_annot_tf_up$ENTREZID, fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Hs.eg.db") %>% 
  left_join(go_annot_tf_up, by = "ENTREZID") %>% 
  filter(!duplicated(ENTREZID)) %>% 
  arrange(SYMBOL)

go_annot_tf_up

### Extract downregulated TFs
go_annot_tf_dn <- go_annot_tf %>% filter(ENTREZID %in% genes_dn)

go_annot_tf_dn <- bitr(go_annot_tf_dn$ENTREZID, fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Hs.eg.db") %>% 
  left_join(go_annot_tf_dn, by = "ENTREZID") %>% 
  filter(!duplicated(ENTREZID)) %>% 
  arrange(SYMBOL)

go_annot_tf_dn
