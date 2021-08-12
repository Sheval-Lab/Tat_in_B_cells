## Load libraries --------------------------------------------------------------
library(tidyverse)
library(org.Hs.eg.db)
library(GSVA)
library(limma)

input_dir <- "output/tables/01_DGEA"
output_dir <- "output/tables/04_GSVA"
fig_dir <- "output/figures/04_GSVA"


## Load normalized counts ------------------------------------------------------
counts_norm <- read_tsv(str_c(input_dir, "processed_counts", "counts_norm.tsv", sep = "/"))
counts_norm$ensID <- str_remove(counts_norm$geneID, ".\\d+$")

### Gene IDs conversion
geneIDs_table <- clusterProfiler::bitr(counts_norm$ensID, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db") %>% 
  filter(!duplicated(ENTREZID), !duplicated(ENSEMBL)) %>% 
  dplyr::rename(ensID = ENSEMBL, entrezID = ENTREZID)

counts_norm_flt <- counts_norm %>% 
  right_join(geneIDs_table, by = "ensID")

counts_norm_mtx <- counts_norm_flt %>% 
  dplyr::select(where(is.numeric)) %>% 
  as.matrix()
rownames(counts_norm_mtx) <- counts_norm_flt$entrezID


## Prepare GO BP annotation as a list ------------------------------------------
go_annot <- select(org.Hs.eg.db, keys = keys(org.Hs.eg.db), columns = "GO")
head(go_annot)
gobp_annot <- go_annot %>% filter(ONTOLOGY == "BP")

genes_by_gobp <- split(gobp_annot$ENTREZID, gobp_annot$GO)
length(genes_by_gobp)


## Run GSVA on GO BP database --------------------------------------------------
gobp_gsva <- gsva(counts_norm_mtx, genes_by_gobp, min.sz = 10, max.sz = 500)

saveRDS(gobp_gsva, str_c(output_dir, "gobp_gsva.rds", sep = "/"))

gobp_gsva <- readRDS(str_c(output_dir, "gobp_gsva.rds", sep = "/"))


## Differential expression at GO BP level --------------------------------------
design <- data.frame(
  sample = str_remove(colnames(counts_norm_mtx), ".\\d"),
  row.names = colnames(counts_norm_mtx)) 

design_flt <- design %>% 
  filter(sample %in% c("Tat", "Cys"))

mod <- model.matrix(~sample, data = design_flt)


fit <- lmFit(gobp_gsva[,c(1:3, 10:12)], mod)
# fit <- lmFit(gobp_gsva[,c(10:12, 7:9)], mod)
fit <- eBayes(fit)
res <- decideTests(fit, p.value = 0.05)

summary(res) # no significant results

tt <- topTable(fit, coef = 2, n = Inf)
