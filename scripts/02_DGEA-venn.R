## Load libraries --------------------------------------------------------------
library(tidyverse)
library(VennDiagram)

input_dir <- "output/tables/01_DGEA/deseq"
output_dir <- "output/tables/02_DGEA-plots"
fig_dir <- "output/figures/02_DGEA-plots"


## Specify log2FoldChange threshold for DEGs -----------------------------------
log2fc_threshold <- log2(1.5)


## Load DGEA results -----------------------------------------------------------
tat_vs_lcl <- read_tsv(str_c(input_dir, "Tat_vs_LCL.LFC.DE.tsv", sep = "/"))
cys_vs_lcl <- read_tsv(str_c(input_dir, "Cys_vs_LCL.LFC.DE.tsv", sep = "/"))


## Prepare list of UP and DOWN-regulated DEGs ----------------------------------
### All genes
tat_cys <- list(
  tat_up = filter(tat_vs_lcl, log2FC >= log2fc_threshold) %>%  pull(geneID), 
  cys_up = filter(cys_vs_lcl, log2FC >= log2fc_threshold) %>%  pull(geneID), 
  tat_dn = filter(tat_vs_lcl, log2FC <= (-log2fc_threshold)) %>%  pull(geneID),
  cys_dn = filter(cys_vs_lcl, log2FC <= (-log2fc_threshold)) %>%  pull(geneID))

### Protein-coding genes
tat_cys_pc <- list(
  tat_up = filter(tat_vs_lcl, log2FC >= log2fc_threshold, gene_type == "protein_coding") %>%  pull(geneID), 
  cys_up = filter(cys_vs_lcl, log2FC >= log2fc_threshold, gene_type == "protein_coding") %>%  pull(geneID), 
  tat_dn = filter(tat_vs_lcl, log2FC <= (-log2fc_threshold), gene_type == "protein_coding") %>%  pull(geneID),
  cys_dn = filter(cys_vs_lcl, log2FC <= (-log2fc_threshold), gene_type == "protein_coding") %>%  pull(geneID))

cat_names <- c(
  expression("RPMI"^{"Tat"}~"vs"~"RPMI"), 
  expression("RPMI"^{"Cys"}~"vs"~"RPMI"))

# cols <- c("#db3a34", "#ef6351", "#f7a399", "#03045e", "#0077b6", "#00b4d8")

cols2 <- c("#F39B7F", "#8491B4")

## Create Venn plot ------------------------------------------------------------
### UP-regulated DEGs ----------------------------------------------------------
#### All genes -----------------------------------------------------------------
#### PNG 
venn.diagram(
  tat_cys[1:2],
  str_c(fig_dir, "UP_venn.png", sep = "/"), 
  imagetype = "png",
  units = "cm", width = 4, height = 4,
  ## Circle's circumference
  # col = cols[c(2, 5)], lwd = c(1, 1),
  col = cols2, lwd = c(1, 1),
  ## Circle's area
  # fill = cols[c(2, 5)], alpha = c(0.5, 0.5),
  fill = cols2, alpha = c(0.5, 0.5),
  ## Category names
  category.names = cat_names,
  cat.fontfamily = "sans",
  cat.cex = c(0.5, 0.5),
  cat.dist = c(0.07, 0.07),
  cat.pos = c(340, 20),
  ## Main title
  main = "Upregulated DEGs", 
  main.cex = 0.6, 
  main.fontfamily = "sans",
  main.pos = c(0.5,0.78),
  ## Set labels
  fontfamily = "sans",
  cex = c(0.7, 0.7, 0.7),
  margin = 0.3)

#### SVG
venn.diagram(
  tat_cys[1:2],
  str_c(fig_dir, "UP_venn.svg", sep = "/"), 
  imagetype = "svg",
  units = "cm", width = 2, height = 2,
  ## Circle's circumference
  # col = cols[c(2, 5)], lwd = c(1, 1),
  col = cols2, lwd = c(1, 1),
  ## Circle's area
  # fill = cols[c(2, 5)], alpha = c(0.5, 0.5),
  fill = cols2, alpha = c(0.5, 0.5),
  ## Category names
  category.names = cat_names,
  cat.fontfamily = "sans",
  cat.cex = c(0.7, 0.7),
  cat.dist = c(0.07, 0.07),
  cat.pos = c(340, 20),
  ## Main title
  main = "Upregulated DEGs", 
  main.cex = 0.8, 
  main.fontfamily = "sans",
  main.pos = c(0.5,0.78),
  ## Set labels
  fontfamily = "sans",
  cex = c(0.8, 0.8, 0.8),
  margin = 0.3)


#### Protein-coding genes ------------------------------------------------------
#### PNG 
venn.diagram(
  tat_cys_pc[1:2],
  str_c(fig_dir, "UP_PC_venn.png", sep = "/"), 
  imagetype = "png",
  units = "cm", width = 4, height = 4,
  ## Circle's circumference
  # col = cols[c(2, 5)], lwd = c(1, 1),
  col = cols2, lwd = c(1, 1),
  ## Circle's area
  # fill = cols[c(2, 5)], alpha = c(0.5, 0.5),
  fill = cols2, alpha = c(0.5, 0.5),
  ## Category names
  category.names = cat_names,
  cat.fontfamily = "sans",
  cat.cex = c(0.5, 0.5),
  cat.dist = c(0.07, 0.07),
  cat.pos = c(20, 340),
  ## Main title
  main = "Upregulated DEGs", 
  main.cex = 0.6, 
  main.fontfamily = "sans",
  main.pos = c(0.5,0.78),
  ## Set labels
  fontfamily = "sans",
  cex = c(0.7, 0.7, 0.7),
  margin = 0.3,
  inverted = TRUE)

#### SVG
venn.diagram(
  tat_cys_pc[1:2],
  str_c(fig_dir, "UP_PC_venn.svg", sep = "/"), 
  imagetype = "svg",
  units = "cm", width = 2, height = 2,
  ## Circle's circumference
  # col = cols[c(2, 5)], lwd = c(1, 1),
  col = cols2, lwd = c(1, 1),
  ## Circle's area
  # fill = cols[c(2, 5)], alpha = c(0.5, 0.5),
  fill = cols2, alpha = c(0.5, 0.5),
  ## Category names
  category.names = cat_names,
  cat.fontfamily = "sans",
  cat.cex = c(0.7, 0.7),
  cat.dist = c(0.07, 0.07),
  cat.pos = c(20, 340),
  ## Main title
  main = "Upregulated DEGs", 
  main.cex = 0.8, 
  main.fontfamily = "sans",
  main.pos = c(0.5,0.78),
  ## Set labels
  fontfamily = "sans",
  cex = c(0.8, 0.8, 0.8),
  margin = 0.3,
  inverted = TRUE)

### DOWN-regulated DEGs --------------------------------------------------------
#### All genes -----------------------------------------------------------------
#### PNG 
venn.diagram(
  tat_cys[3:4],
  str_c(fig_dir, "DOWN_venn.png", sep = "/"), 
  imagetype = "png",
  units = "cm", width = 4, height = 4,
  ## Circle's circumference
  # col = cols[c(1, 4)], lwd = c(1, 1),
  col = cols2, lwd = c(1, 1),
  ## Circle's area
  # fill = cols[c(1, 4)], alpha = c(0.4, 0.4),
  fill = cols2, alpha = c(0.5, 0.5),
  ## Category names
  category.names = cat_names,
  cat.fontfamily = "sans",
  cat.cex = c(0.5, 0.5),
  cat.dist = c(0.06, 0.08),
  cat.pos = c(340, 20),
  ## Main title
  main = "Downregulated DEGs", 
  main.cex = 0.6, 
  main.fontfamily = "sans",
  main.pos = c(0.5,0.78),
  ## Set labels
  fontfamily = "sans",
  cex = c(0.7, 0.7, 0.7),
  margin = 0.3)

#### SVG 
venn.diagram(
  tat_cys[3:4],
  str_c(fig_dir, "DOWN_venn.svg", sep = "/"), 
  imagetype = "svg",
  units = "cm", width = 2, height = 2,
  ## Circle's circumference
  # col = cols[c(1, 4)], lwd = c(1, 1),
  col = cols2, lwd = c(1, 1),
  ## Circle's area
  # fill = cols[c(1, 4)], alpha = c(0.4, 0.4),
  fill = cols2, alpha = c(0.5, 0.5),
  ## Category names
  category.names = cat_names,
  cat.fontfamily = "sans",
  cat.cex = c(0.7, 0.7),
  cat.dist = c(0.06, 0.08),
  cat.pos = c(340, 20),
  ## Main title
  main = "Downregulated DEGs", 
  main.cex = 0.8, 
  main.fontfamily = "sans",
  main.pos = c(0.5,0.78),
  ## Set labels
  fontfamily = "sans",
  cex = c(0.8, 0.8, 0.8),
  margin = 0.3)
  
#### Protein-coding genes ------------------------------------------------------
#### PNG 
venn.diagram(
  tat_cys_pc[3:4],
  str_c(fig_dir, "DOWN_PC_venn.png", sep = "/"), 
  imagetype = "png",
  units = "cm", width = 4, height = 4,
  ## Circle's circumference
  # col = cols[c(1, 4)], lwd = c(1, 1),
  col = cols2, lwd = c(1, 1),
  ## Circle's area
  # fill = cols[c(1, 4)], alpha = c(0.4, 0.4),
  fill = cols2, alpha = c(0.5, 0.5),
  ## Category names
  category.names = cat_names,
  cat.fontfamily = "sans",
  cat.cex = c(0.5, 0.5),
  cat.dist = c(0.06, 0.08),
  cat.pos = c(340, 20),
  ## Main title
  main = "Downregulated DEGs", 
  main.cex = 0.6, 
  main.fontfamily = "sans",
  main.pos = c(0.5,0.78),
  ## Set labels
  fontfamily = "sans",
  cex = c(0.7, 0.7, 0.7),
  margin = 0.3)

#### SVG 
venn.diagram(
  tat_cys_pc[3:4],
  str_c(fig_dir, "DOWN_PC_venn.svg", sep = "/"), 
  imagetype = "svg",
  units = "cm", width = 2, height = 2,
  ## Circle's circumference
  # col = cols[c(1, 4)], lwd = c(1, 1),
  col = cols2, lwd = c(1, 1),
  ## Circle's area
  # fill = cols[c(1, 4)], alpha = c(0.4, 0.4),
  fill = cols2, alpha = c(0.5, 0.5),
  ## Category names
  category.names = cat_names,
  cat.fontfamily = "sans",
  cat.cex = c(0.7, 0.7),
  cat.dist = c(0.06, 0.08),
  cat.pos = c(340, 20),
  ## Main title
  main = "Downregulated DEGs", 
  main.cex = 0.8, 
  main.fontfamily = "sans",
  main.pos = c(0.5,0.78),
  ## Set labels
  fontfamily = "sans",
  cex = c(0.8, 0.8, 0.8),
  margin = 0.3)
