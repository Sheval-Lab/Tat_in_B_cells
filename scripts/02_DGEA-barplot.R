## Load libraries --------------------------------------------------------------
library(tidyverse)

input_dir <- "output/tables/01_DGEA/deseq"
output_dir <- "output/tables/02_DGEA-plots"
fig_dir <- "output/figures/02_DGEA-plots"


## Load DGEA results -----------------------------------------------------------
tat_vs_lcl <- read_tsv(str_c(input_dir, "Tat_vs_LCL.LFC.DE.tsv", sep = "/"))
cys_vs_lcl <- read_tsv(str_c(input_dir, "Cys_vs_LCL.LFC.DE.tsv", sep = "/"))
gfp_vs_lcl <- read_tsv(str_c(input_dir, "GFP_vs_LCL.LFC.DE.tsv", sep = "/"))


## Calculate DEGs numbers ------------------------------------------------------
degs <- rbind(
  ### Tat vs LCL
  transmute(
    tat_vs_lcl, geneID, log2FC, 
    gene_type = ifelse(gene_type == "protein_coding", "protein_coding", "other"),
    comparison = "tat_vs_lcl"),
  ### Cys vs LCL
  transmute(
    cys_vs_lcl, geneID, log2FC, 
    gene_type = ifelse(gene_type == "protein_coding", "protein_coding", "other"),
    comparison = "cys_vs_lcl"),
  ### GFP vs LCL
  transmute(
    gfp_vs_lcl, geneID, log2FC, 
    gene_type = ifelse(gene_type == "protein_coding", "protein_coding", "other"),
    comparison = "gfp_vs_lcl"))

### Calculate number of all types of DEGs
degs_number_all <- degs %>% 
  mutate(direction = ifelse(log2FC > 0, "pos", "neg")) %>% 
  group_by(comparison, direction) %>% 
  summarise(n = n()) %>% 
  mutate(type = "all") %>% 
  ungroup()

### Calculate number of protein-coding DEGs
degs_number_pc <- degs %>% 
  filter(gene_type == "protein_coding") %>% 
  mutate(direction = ifelse(log2FC > 0, "pos", "neg")) %>% 
  group_by(comparison, direction) %>% 
  summarise(n = n()) %>% 
  mutate(type = "protein_coding") %>% 
  ungroup()

### Combine dataframes with DEGs numbers
degs_number <- rbind(degs_number_all, degs_number_pc) %>% 
  mutate(
    n = ifelse(direction == "neg", -n, n),
    text_position = ifelse(direction == "neg", -100, 100),
    comparison = factor(comparison, 
                        levels = c("gfp_vs_lcl", "tat_vs_lcl", "cys_vs_lcl")))

write_tsv(degs_number, str_c(output_dir, "degs_number.tsv", sep = "/"))

## Plot DEGs number as barplot -------------------------------------------------
degs_number <- read_tsv(str_c(output_dir, "degs_number.tsv", sep = "/")) %>% 
  mutate(comparison = factor(comparison, 
                             levels = c("gfp_vs_lcl", "tat_vs_lcl", "cys_vs_lcl")))

cols <- c("pos" = "red", "neg" = "blue")
facet_names <- c(all = "All genes", protein_coding = "Protein-coding genes")
sample_names_format <- c(
  expression("RPMI"^{"EGFP"}~"vs"~"RPMI"), 
  expression("RPMI"^{"Tat"}~"vs"~"RPMI"),  
  expression("RPMI"^{"Cys"}~"vs"~"RPMI"))
y_axis <- c(-600, -400, -200, 0, 200, 400)

ggplot(degs_number, aes(x = comparison, y = n, fill = direction)) + 
  geom_bar(stat="identity", position="identity", alpha = 0.7) + 
  labs(x="", y="# DEGs") +
  geom_text(data = degs_number, aes(x = comparison, y = text_position, label = abs(n)), size = 3) +
  facet_wrap(. ~ type, labeller = labeller(type = facet_names)) +
  ## Edit color scheme
  scale_fill_manual(name = "", breaks = c("pos", "neg"), 
                    labels = c("Upregulated", "Downregulated"), values = cols) +
  ## Edit sample names
  scale_x_discrete(labels = rep(sample_names_format, times = 2)) +
  ## Edit Y-axis
  scale_y_continuous(breaks = y_axis, labels = abs(y_axis)) +
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


ggsave(str_c(fig_dir, "DEGs_barplot.png", sep = "/"), units = "cm", width = 7, height = 10, scaling = 2/3)
ggsave(str_c(fig_dir, "DEGs_barplot.pdf", sep = "/"), units = "cm", width = 7, height = 10, scaling = 2/3)

