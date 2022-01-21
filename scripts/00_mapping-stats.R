## Load libraries --------------------------------------------------------------
library(tidyverse)


input_dir <- "data/summaries"
output_dir <- "output/tables/01_DGEA/processed_counts"
fig_dir <- "output/figures/01_DGEA"


## Load data -------------------------------------------------------------------
### Import hisat2 summaries summarized :-)
hisat2_summary_file <- str_c(input_dir, "hisat2_summary.tsv", sep = "/")
hisat2_summary <- read_tsv(hisat2_summary_file)


### Import gene count summary
counts_summary_file <- str_c(output_dir, "counts_summary.tsv", sep = "/")
counts_summary <- read_tsv(counts_summary_file)


### Sample names
sample_names_format_all <- c(
  expression("RPMI"^{"Tat"}~1), expression("RPMI"^{"Tat"}~2), expression("RPMI"^{"Tat"}~3),
  expression("RPMI"^{"Cys"}~1), expression("RPMI"^{"Cys"}~2), expression("RPMI"^{"Cys"}~3),  
  expression("RPMI"^{"EGFP"}~1), expression("RPMI"^{"EGFP"}~2), expression("RPMI"^{"EGFP"}~3),
  expression("RPMI"~1), expression("RPMI"~2), expression("RPMI"~3))

sample_names <- c(
  "Tat.1", "Tat.2", "Tat.3",
  "Cys.1", "Cys.2", "Cys.3",
  "GFP.1", 'GFP.2', "GFP.3",
  "LCL.1", "LCL.2", "LCL.3")


## Process summaries -----------------------------------------------------------
## Combine hisat2 and gene count summaries
summary_cmb <- hisat2_summary %>% 
  left_join(counts_summary, by = "sample") %>% 
  mutate(
    `__ambiguous` = one_al - `__no_feature` - `__unique`)

summary_cmb_perc <- summary_cmb %>% 
  mutate(
    # mapping to gene annotation
    `__unique` = (`__unique` / one_al) * 100, 
    `__ambiguous` = (`__ambiguous` / one_al) * 100, 
    `__no_feature` = (`__no_feature` / one_al) * 100, 
    # mapping to genome
    mult_al = (mult_al / total_reads) * 100,
    one_al = (one_al / total_reads) * 100,
    zero_al = (zero_al / total_reads) * 100,
    total_reads = 100)


## Make plots ------------------------------------------------------------------
### Gene count summary
summary_gc_lng_ <- summary_cmb %>% 
  select(sample, `__no_feature`, `__ambiguous`, `__unique`) %>% 
  pivot_longer(cols = -sample, names_to = "read_group", values_to = "count") %>% 
  mutate(read_group = factor(read_group, levels = c("__unique", "__ambiguous", "__no_feature")))

summary_gc_perc_lng <- summary_cmb_perc %>% 
  select(sample, `__no_feature`, `__ambiguous`, `__unique`) %>% 
  pivot_longer(cols = -sample, names_to = "read_group", values_to = "percent") %>% 
  mutate(read_group = factor(read_group, levels = c("__unique", "__ambiguous", "__no_feature")))


summary_gc_lng <- summary_gc_lng_ %>% 
  left_join(summary_gc_perc_lng, by = c("sample", "read_group")) %>% 
  group_by(sample) %>% 
  mutate(total = sum(count)) %>% 
  ungroup() %>% 
  mutate(
    sample = factor(sample, levels = sample_names),
    pos = case_when(
      read_group == "__unique" ~ total * 0.5,
      read_group == "__ambiguous" ~ total * 0.12,
      read_group == "__no_feature" ~ total * 0.025))


cols4 <- c("#CC79A7", "#D55E00", "#009E73", "#0072B2")

ggplot(summary_gc_lng, aes(x = sample, y = count, fill = read_group)) +
  geom_col(width = 0.65, alpha = 0.7) +
  geom_text(aes(x = sample, y = pos, label = paste0(round(percent, 1), "%")), size = 3) +
  scale_x_discrete(labels = sample_names_format_all) +
  scale_y_continuous(labels = scales::format_format(big.mark = " ", decimal.mark = ".", scientific = FALSE)) +
  scale_fill_manual(values = rev(cols4[1:3]), labels = c("Uniquely mapped", "Ambiguously mapped", "Unmapped")) +
  labs(x = "", y = "# reads", fill = "") +
  ## Edit theme
  theme(text = element_text(size = 12, color = "black"), 
        axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
        axis.text.y = element_text(color = "black"),
        panel.background = element_rect(fill = "white", colour = "white", size = 0.5, linetype = "solid"), 
        axis.line.x.bottom = element_line(color = "black"),
        axis.line.y.left = element_line(color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "top",
        aspect.ratio = 2/3)


ggsave(str_c(fig_dir, "gene_counts_stat.png", sep = "/"), units = "cm", width = 15, scaling = 3/4)
ggsave(str_c(fig_dir, "gene_counts_stat.pdf", sep = "/"), units = "cm", width = 15, scale = 4/3)

### Mapping and count summary
summary_mc_lng_ <- summary_cmb %>% 
  select(sample, total_reads, one_al, `__unique`) %>% 
  pivot_longer(cols = -sample, names_to = "read_group", values_to = "count") %>% 
  mutate(read_group = factor(read_group, levels = c("total_reads", "one_al", "__unique")))

summary_mc_perc_lng <- summary_cmb_perc %>% 
  select(sample, total_reads, one_al, `__unique`) %>% 
  pivot_longer(cols = -sample, names_to = "read_group", values_to = "percent") %>% 
  mutate(read_group = factor(read_group, levels = c("total_reads", "one_al", "__unique")))

summary_mc_lng <- summary_mc_lng_ %>% 
  left_join(summary_mc_perc_lng, by = c("sample", "read_group")) %>% 
  mutate(
    sample = factor(sample, levels = sample_names),
    pos = count)


cols4 <- c("#CC79A7", "#D55E00", "#009E73", "#0072B2")

ggplot(summary_mc_lng, aes(x = sample, y = count, fill = read_group)) +
  geom_col(width = 0.65, alpha = 0.7, position = position_dodge()) +
  # geom_text(aes(x = sample, y = pos, label = paste0(round(percent, 1), "%")), size = 3) +
  scale_x_discrete(labels = sample_names_format_all) +
  scale_y_continuous(labels = scales::format_format(big.mark = " ", decimal.mark = ".", scientific = FALSE)) +
  scale_fill_manual(
    values = rev(cols4[1:3]), 
    labels = c("Total reads", "Uniquely aligned to genome", "Uniquely mapped to GENCODE annotation")) +
  labs(x = "", y = "# reads", fill = "") +
  ## Edit theme
  theme(text = element_text(size = 12, color = "black"), 
        axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
        axis.text.y = element_text(color = "black"),
        panel.background = element_rect(fill = "white", colour = "white", size = 0.5, linetype = "solid"), 
        axis.line.x.bottom = element_line(color = "black"),
        axis.line.y.left = element_line(color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "top",
        aspect.ratio = 2/3)

ggsave(str_c(fig_dir, "S1_gene_map_counts_stat.png", sep = "/"), units = "cm", width = 15, scaling = 3/4)
ggsave(str_c(fig_dir, "S1_gene_map_counts_stat.pdf", sep = "/"), units = "cm", width = 15, scale = 4/3)


