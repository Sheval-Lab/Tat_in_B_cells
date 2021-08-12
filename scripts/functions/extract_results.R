#' extract_results function extracts results of DESeq2 DGEA:
#' transforms the results into a dataframe,
#' adds gene annotation (gene_name, gene_type),
#' and adds normalized counts.
#' It saves extracted results (with and without LFC shrinkage)
#' into 'deseq' subdirectory of the provided output directory 

extract_results <- function(dds, numerator, denominator, gene_annotation, log2FC_threshold){
  ### Extract normalized counts
  counts_norm <- counts(dds, normalized=TRUE) %>% 
    as.data.frame() %>% 
    mutate(geneID = rownames(.)) %>% 
    select(geneID, everything())
  
  samples2extract <- 
    c(TRUE,
      colnames(counts(dds, normalized=TRUE)) %>% 
        str_detect(str_c(as.character(numerator), as.character(denominator), sep = "|")))
  
  counts_norm <- counts_norm[, samples2extract]
  
  ### Extract results (without LFC shrinkage)
  res <- results(dds, contrast = c("sample", as.character(numerator), as.character(denominator)), alpha = 0.05)
  res_df <- res %>% as.data.frame() %>% 
    mutate(
      geneID = rownames(.), 
      ensID = str_remove(geneID, ".\\d+$"),
      FC = 2^log2FoldChange) %>% 
    filter(!is.na(padj)) %>%
    arrange(padj) %>%
    left_join(., gene_annotation, by = "geneID") %>% 
    left_join(., counts_norm, by = "geneID") %>% 
    select(geneID, ensID, gene_name, gene_type,
           padj, log2FC = log2FoldChange, FC, 
           everything())
  
  ### Extract results (with LFC shrinkage)
  res_LFC <- lfcShrink(dds, coef = str_c("sample", as.character(numerator), "vs", as.character(denominator), sep = "_"), 
                       res = res, type="apeglm")
  res_LFC_df <- res_LFC %>% as.data.frame() %>% 
    mutate(
      geneID = rownames(.), 
      ensID = str_remove(geneID, ".\\d+$"),
      FC = 2^log2FoldChange) %>% 
    filter(!is.na(padj)) %>%
    arrange(padj) %>%
    left_join(., gene_annotation, by = "geneID") %>% 
    left_join(., counts_norm, by = "geneID") %>% 
    select(geneID, ensID, gene_name, gene_type,
           padj, log2FC = log2FoldChange, FC, 
           everything())
  
  ### Save tables
  res_file <- str_c(as.character(numerator), "vs", as.character(denominator), sep = "_") %>% 
    str_c(., ".tsv") %>% 
    str_c(output_dir, "deseq", ., sep = "/")
  
  res_DE_file <- str_c(as.character(numerator), "vs", as.character(denominator), sep = "_") %>% 
    str_c(., ".DE.tsv") %>% 
    str_c(output_dir, "deseq", ., sep = "/")
  
  res_LFC_file <- str_c(as.character(numerator), "vs", as.character(denominator), sep = "_") %>% 
    str_c(., ".LFC.tsv") %>% 
    str_c(output_dir, "deseq", ., sep = "/")
  
  res_LFC_DE_file <- str_c(as.character(numerator), "vs", as.character(denominator), sep = "_") %>% 
    str_c(., ".LFC.DE.tsv") %>% 
    str_c(output_dir, "deseq", ., sep = "/")
  
  write_tsv(res_df, res_file)
  write_tsv(res_LFC_df, res_LFC_file)
  
  res_df %>% filter(padj < 0.05, abs(log2FC) >= as.double(log2FC_threshold)) %>% 
    write_tsv(., res_DE_file)
  
  res_LFC_df %>% filter(padj < 0.05, abs(log2FC) >= as.double(log2FC_threshold)) %>% 
    write_tsv(., res_LFC_DE_file)
  
}