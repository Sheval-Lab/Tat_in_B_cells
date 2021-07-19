# exclude_nested_clusters function
exclude_nested_clusters <- function(result_object, max_unique_genes){
  df <- result_object@result %>% 
    mutate(
      genes = str_split(geneID, "/"),
      cluster_size = length(genes)) %>% 
    arrange(cluster_size) 
    
  
  ## Save indexes of clusters that are 'nested' into other clusters
  clusters2exclude <- c()
  
  for (i in 1:nrow(df)){
    combined <- c(df[-c(clusters2exclude, i), "genes"]) %>% flatten_chr()
    
    cluster_i <- df[i, "genes"] %>% flatten_chr()
    cluster_size <- length(cluster_i)
    
    size_of_nested_subcluster <- cluster_i %>% magrittr::is_in(combined) %>% sum()
    
    if ((cluster_size - size_of_nested_subcluster) <= max_unique_genes){
      clusters2exclude <- append(clusters2exclude, as.integer(i))
    }
  }
  
  result_object@result <- df[-clusters2exclude,] %>% dplyr::select(-genes:cluster_size)
  
  return(result_object)
}


