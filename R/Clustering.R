#' Clustering of a cocktail list using DBSCAN on umap projection
#' @export
clustering_genetic_algorithm <- function(cocktail_list,
                                         patient_data,
                                         node_column,
                                         target_column,
                                         tree,
                                         depth_column,
                                         upper_bound_column = R_NilValue,
                                         name_column = R_NilValue,
                                         score_type = "wilcoxon",
                                         umap_config= R_NilValue){
  requireNamespace("umap")
  requireNamespace("dbscan")
  
  divergence <- get_dissimilarity_of_list(genetic_results,
                                          patient_data,
                                          node_column,
                                          target_column,
                                          tree,
                                          depth_column,
                                          upper_bound_column,
                                          name_column,
                                          score_type)
  
  if(is.null(umap_config)){
    umap_config = umap::umap.defaults
  }
  
  umap_results <- umap::umap(divergence, config = umap_config)
  layout <- umap_results$layout
  dbscan_results <- dbscan::dbscan(layout)
  return (data.frame(cocktails = genetic_results,
                     UMAP1 = umap_results$layout[,1],
                     UMAP2 = umap_results$layout[,2],
                     cluster = dbscan_results$cluster))
}