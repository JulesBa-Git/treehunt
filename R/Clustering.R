#' Clustering of a cocktail list using DBSCAN on umap projection
#' @export
clustering_genetic_algorithm <- function(cocktail_list,
                                         patient_data,
                                         node_column,
                                         target_column,
                                         tree,
                                         depth_column,
                                         upper_bound_column = NULL,
                                         name_column = NULL,
                                         score_type = "wilcoxon",
                                         umap_config= NULL,
                                         eps_dbscan = .2,
                                         min_pts_dbscan = 5){
  requireNamespace("umap")
  requireNamespace("dbscan")
  
  divergence <- get_dissimilarity_of_list(cocktail_list,
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
  umap_config$input <- 'dist'
  umap_results <- umap::umap(divergence, config = umap_config)
  layout <- umap_results$layout
  dbscan_results <- dbscan::dbscan(layout, eps_dbscan, min_pts_dbscan)
  return (data.frame(cocktails = I(cocktail_list),
                     UMAP1 = umap_results$layout[,1],
                     UMAP2 = umap_results$layout[,2],
                     cluster = dbscan_results$cluster))
}