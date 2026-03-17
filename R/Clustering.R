#' Clustering of a Cocktail List using DBSCAN on UMAP Projection
#'
#' Projects cocktail dissimilarities into 2D space using UMAP and identifies
#' clusters using DBSCAN.
#'
#' @param cocktails A list of 1-based integer vectors or a data frame containing them.
#' @param umap_config Optional configuration list for UMAP. Defaults to `umap.defaults`.
#' @param eps_dbscan Distance parameter for DBSCAN. Default is 0.2.
#' @param min_pts_dbscan Minimum points to form a cluster in DBSCAN. Default is 5.
#' @inheritParams process_ga_scores
#' @export
clustering_genetic_algorithm <- function(cocktails,
                                         patient_data,
                                         node_column,
                                         target_column,
                                         tree_df,
                                         depth_column,
                                         upper_bound_column = NULL,
                                         name_column = NULL,
                                         score_type = "wilcoxon",
                                         umap_config = NULL,
                                         eps_dbscan = 0.2,
                                         min_pts_dbscan = 5) {
  
  requireNamespace("umap", quietly = TRUE)
  requireNamespace("dbscan", quietly = TRUE)
  
  # Work with either a dataframe outputed
  working_list <- if (is.data.frame(cocktails)) {
    if ("vec_cocktails" %in% names(cocktails)) cocktails$vec_cocktail else Filter(is.list, cocktails)[[1]]
  } else {
    cocktails
  }
  
  # Calculate Dissimilarity Matrix
  divergence <- get_dissimilarity_of_list(
    working_list,
    patient_data,
    node_column,
    target_column,
    tree_df,
    depth_column,
    upper_bound_column,
    name_column,
    score_type
  )
  
  # Dimensionality Reduction
  if (is.null(umap_config)) umap_config <- umap::umap.defaults
  umap_config$input <- 'dist'
  
  umap_results <- umap::umap(divergence, config = umap_config)
  layout <- umap_results$layout
  
  # Clustering
  dbscan_results <- dbscan::dbscan(layout, eps = eps_dbscan, minPts = min_pts_dbscan)
  
  # Prepare Output
  res_df <- data.frame(
    UMAP1 = layout[,1],
    UMAP2 = layout[,2],
    cluster = dbscan_results$cluster
  )
  
  # If input was a dataframe, merge. Otherwise, return with cocktails list.
  if (is.data.frame(cocktails)) {
    return(cbind(cocktails, res_df))
  } else {
    return(cbind(data.frame(cocktails = I(working_list)), res_df))
  }
}