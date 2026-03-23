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

#' Plot GA Clustering Results
#'
#' This wrapper runs the clustering pipeline and produces a UMAP scatter plot 
#' where points are colored by their identified cluster.
#'
#' @param cocktails A list of 1-based integer vectors or a data frame containing them.
#' @param ... Arguments passed directly to \code{\link{clustering_genetic_algorithm}}.
#' @param point_size Numeric. Size of the points in the scatter plot. Default is 2.
#' @param alpha Numeric. Transparency of the points. Default is 0.7.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{data}: The data frame with UMAP coordinates and cluster assignments.
#'   \item \code{plot}: A ggplot2 object.
#' }
#' @export
plot_ga_clusters <- function(cocktails, 
                             point_size = 2, 
                             alpha = 0.7, 
                             ...) {
  
  requireNamespace("ggplot2", quietly = TRUE)
  
  # Run the clustering logic
  clustered_data <- clustering_genetic_algorithm(cocktails = cocktails, ...)
  clustered_data$cluster <- as.factor(clustered_data$cluster)
  
  # Plot
  p <- ggplot2::ggplot(clustered_data, 
                       ggplot2::aes(x = UMAP1, y = UMAP2, color = cluster)) +
    ggplot2::geom_point(size = point_size, alpha = alpha) +
    ggplot2::theme_minimal() +
    #ggplot2::scale_color_viridis_d(option = "plasma") +
    ggplot2::labs(
      title = "Genetic Algorithm Cocktail Clusters",
      subtitle = paste("UMAP Projection with DBSCAN Clustering"),
      x = "UMAP Dimension 1",
      y = "UMAP Dimension 2",
      color = "Cluster"
    ) +
    ggplot2::theme(
      legend.position = "right",
      plot.title = ggplot2::element_text(face = "bold")
    )
  
  return(list(data = clustered_data, plot = p))
}