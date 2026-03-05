#' Run Multiple Genetic Algorithm Replicates from Configuration
#'
#' @param config_path Path to a JSON file containing GA hyperparameters.
#' @param patient_data The data frame containing patient/node observations.
#' @param tree_df The data frame defining the tree structure.
#' @param replicates Number of times to run each configuration. Default is 5.
#' @param output_dir Directory where result JSONs will be saved.
#' 
#' @export
run_ga_batch <- function(config_path, 
                         patient_data, 
                         tree_df, 
                         seed_population = NULL,
                         replicates = 5, 
                         output_dir = "results") {
  
  if (!dir.exists(output_dir)) dir.create(output_dir)
  
  # Load configurations
  configs <- jsonlite::fromJSON(config_path)
  
  # Ensure configs is iterable (handles single-row JSONs)
  if (is.list(configs) && !is.data.frame(configs)) {
    configs <- as.data.frame(configs)
  }
  
  for (i in seq_len(nrow(configs))) {
    conf <- configs[i, ]
    
    conf_name <- if ("name" %in% names(conf)) conf$name else paste0("config_", i)
    
    message(sprintf("===> Starting Configuration: %s (%d replicates)", conf_name, replicates))
    
    batch_results <- list()
    
    for (r in seq_len(replicates)) {
      message(sprintf("  Replicate %d/%d...", r, replicates))
      
      res <- run_genetic_algorithm_df_tree(
        patient_data        = patient_data,
        node_column         = conf$node_column,   # Name or Index
        target_column       = conf$target_column, # Name or Index
        tree                = tree_df,
        depth_column        = conf$depth_column,  # Name or Index
        upper_bound_column  = if (!is.null(conf$upper_bound_column)) conf$upper_bound_column else NULL,
        name_column         = if (!is.null(conf$name_column)) conf$name_column else NULL,
        seed_population     = seed_population, 
        population_size     = as.integer(conf$population_size),
        epochs              = as.integer(conf$epochs),
        mutation_rate       = as.numeric(conf$mutation_rate),
        prob_mutation_type1 = as.numeric(conf$prob_mutation_type1),
        crossover_rate      = as.numeric(conf$crossover_rate),
        elite_count         = as.integer(ifelse(is.null(conf$elite_count), 0, conf$elite_count)),
        tournament_size     = as.integer(ifelse(is.null(conf$tournament_size), 3, conf$tournament_size)),
        alpha               = as.numeric(conf$alpha),
        score_type          = as.character(conf$score_type),
        diversity           = as.logical(conf$diversity),
        verbose             = as.logical(conf$verbose)
      )
      
      res$metadata <- list(
        replicate = r,
        config_name = conf_name,
        timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
      )
      
      batch_results[[r]] <- res
    }
    
    output_file <- file.path(output_dir, paste0("results_", conf_name, ".json"))
    jsonlite::write_json(batch_results, output_file, simplifyVector = TRUE, pretty = TRUE)
    message(sprintf("Done! Saved to %s\n", output_file))
  }
}