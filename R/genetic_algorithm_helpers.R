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

#' Aggregate GA Results from Multiple JSON Files returned by run_ga_batch
#'
#' @param folder_path Path to the folder containing result JSON files.
#' @return A data frame of unique cocktails with their scores.
aggregate_ga_results <- function(folder_path = "results") {
  
  files <- list.files(path = folder_path, pattern = "\\.json$", full.names = TRUE)
  
  if (length(files) == 0) stop("No JSON files found in the specified folder.")
  
  # process each file
  all_data <- map_df(files, function(file) {
    
    raw_data <- fromJSON(file)
    
    
    replicate_list <- lapply(seq_len(nrow(raw_data)), function(i) {
      
      cocktails <- raw_data$final_population[[i]]
      scores    <- raw_data$final_scores[[i]] 
      
      cocktail_strings <- sapply(cocktails, function(x) paste(sort(x), collapse = ","))
      
      data.frame(
        cocktail = cocktail_strings,
        score = scores,
        config = raw_data$metadata$config_name[[i]],
        stringsAsFactors = FALSE
      )
    })
    
    bind_rows(replicate_list)
  })
  
  # aggregate unique cocktails
  final_summary <- all_data %>%
    group_by(cocktail) %>%
    summarise(
      score = first(score),
      occurrence_count = n(),
      found_in_configs = paste(unique(config), collapse = "; ")
    ) %>%
    arrange(desc(score))
  
  return(final_summary)
}

#' Map Cocktail Indices to Names
#' TODO : Results of GA are 0 based index for now, they will be change in a near future
#' when this is done, remove the cpp to R index conversion
#' 
#' @param aggregated_df The data frame returned by aggregate_ga_results.
#' @param new_tree The data frame containing the tree structure and $Name column.
#' @return The same data frame with a new 'cocktail_names' column.
map_cocktail_names <- function(aggregated_df, new_tree) {
  
  node_names <- as.character(new_tree$Name)
  
  aggregated_df %>%
    mutate(cocktail_names = sapply(cocktail, function(c_string) {
      # split string "0,4,10" into individual characters
      cpp_indices <- as.numeric(unlist(strsplit(c_string, ",")))
      
      # convert cpp (0-based) to R (1-based) index
      r_indices <- cpp_indices + 1
      
      # safety check
      r_indices <- r_indices[r_indices > 0 & r_indices <= length(node_names)]
      
      # extract names and join them back
      names_vector <- node_names[r_indices]
      return(paste(names_vector, collapse = " + "))
    })) %>%
    # Move the new column to a more visible position
    select(cocktail, cocktail_names, everything())
}

#' Process and Attach GA Scores and Statistics (Generic)
#'
#' @param df The data frame from map_cocktail_names.
#' @param patient_data The dataset used for scoring.
#' @param tree_df The tree structure data frame.
#' @param node_id_col Column name in patient_data representing node indices.
#' @param target_col Column name for the metric (e.g., "QT_c").
#' @param depth_col Column name for tree depth.
#' @param upper_bound_col Column name for tree upper bounds.
#' @param score_type The scoring method (e.g., "Wilcoxon", "RR", "phyper").
#' @param ... Additional arguments passed to compute_score.
process_ga_scores <- function(df, 
                              patient_data, 
                              tree_df, 
                              node_column,
                              target_column,
                              depth_column,
                              upper_bound_column,
                              score_type = "Wilcoxon",
                              ...) {
  
  # prepare indices (0-based to 1-based conversion) TODO : remove this once ga
  # handle the return of 1-based index
  indices_list <- lapply(strsplit(df$cocktail, ","), function(x) as.integer(x) + 1)
  
  # compute_score
  raw_scores <- compute_score(
    indices_list, 
    patient_data, 
    node_column, 
    target_column, 
    tree_df, 
    depth_column, 
    upper_bound_column,
    score_type = score_type,
    ...
  )
  
  # extract and calculate statistics
  df$taker_count <- raw_scores$`number of takers`
  df$scores      <- unlist(raw_scores$scores)
  
  # Helper for safe stat calculation
  get_stat <- function(dist, func) {
    sapply(dist, function(v) {
      if (length(v) == 0 || all(is.na(v))) return(0)
      func(v)
    })
  }
  
  dist <- raw_scores$QT_diff_distribution
  df$QT_mean   <- get_stat(dist, mean)
  df$QT_median <- get_stat(dist, median)
  df$QT_min    <- get_stat(dist, min)
  df$QT_max    <- get_stat(dist, max)
  
  return(df)
}

#' Complete GA Results Analysis Pipeline (Generic)
#'
#' @param folder_path Path to JSON results.
#' @param patient_data The data frame containing patient observations.
#' @param tree_df The tree structure data frame.
#' @param min_score Minimum score threshold to keep a cocktail. Default 0.
#' @param ... Arguments passed down to process_ga_scores (and then to compute_score).
#' @export
summarize_ga_pipeline <- function(folder_path, 
                                  patient_data, 
                                  tree_df, 
                                  min_score = 0,
                                  ...) {
  
  message("Step 1: Aggregating JSON files...")
  results <- aggregate_ga_results(folder_path)
  
  message("Step 2: Mapping cocktail names...")
  results <- map_cocktail_names(results, tree_df)
  
  # Filter early
  results <- results[results$score > min_score, ]
  
  if (nrow(results) == 0) {
    warning("No cocktails found above the minimum score threshold.")
    return(results)
  }
  
  message("Step 3: Computing detailed scores and distributions...")
  results <- process_ga_scores(results, patient_data, tree_df, ...)
  
  message("Success!")
  return(results)
}