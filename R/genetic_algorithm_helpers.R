#' Run Multiple Genetic Algorithm Replicates from Configuration
#'
#' @param config_path Path to a JSON file containing GA hyperparameters.
#' @param patient_data The data frame containing patient/node observations.
#' @param tree_df The data frame defining the tree structure.
#' @param seed_population Optional initial combinations, using one-based tree indices.
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
        id_column           = if (!is.null(conf$id_column)) conf$id_column else NULL,
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
#' @export
aggregate_ga_results <- function(folder_path = "results") {
  files <- list.files(folder_path, pattern = "\\.json$", full.names = TRUE)
  if (!length(files)) stop("No JSON files found in the specified folder.")
  empty <- data.frame(cocktail = character(), score = numeric(),
                      occurrence_count = integer(), found_in_configs = character())
  rows <- unlist(lapply(files, function(file) {
    runs <- jsonlite::fromJSON(file, simplifyVector = FALSE)
    lapply(runs, function(run) {
      cocktails <- run$final_population
      scores <- as.numeric(unlist(run$final_scores, use.names = FALSE))
      if (length(cocktails) != length(scores)) {
        stop("Each saved combination must have one score in ", basename(file), call. = FALSE)
      }
      if (!length(cocktails)) return(NULL)
      keys <- vapply(cocktails, function(x) paste(sort(unlist(x)), collapse = ","), character(1))
      config <- run$metadata$config_name
      if (is.null(config)) config <- NA_character_
      config <- as.character(unlist(config, use.names = FALSE))[1L]
      data.frame(cocktail = keys, score = scores, config = config)
    })
  }), recursive = FALSE)
  all_data <- do.call(rbind, rows)
  if (is.null(all_data) || !nrow(all_data)) return(empty)
  keys <- sort(unique(all_data$cocktail))
  summary <- do.call(rbind, lapply(keys, function(key) {
    group <- all_data[all_data$cocktail == key, , drop = FALSE]
    data.frame(cocktail = key, score = group$score[1L],
               occurrence_count = nrow(group),
               found_in_configs = paste(unique(group$config), collapse = "; "))
  }))
  summary <- summary[order(-summary$score, summary$cocktail), , drop = FALSE]
  rownames(summary) <- NULL
  summary
}

#' Map Cocktail Indices to Names
#' Saved GA combinations use zero-based node indices; tree rows use one-based indices.
#' 
#' @param aggregated_df The data frame returned by aggregate_ga_results.
#' @param new_tree Tree data frame with `Name` and `Code` columns.
#' @return The data frame with `cocktail_names` and `cocktail_codes` columns.
#' @export
map_cocktail_names <- function(aggregated_df, new_tree) {
  node_names <- as.character(new_tree$Name)
  node_codes <- as.character(new_tree$Code)
  indices <- lapply(strsplit(aggregated_df$cocktail, ","), function(x) {
    rows <- as.integer(x) + 1L
    rows[!is.na(rows) & rows > 0 & rows <= length(node_names)]
  })
  aggregated_df$cocktail_names <- vapply(indices, function(rows) {
    paste(node_names[rows], collapse = " | ")
  }, character(1))
  aggregated_df$cocktail_codes <- vapply(indices, function(rows) {
    paste(node_codes[rows], collapse = " | ")
  }, character(1))
  first <- c("cocktail", "cocktail_names", "cocktail_codes")
  aggregated_df[, c(first, setdiff(names(aggregated_df), first)), drop = FALSE]
}

#' Process and Attach GA Scores and Statistics (Generic)
#'
#' @param df The data frame from map_cocktail_names.
#' @param patient_data The dataset used for scoring.
#' @param tree_df The tree structure data frame.
#' @param node_column Column in `patient_data` containing zero-based node indices.
#' @param target_column Column containing the outcome (for example, `"QT_c"`).
#' @param depth_column Column in `tree_df` containing node depths.
#' @param upper_bound_column Column in `tree_df` containing zero-based subtree bounds.
#' @param score_type The scoring method (e.g., "Wilcoxon", "RR", "phyper").
#' @param ... Additional arguments passed to compute_score.
#' @param id_column Optional observation-unit identifier required by patient-level
#'   continuous-outcome scores.
#' @param name_column Optional node-label column in `tree_df`.
#' @export
process_ga_scores <- function(df, 
                              patient_data, 
                              tree_df, 
                              node_column,
                              target_column,
                              depth_column,
                              upper_bound_column,
                              score_type = "composite",
                              ...,
                              id_column = NULL,
                              name_column = NULL) {
  
  # The GA returns zero-based candidates, whereas compute_score() deliberately
  # accepts one-based R row indices at its public boundary.
  indices_list <- lapply(strsplit(df$cocktail, ","), function(x) as.integer(x) + 1)
  
  # compute_score
  raw_scores <- compute_score(
    cocktail_list = indices_list,
    patient_data = patient_data,
    node_column = node_column,
    target_column = target_column,
    tree = tree_df,
    depth_column = depth_column,
    id_column = id_column,
    upper_bound_column = upper_bound_column,
    name_column = name_column,
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
  df$QT_median <- get_stat(dist, stats::median)
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

#' Filter a list of cocktails accoring to the mean depth of nodes inside 
#' each cocktail. A cocktail is kept if the mean depth of his nodes is greater 
#' than \code{mean_depth_cocktail} (default = 3). This function aim is to provide
#' more interpretable cocktail to the user.
#' 
#' @param results results returned by function \code{summarize_ga_pipeline}. Or a 
#' list of string of 0 indexed vector of tree nodes. (If nodes are 1-indexed, make sure that the parameter
#' one_index is set to TRUE).
#' @param tree_depth The tree structure depth.
#' @param mean_depth_cocktail The desired minimal mean depth of nodes in a cocktail.
#' @param one_index A boolean telling if the cocktails are 1-indexed or 0-indexed. Default is FALSE
#' 
#' @return The filtered list of cocktail
#' @export
filter_out_cocktails <- function(results, tree_depth, mean_depth_cocktail = 3,
                                 one_index = FALSE) {
  if (is.character(results)) results <- data.frame(cocktail = results)
  indices <- lapply(strsplit(results$cocktail, ","), as.integer)
  if (!one_index) indices <- lapply(indices, function(x) x + 1L)
  results$tmp_idx <- indices
  depths <- vapply(indices, function(x) mean(tree_depth[x]), numeric(1))
  results[!is.na(depths) & depths >= mean_depth_cocktail, , drop = FALSE]
}
