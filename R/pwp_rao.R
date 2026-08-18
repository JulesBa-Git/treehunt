#' Fit the reduced PWP model and prepare a Rao-score context
#'
#' The reduced Cox model is fitted once. Its linear predictor, nuisance design
#' matrix and model-based nuisance variance are stored in a compact context used
#' by the optimized C++ score engine. Patient-clustered variance is retained in
#' the fitted model for reporting, while the high-throughput Rao statistic uses
#' the exact model-based efficient information.
#'
#' @param data Interval-level data frame.
#' @param status_column Binary endpoint column.
#' @param time_column Gap-time column.
#' @param strata_column PWP event-order stratum.
#' @param id_column Patient identifier used for robust clustering.
#' @param covariates Character vector of reduced-model covariates.
#' @param ties Cox tied-time method. The optimized engine currently requires
#'   `"efron"`.
#' @param positive_only If `TRUE`, the search fitness is `max(0, signed_z)`;
#'   full signed statistics are still returned during detailed evaluation.
#' @param min_covered_patients Minimum number of distinct exposed patients for
#'   a candidate to receive non-zero search fitness.
#' @param min_covered_events Minimum number of exposed endpoint events for a
#'   candidate to receive non-zero search fitness.
#' @return An object of class `pwp_rao_context`.
#' @export
fit_pwp_rao_context <- function(
    data,
    status_column,
    time_column = "gap_time",
    strata_column = "event_stratum_25",
    id_column = "id",
    covariates = c("age_start", "index_emergency", "year_start"),
    ties = "efron",
    positive_only = TRUE,
    min_covered_patients = 1L,
    min_covered_events = 1L) {
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("The 'survival' package is required.", call. = FALSE)
  }
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame.", call. = FALSE)
  }
  if (!identical(tolower(ties), "efron")) {
    stop("The optimized PWP Rao engine currently supports ties = 'efron' only.",
         call. = FALSE)
  }
  min_covered_patients <- as.integer(min_covered_patients)
  min_covered_events <- as.integer(min_covered_events)
  if (length(min_covered_patients) != 1L || is.na(min_covered_patients) ||
      min_covered_patients < 1L || length(min_covered_events) != 1L ||
      is.na(min_covered_events) || min_covered_events < 1L) {
    stop("Coverage thresholds must be single positive integers.", call. = FALSE)
  }

  required <- unique(c(
    time_column, status_column, strata_column, id_column, covariates
  ))
  missing <- setdiff(required, names(data))
  if (length(missing)) {
    stop("Missing PWP columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  if (any(make.names(required) != required)) {
    stop("PWP model column names must be syntactically valid R names.", call. = FALSE)
  }

  columns <- lapply(required, function(name) data[[name]])
  names(columns) <- required
  keep <- Reduce(`&`, lapply(columns, function(x) !is.na(x)))
  numeric_columns <- unique(c(time_column, status_column, strata_column,
                              id_column, covariates))
  for (name in numeric_columns) {
    if (!is.numeric(data[[name]]) && !is.integer(data[[name]])) {
      stop("PWP column '", name, "' must be numeric or integer.", call. = FALSE)
    }
    keep <- keep & is.finite(as.numeric(data[[name]]))
  }
  row_index <- which(keep)
  if (!length(row_index)) {
    stop("No complete PWP analysis row remains.", call. = FALSE)
  }

  model_data <- as.data.frame(
    lapply(required, function(name) data[[name]][row_index]),
    optional = TRUE
  )
  names(model_data) <- required
  model_data[[status_column]] <- as.integer(model_data[[status_column]])
  model_data[[strata_column]] <- as.integer(model_data[[strata_column]])
  model_data[[id_column]] <- as.integer(model_data[[id_column]])

  if (any(!model_data[[status_column]] %in% 0:1)) {
    stop("The endpoint must contain only 0 and 1.", call. = FALSE)
  }
  if (any(model_data[[time_column]] <= 0)) {
    stop("All PWP gap times must be strictly positive.", call. = FALSE)
  }
  if (!sum(model_data[[status_column]])) {
    stop("The selected endpoint contains no event.", call. = FALSE)
  }

  rhs <- c(
    covariates,
    sprintf("strata(%s)", strata_column),
    sprintf("cluster(%s)", id_column)
  )
  formula <- stats::as.formula(sprintf(
    "Surv(%s, %s) ~ %s",
    time_column, status_column, paste(rhs, collapse = " + ")
  ))
  environment(formula) <- asNamespace("survival")

  fit <- survival::coxph(
    formula = formula,
    data = model_data,
    ties = "efron",
    na.action = stats::na.fail,
    x = TRUE,
    y = FALSE,
    model = FALSE
  )
  if (any(!is.finite(stats::coef(fit)))) {
    stop("The reduced PWP model produced non-finite coefficients.", call. = FALSE)
  }

  nuisance_variance <- if (!is.null(fit$naive.var)) fit$naive.var else fit$var
  x <- unclass(fit$x)
  if (is.null(dim(x))) {
    x <- matrix(x, ncol = 1L)
  }
  storage.mode(x) <- "double"
  nuisance_variance <- as.matrix(nuisance_variance)
  storage.mode(nuisance_variance) <- "double"

  score_data <- list(
    time = as.numeric(model_data[[time_column]]),
    status = as.integer(model_data[[status_column]]),
    strata = as.integer(model_data[[strata_column]]),
    id = as.integer(model_data[[id_column]]),
    risk = exp(as.numeric(fit$linear.predictors)),
    x = x,
    nuisance_variance = nuisance_variance,
    positive_only = isTRUE(positive_only),
    min_covered_patients = min_covered_patients,
    min_covered_events = min_covered_events
  )

  compact_fit <- fit
  compact_fit$x <- NULL
  compact_fit$y <- NULL
  compact_fit$model <- NULL

  out <- list(
    fit = compact_fit,
    score_data = score_data,
    row_index = row_index,
    specification = list(
      formula = formula,
      status_column = status_column,
      time_column = time_column,
      strata_column = strata_column,
      id_column = id_column,
      covariates = covariates,
      ties = "efron",
      positive_only = isTRUE(positive_only),
      min_covered_patients = min_covered_patients,
      min_covered_events = min_covered_events,
      original_rows = nrow(data),
      analysis_rows = length(row_index),
      omitted_rows = nrow(data) - length(row_index),
      events = sum(model_data[[status_column]])
    )
  )
  class(out) <- "pwp_rao_context"
  out
}

#' @rdname fit_pwp_rao_context
#' @param x A `pwp_rao_context` to print.
#' @param ... Additional print arguments (currently ignored).
#' @export
print.pwp_rao_context <- function(x, ...) {
  spec <- x$specification
  cat("PWP Rao score context\n")
  cat("  Endpoint:", spec$status_column, "\n")
  cat("  Analysis rows:", format(spec$analysis_rows, big.mark = ","), "\n")
  cat("  Events:", format(spec$events, big.mark = ","), "\n")
  cat("  Omitted rows:", format(spec$omitted_rows, big.mark = ","), "\n")
  cat("  Ties:", spec$ties, "\n")
  cat("  Minimum coverage:", spec$min_covered_patients, "patients and",
      spec$min_covered_events, "events\n")
  invisible(x)
}

.pwp_analysis_data <- function(data, context) {
  if (!inherits(context, "pwp_rao_context")) {
    stop("'context' must be produced by fit_pwp_rao_context().", call. = FALSE)
  }
  if (nrow(data) != context$specification$original_rows) {
    stop("The data no longer have the row count used to fit the PWP context.",
         call. = FALSE)
  }
  # The explicit comma also works when a serialized data.table is loaded before
  # its namespace (and therefore its `[` method) has been registered.
  analysis_data <- data[context$row_index, , drop = FALSE]
  spec <- context$specification
  identity_fields <- list(
    time = as.numeric(analysis_data[[spec$time_column]]),
    status = as.integer(analysis_data[[spec$status_column]]),
    strata = as.integer(analysis_data[[spec$strata_column]]),
    id = as.integer(analysis_data[[spec$id_column]])
  )
  for (field in names(identity_fields)) {
    if (!identical(identity_fields[[field]], context$score_data[[field]])) {
      stop(
        "The PWP data rows or values no longer match the fitted context (",
        field, "). Refit the context on these data.", call. = FALSE
      )
    }
  }
  analysis_data
}

.pwp_tree <- function(tree, depth_column, upper_bound_column, name_column) {
  tree <- as.data.frame(tree, stringsAsFactors = FALSE)
  required <- c(depth_column, upper_bound_column, name_column)
  required <- required[!vapply(required, is.null, logical(1))]
  missing <- setdiff(required, names(tree))
  if (length(missing)) {
    stop("Missing tree columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  tree[[depth_column]] <- as.integer(tree[[depth_column]])
  if (!is.null(upper_bound_column)) {
    tree[[upper_bound_column]] <- as.integer(tree[[upper_bound_column]])
    if (any(tree[[upper_bound_column]] < 0L |
            tree[[upper_bound_column]] >= nrow(tree))) {
      stop("Tree upper bounds must be 0-based indexes in [0, nrow(tree)-1].",
           call. = FALSE)
    }
  }
  tree
}

#' Run the genetic algorithm with the PWP Rao fitness
#'
#' @param data Original interval-level data used to fit `context`.
#' @param context A `pwp_rao_context`.
#' @param tree Tree data frame.
#' @param node_column List-column of 0-based tree indexes.
#' @param depth_column,upper_bound_column,name_column Tree column names.
#' @param seed_population Optional list of 1-based node vectors, following the
#'   historical treehunt seed convention.
#' @param population_size Number of solutions per generation.
#' @param epochs Number of generations.
#' @param mutation_rate Probability of mutation.
#' @param prob_mutation_type1 Probability of the first mutation operator.
#' @param crossover_rate Probability of crossover.
#' @param elite_count Number of elite solutions copied to the next generation.
#' @param tournament_size Number of solutions sampled per tournament.
#' @param alpha Add/remove mutation bias. Higher values favor adding nodes and
#'   therefore tend to produce longer combinations.
#' @param diversity Whether to apply the engine's diversity rule.
#' @param seed Integer seed for the C++ random-number generator.
#' @param verbose Whether to print progress.
#' @export
run_pwp_genetic_algorithm <- function(
    data,
    context,
    tree,
    node_column = "icd10_indices",
    depth_column = "Depth",
    upper_bound_column = "UpperBound",
    name_column = "Code",
    seed_population = NULL,
    population_size = 100L,
    epochs = 1000L,
    mutation_rate = 0.1,
    prob_mutation_type1 = 0.2,
    crossover_rate = 0.8,
    elite_count = 2L,
    tournament_size = 3L,
    alpha = 1,
    diversity = FALSE,
    seed = 1L,
    verbose = FALSE) {
  analysis_data <- .pwp_analysis_data(data, context)
  tree <- .pwp_tree(tree, depth_column, upper_bound_column, name_column)
  run_genetic_algorithm_pwp_cpp(
    patient_data = analysis_data,
    node_column = node_column,
    target_column = context$specification$status_column,
    id_column = context$specification$id_column,
    tree = tree,
    depth_column = depth_column,
    upper_bound_column = upper_bound_column,
    name_column = name_column,
    score_context = context$score_data,
    seed_population = seed_population,
    population_size = as.integer(population_size),
    epochs = as.integer(epochs),
    mutation_rate = mutation_rate,
    prob_mutation_type1 = prob_mutation_type1,
    crossover_rate = crossover_rate,
    elite_count = as.integer(elite_count),
    tournament_size = as.integer(tournament_size),
    alpha = alpha,
    diversity = diversity,
    seed = as.integer(seed),
    verbose = verbose
  )
}

#' Run MCMC with the PWP Rao fitness
#'
#' @inheritParams run_pwp_genetic_algorithm
#' @param epochs Number of MCMC iterations.
#' @param temperature Metropolis-Hastings temperature.
#' @param n_results Maximum number of retained top solutions.
#' @param cocktail_size Initial combination size.
#' @param prob_type1 Probability of the first proposal operator.
#' @param beta Minimum number of distinct covered patients for filtered output.
#' @param max_score Upper limit used by score-distribution bins.
#' @export
run_pwp_mcmc <- function(
    data,
    context,
    tree,
    node_column = "icd10_indices",
    depth_column = "Depth",
    upper_bound_column = "UpperBound",
    name_column = "Code",
    epochs = 100000L,
    temperature = 1,
    n_results = 20L,
    cocktail_size = 2L,
    prob_type1 = 0.01,
    beta = 20L,
    max_score = 50,
    seed = 1L,
    verbose = FALSE) {
  analysis_data <- .pwp_analysis_data(data, context)
  tree <- .pwp_tree(tree, depth_column, upper_bound_column, name_column)
  run_mcmc_pwp_cpp(
    patient_data = analysis_data,
    node_column = node_column,
    target_column = context$specification$status_column,
    id_column = context$specification$id_column,
    tree = tree,
    depth_column = depth_column,
    upper_bound_column = upper_bound_column,
    name_column = name_column,
    score_context = context$score_data,
    epochs = as.integer(epochs),
    temperature = temperature,
    n_results = as.integer(n_results),
    cocktail_size = as.integer(cocktail_size),
    prob_type1 = prob_type1,
    beta = as.integer(beta),
    max_score = max_score,
    seed = as.integer(seed),
    verbose = verbose
  )
}

#' Evaluate full Rao statistics for combinations
#'
#' @inheritParams run_pwp_genetic_algorithm
#' @param combinations List of node-index vectors.
#' @param index_base Either 0 or 1.
#' @export
score_pwp_combinations <- function(
    combinations,
    data,
    context,
    tree,
    index_base = 0L,
    node_column = "icd10_indices",
    depth_column = "Depth",
    upper_bound_column = "UpperBound",
    name_column = "Code") {
  analysis_data <- .pwp_analysis_data(data, context)
  tree <- .pwp_tree(tree, depth_column, upper_bound_column, name_column)
  details <- score_pwp_combinations_cpp(
    combinations = combinations,
    index_base = as.integer(index_base),
    patient_data = analysis_data,
    node_column = node_column,
    target_column = context$specification$status_column,
    id_column = context$specification$id_column,
    tree = tree,
    depth_column = depth_column,
    upper_bound_column = upper_bound_column,
    name_column = name_column,
    score_context = context$score_data
  )
  details$combination <- I(combinations)
  details
}

#' Decode treehunt node indexes
#'
#' @param combinations List of node-index vectors.
#' @param tree Tree data frame.
#' @param index_base Index base of `combinations`.
#' @param code_column Column containing node labels.
#' @export
decode_treehunt_solutions <- function(combinations, tree, index_base = 0L,
                                      code_column = "Code") {
  if (!code_column %in% names(tree)) {
    stop("Unknown code column.", call. = FALSE)
  }
  offset <- if (identical(as.integer(index_base), 0L)) 1L else 0L
  lapply(combinations, function(nodes) {
    rows <- as.integer(nodes) + offset
    if (any(rows < 1L | rows > nrow(tree))) {
      stop("A solution index is outside the tree.", call. = FALSE)
    }
    as.character(tree[[code_column]][rows])
  })
}

#' Refit selected combinations with patient-clustered Cox models
#'
#' This step is intentionally separate from high-throughput screening. It
#' estimates hazard ratios and robust confidence intervals for a small set of
#' selected combinations.
#'
#' @inheritParams score_pwp_combinations
#' @param keep_fits Retain fitted `coxph` objects in a list-column.
#' @export
refit_pwp_combinations <- function(
    combinations,
    data,
    context,
    tree,
    index_base = 0L,
    node_column = "icd10_indices",
    depth_column = "Depth",
    upper_bound_column = "UpperBound",
    name_column = "Code",
    keep_fits = FALSE) {
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("The 'survival' package is required.", call. = FALSE)
  }
  analysis_data <- .pwp_analysis_data(data, context)
  tree <- .pwp_tree(tree, depth_column, upper_bound_column, name_column)
  spec <- context$specification
  model_columns <- unique(c(
    spec$time_column, spec$status_column, spec$strata_column,
    spec$id_column, spec$covariates
  ))
  model_data <- as.data.frame(lapply(model_columns, function(name) {
    analysis_data[[name]]
  }), optional = TRUE)
  names(model_data) <- model_columns

  rhs <- c(
    spec$covariates, ".combo_flag",
    sprintf("strata(%s)", spec$strata_column),
    sprintf("cluster(%s)", spec$id_column)
  )
  formula <- stats::as.formula(sprintf(
    "Surv(%s, %s) ~ %s", spec$time_column, spec$status_column,
    paste(rhs, collapse = " + ")
  ))
  environment(formula) <- asNamespace("survival")

  rows <- vector("list", length(combinations))
  fits <- vector("list", length(combinations))
  for (i in seq_along(combinations)) {
    model_data$.combo_flag <- as.integer(pwp_combination_flag_cpp(
      combination = as.integer(combinations[[i]]),
      index_base = as.integer(index_base),
      patient_data = analysis_data,
      node_column = node_column,
      target_column = spec$status_column,
      id_column = spec$id_column,
      tree = tree,
      depth_column = depth_column,
      upper_bound_column = upper_bound_column,
      name_column = name_column
    ))
    fit <- survival::coxph(
      formula, data = model_data, ties = spec$ties,
      na.action = stats::na.fail
    )
    coefficient_index <- match(".combo_flag", names(stats::coef(fit)))
    if (is.na(coefficient_index)) {
      stop("The candidate indicator was dropped from the fitted model.",
           call. = FALSE)
    }
    coefficient <- stats::coef(fit)[[coefficient_index]]
    # coxph does not consistently retain dimension names on the robust
    # variance matrix.  Its rows and columns nevertheless follow coef(fit).
    variance <- fit$var[coefficient_index, coefficient_index]
    standard_error <- sqrt(variance)
    z <- coefficient / standard_error
    rows[[i]] <- data.frame(
      estimate = coefficient,
      robust_se = standard_error,
      hazard_ratio = exp(coefficient),
      conf_low = exp(coefficient - 1.96 * standard_error),
      conf_high = exp(coefficient + 1.96 * standard_error),
      robust_z = z,
      robust_p_value = 2 * stats::pnorm(abs(z), lower.tail = FALSE),
      exposed_intervals = sum(model_data$.combo_flag),
      exposed_events = sum(model_data$.combo_flag * model_data[[spec$status_column]])
    )
    if (isTRUE(keep_fits)) {
      fits[[i]] <- fit
    }
  }
  out <- do.call(rbind, rows)
  out$combination <- I(combinations)
  if (isTRUE(keep_fits)) {
    out$fit <- I(fits)
  }
  out
}
