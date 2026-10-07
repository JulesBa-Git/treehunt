test_that("process_ga_scores binds compute_score columns by name", {
  tree <- data.frame(
    Depth = c(1L, 1L),
    UpperBound = c(0L, 1L),
    Name = c("A", "B")
  )
  observations <- data.frame(
    patient_id = c(1L, 1L, 2L, 2L, 3L, 3L),
    delta = c(1, 2, 3, 4, 5, 6)
  )
  observations$nodes <- list(0L, 0L, 1L, 1L, 0L, 1L)
  candidates <- data.frame(cocktail = "0")

  result <- process_ga_scores(
    df = candidates,
    patient_data = observations,
    tree_df = tree,
    node_column = "nodes",
    target_column = "delta",
    depth_column = "Depth",
    upper_bound_column = "UpperBound",
    score_type = "wilcoxon",
    id_column = "patient_id",
    name_column = "Name"
  )

  expect_equal(result$taker_count, 2)
  expect_true(is.finite(result$scores))
})

test_that("saved GA results aggregate without attached tidyverse packages", {
  folder <- tempfile("treehunt-results-")
  dir.create(folder)
  on.exit(unlink(folder, recursive = TRUE))
  runs <- list(
    list(final_population = list(c(0L, 2L), 1L), final_scores = c(3, 2),
         metadata = list(config_name = "A")),
    list(final_population = list(c(2L, 0L)), final_scores = 3,
         metadata = list(config_name = "B"))
  )
  jsonlite::write_json(runs, file.path(folder, "results.json"))
  result <- aggregate_ga_results(folder)
  expect_equal(result$cocktail, c("0,2", "1"))
  expect_equal(result$score, c(3, 2))
  expect_equal(result$occurrence_count, c(2L, 1L))
  expect_equal(result$found_in_configs, c("A; B", "A"))
  tree <- data.frame(Name = c("Root", "A", "B"), Code = c("R", "A1", "B1"))
  named <- map_cocktail_names(result, tree)
  expect_equal(named$cocktail_names, c("Root | B", "A"))
  expect_equal(named$cocktail_codes, c("R | B1", "A1"))
  expect_equal(filter_out_cocktails(named, c(1, 2, 3), 2.1)$cocktail, character())
  expect_equal(filter_out_cocktails(named, c(1, 2, 3), 2)$cocktail, result$cocktail)
  expect_equal(filter_out_cocktails(c("1,3", "2"), c(1, 2, 3), 2, one_index = TRUE)$cocktail,
               c("1,3", "2"))
})

test_that("GA batches require an explicit directory without writing by default", {
  work <- tempfile("treehunt-batch-work-")
  dir.create(work)
  old_wd <- setwd(work)
  on.exit({
    setwd(old_wd)
    unlink(work, recursive = TRUE)
  })

  expect_error(run_ga_batch("missing.json", NULL, NULL),
               "Supply 'output_dir' explicitly")
  for (path in list(NULL, character(), "", "  ", NA_character_, c("a", "b"), 1)) {
    expect_error(run_ga_batch("missing.json", NULL, NULL, output_dir = path),
                 "Supply 'output_dir' explicitly")
  }
  expect_length(list.files(work, all.files = TRUE, no.. = TRUE), 0L)
})

test_that("GA batches write readable replicates to the requested directory", {
  work <- tempfile("treehunt-batch-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE))
  config_path <- file.path(work, "config.json")
  output_dir <- file.path(work, "explicit-output")
  config <- data.frame(
    name = "toy", node_column = "nodes", target_column = "outcome",
    depth_column = "Depth", population_size = 4L, epochs = 2L,
    mutation_rate = 0.1, prob_mutation_type1 = 0.2, crossover_rate = 0.8,
    elite_count = 0L, tournament_size = 2L, alpha = 1,
    score_type = "hypergeometric", diversity = FALSE, verbose = FALSE
  )
  jsonlite::write_json(config, config_path)
  observations <- data.frame(outcome = c(1L, 0L, 1L, 0L))
  observations$nodes <- list(1L, 2L, c(1L, 2L), 2L)
  tree <- data.frame(Depth = c(1L, 2L, 2L))

  result <- withVisible(suppressMessages(run_ga_batch(
    config_path, observations, tree,
    seed_population = list(2L, 3L, c(2L, 3L), 2L),
    replicates = 2L, output_dir = output_dir
  )))
  expect_null(result$value)
  expect_false(result$visible)
  expect_equal(list.files(output_dir), "results_toy.json")
  saved <- jsonlite::fromJSON(file.path(output_dir, "results_toy.json"),
                              simplifyVector = FALSE)
  expect_length(saved, 2L)
  expect_equal(vapply(saved, function(x) {
    as.integer(unlist(x$metadata$replicate))
  }, integer(1)), 1:2)
  expect_true(all(vapply(saved, function(x) {
    length(x$final_population) == 4L && length(x$final_scores) == 4L &&
      identical(unlist(x$metadata$config_name, use.names = FALSE), "toy")
  }, logical(1))))
  aggregated <- aggregate_ga_results(output_dir)
  expect_equal(sum(aggregated$occurrence_count), 8L)
  expect_true(all(is.finite(aggregated$score)))
})
