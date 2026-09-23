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
