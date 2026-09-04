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
