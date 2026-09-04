make_compute_score_fixture <- function() {
  tree <- data.frame(
    Depth = c(1L, 2L, 2L),
    UpperBound = c(2L, 1L, 2L),
    Name = c("root", "A", "B")
  )
  observations <- data.frame(outcome = c(1L, 0L))
  observations$nodes <- list(1L, 2L)
  list(tree = tree, observations = observations)
}

test_that("compute_score keeps its one-based candidate interface", {
  fixture <- make_compute_score_fixture()

  result <- compute_score(
    cocktail_list = list(1L),
    patient_data = fixture$observations,
    node_column = "nodes",
    target_column = "outcome",
    tree = fixture$tree,
    depth_column = "Depth",
    upper_bound_column = "UpperBound",
    name_column = "Name",
    score_type = "hypergeometric"
  )

  expect_equal(result$solutions, list(1L))
  expect_equal(result$`number of takers`, 2)
})

test_that("compute_score rejects invalid one-based tree rows before conversion", {
  fixture <- make_compute_score_fixture()
  score <- function(cocktails) {
    compute_score(
      cocktail_list = cocktails,
      patient_data = fixture$observations,
      node_column = "nodes",
      target_column = "outcome",
      tree = fixture$tree,
      depth_column = "Depth",
      upper_bound_column = "UpperBound"
    )
  }

  for (bad in list(0L, -1L, 4L, NA_integer_, c(1L, 4L))) {
    expect_error(score(list(bad)), "indices must lie between 1 and 3")
  }
})

test_that("continuous compute_score calls require an identifier", {
  fixture <- make_compute_score_fixture()
  fixture$observations$outcome <- c(1.5, 2.5)

  expect_error(
    compute_score(
      cocktail_list = list(1L),
      patient_data = fixture$observations,
      node_column = "nodes",
      target_column = "outcome",
      tree = fixture$tree,
      depth_column = "Depth"
    ),
    "id_column is required"
  )
})
