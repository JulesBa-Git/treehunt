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


test_that("hypergeometric scores agree with Fisher tails including zero events", {
  fixture <- make_compute_score_fixture()
  result <- compute_score(list(1L, 2L, 3L), fixture$observations,
                          "nodes", "outcome", fixture$tree, "Depth",
                          score_type = "hypergeometric")
  expect_equal(result$scores, c(0, log(2), 0))
  expect_true(all(is.finite(result$scores)))
})

test_that("patient-level scores use medians and are invariant to row order", {
  medians <- c(11:14, -1:-4)
  dat <- data.frame(id = rep(1:8, each = 3),
                    outcome = unlist(lapply(medians, function(m) c(1000, -1000, m))))
  dat$nodes <- c(rep(list(0L), 12), rep(list(1L), 12))
  tree <- data.frame(Depth = c(1L, 1L))
  score <- function(data, type) compute_score(list(1L), data, "nodes", "outcome", tree,
                                              "Depth", id_column = "id", score_type = type)
  result <- score(dat, "wilcoxon")
  expect_equal(sort(result$QT_diff_distribution[[1]]), as.numeric(11:14))
  reordered <- dat[c(seq(3, 24, 3), seq(1, 24, 3), seq(2, 24, 3)), ]
  for (type in c("wilcoxon", "residuals")) {
    expect_equal(score(dat, type)$scores, score(reordered, type)$scores, tolerance = 1e-12)
  }
})
