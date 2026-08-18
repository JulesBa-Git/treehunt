test_that("optimized PWP Rao statistic matches coxph score test", {
  skip_if_not_installed("survival")
  set.seed(42)
  n_patients <- 300L
  dat <- data.frame(
    id = rep(seq_len(n_patients), each = 2L),
    event_stratum_25 = rep(1:2, n_patients),
    gap_time = sample(1:20, 2L * n_patients, replace = TRUE),
    status_any = rbinom(2L * n_patients, 1, 0.65),
    age_start = rnorm(2L * n_patients, 70, 8),
    index_emergency = rbinom(2L * n_patients, 1, 0.25),
    year_start = sample(2018:2022, 2L * n_patients, replace = TRUE)
  )
  dat$icd10_indices <- lapply(seq_len(nrow(dat)), function(i) {
    nodes <- integer()
    if (runif(1) < stats::plogis(-1 + dat$index_emergency[i])) {
      nodes <- c(nodes, 2L)
    }
    if (runif(1) < 0.25) {
      nodes <- c(nodes, 4L)
    }
    nodes
  })
  tree <- data.frame(
    Code = c("ROOT_A", "GROUP_A", "A", "ROOT_B", "B"),
    Depth = c(1L, 2L, 3L, 1L, 2L),
    UpperBound = c(2L, 2L, 2L, 4L, 4L)
  )

  context <- fit_pwp_rao_context(dat, "status_any")
  optimized <- score_pwp_combinations(
    list(2L), dat, context, tree, index_base = 0L
  )

  analysis <- dat[context$row_index, , drop = FALSE]
  analysis$.combo_flag <- as.integer(vapply(
    analysis$icd10_indices, function(x) 2L %in% x, logical(1)
  ))
  reference <- survival::coxph(
    survival::Surv(gap_time, status_any) ~
      age_start + index_emergency + year_start + .combo_flag +
      survival::strata(event_stratum_25) + survival::cluster(id),
    data = analysis,
    ties = "efron",
    init = c(stats::coef(context$fit), 0),
    control = survival::coxph.control(iter.max = 0)
  )

  expect_equal(optimized$rao_chisq, unname(reference$score), tolerance = 1e-7)
  expect_equal(optimized$covered_intervals, sum(analysis$.combo_flag))
  expect_equal(
    optimized$covered_patients,
    length(unique(analysis$id[analysis$.combo_flag == 1L]))
  )

  refit <- refit_pwp_combinations(
    list(2L), dat, context, tree, index_base = 0L
  )
  expect_true(is.finite(refit$hazard_ratio))
  expect_true(is.finite(refit$robust_se) && refit$robust_se > 0)

  strict_context <- context
  strict_context$score_data$min_covered_patients <- n_patients + 1L
  strict_context$specification$min_covered_patients <- n_patients + 1L
  filtered <- score_pwp_combinations(
    list(2L), dat, strict_context, tree, index_base = 0L
  )
  expect_false(filtered$eligible)
  expect_equal(filtered$fitness, 0)
  expect_equal(filtered$signed_z, optimized$signed_z, tolerance = 0)
})

test_that("PWP genetic algorithm is reproducible with a fixed C++ seed", {
  skip_if_not_installed("survival")
  set.seed(7)
  dat <- data.frame(
    id = rep(1:80, each = 2),
    event_stratum_25 = rep(1:2, 80),
    gap_time = sample(1:10, 160, replace = TRUE),
    status_any = rbinom(160, 1, 0.6),
    age_start = rnorm(160, 70, 5),
    index_emergency = rbinom(160, 1, 0.2),
    year_start = sample(2019:2022, 160, replace = TRUE)
  )
  dat$icd10_indices <- lapply(seq_len(nrow(dat)), function(i) {
    sample(c(2L, 4L), sample(1:2, 1))
  })
  tree <- data.frame(
    Code = c("ROOT_A", "GROUP_A", "A", "ROOT_B", "B"),
    Depth = c(1L, 2L, 3L, 1L, 2L),
    UpperBound = c(2L, 2L, 2L, 4L, 4L)
  )
  context <- fit_pwp_rao_context(dat, "status_any")

  run_once <- function() run_pwp_genetic_algorithm(
    dat, context, tree, population_size = 10L, epochs = 5L,
    elite_count = 1L, seed = 123L
  )
  first <- suppressWarnings(run_once())
  second <- suppressWarnings(run_once())
  expect_identical(first$final_population, second$final_population)
  expect_equal(first$final_scores, second$final_scores, tolerance = 0)
})
