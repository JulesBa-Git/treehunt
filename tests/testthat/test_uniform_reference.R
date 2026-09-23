test_that("weight ESS has its defining values and is scale invariant", {
  expect_equal(weight_ess(rep(1, 5)), 5)
  expect_equal(weight_ess(c(1, 0, 0)), 1)
  expect_equal(weight_ess(c(1, rep(0.25, 4))), 3.2)
  expect_equal(weight_ess(c(1, 2, 3)), weight_ess(c(1, 2, 3) * 1e250))
  expect_equal(weight_ess(c(1, 2, 3)), weight_ess(log(c(1, 2, 3)) - 10000, log = TRUE))
  expect_equal(weight_ess(c(-Inf, 0), log = TRUE), 1)
  expect_true(is.na(weight_ess(numeric())))
})

test_that("invalid importance weights are rejected", {
  for (x in list(c(0, 0), -1, NA_real_, Inf, NaN, "1", matrix(1))) {
    expect_error(weight_ess(x))
  }
  expect_error(weight_ess(c(-Inf, -Inf), log = TRUE))
  expect_error(weight_ess(Inf, log = TRUE))
  expect_error(weight_ess(1, log = NA))
})

test_that("inverse target weights recover a known uniform reference", {
  scores <- c(log(4), 0, log(4), log(4), log(4))
  ref <- uniform_score_reference(scores)
  expect_equal(ref$distribution$score, c(0, log(4)))
  expect_equal(ref$distribution$probability, c(0.5, 0.5))
  expect_equal(ref$distribution$cdf, c(0.5, 1))
  expect_equal(ref$distribution$upper_tail, c(1, 0.5))
  expect_equal(ref$weight_ess, 3.2)
  expect_equal(ref$n_samples, 5)
  expect_equal(sum(ref$weights), 1)
  # A different temperature changes the sampling weights, not the score order.
  expect_equal(uniform_score_reference(scores * 2, temperature = 2)$weights, ref$weights)
  # Equal score values are aggregated; visits are never deduplicated first.
  expect_equal(uniform_score_reference(c(0, 0, 0))$weight_ess, 3)
})

test_that("caps, extreme scores, empty input and tail ties are handled", {
  ref <- uniform_score_reference(c(1000, 1001, 1001))
  expect_equal(ref$distribution$probability, c(1, 2 / exp(1)) / (1 + 2 / exp(1)))
  expect_equal(uniform_score_reference(c(0, 1e300), temperature = 1e-300)$weight_ess, 1)
  expect_equal(uniform_score_reference(c(1e300, 1e300), temperature = 1e-300)$weight_ess, 2)
  capped <- uniform_score_reference(c(0, 2, 1000, Inf), max_score = 2)
  expect_equal(capped$weights, c(1, rep(exp(-2), 3)) / (1 + 3 * exp(-2)))
  expect_equal(capped$distribution$score, c(0, 2))
  expect_true(is.na(uniform_score_reference(numeric())$weight_ess))
  expect_equal(nrow(uniform_score_reference(numeric())$distribution), 0)
  for (t in list(0, -1, Inf, NaN, NA_real_, c(1, 2))) {
    expect_error(uniform_score_reference(1, temperature = t), "temperature")
  }
  for (x in list(NA_real_, Inf, -Inf, matrix(1), "1")) {
    expect_error(uniform_score_reference(x))
  }
  expect_error(uniform_score_reference(1, max_score = 0), "max_score")
})

reference_fixture <- function() {
  dat <- data.frame(outcome = c(1L, 1L, 0L, 0L, 0L, 0L))
  dat$nodes <- list(0L, 0L, c(0L, 1L), 1L, 2L, 2L)
  list(dat = dat, tree = data.frame(Depth = rep(1L, 3)))
}

expect_stream_matches_trace <- function(result, temperature, beta) {
  trace <- result$trace
  for (filtered in c(FALSE, TRUE)) {
    keep <- if (filtered) trace$covered_patients > beta else rep(TRUE, nrow(trace))
    scores <- trace$score[keep]
    ref <- if (filtered) result$uniform_reference_filtered else result$uniform_reference
    expect_equal(ref$n_samples, length(scores))
    if (!length(scores)) {
      expect_true(is.na(ref$weight_ess))
      expect_true(all(is.na(ref$distribution$probability)))
      next
    }
    weights <- exp(-(scores - min(scores)) / temperature)
    weights <- weights / sum(weights)
    bins <- vapply(scores, function(s) max(which(ref$distribution$lower <= s)), integer(1))
    expected <- vapply(seq_len(nrow(ref$distribution)), function(j) sum(weights[bins == j]), numeric(1))
    expect_equal(ref$distribution$probability, expected, tolerance = 1e-12)
    expect_equal(ref$weight_ess, 1 / sum(weights^2), tolerance = 1e-10)
    expect_equal(ref$distribution$upper_tail, rev(cumsum(rev(expected))), tolerance = 1e-12)
    expect_equal(ref$distribution$cdf, cumsum(expected), tolerance = 1e-12)
  }
}

test_that("streaming weights agree with exact post-processing and retain rejects", {
  fixture <- reference_fixture()
  result <- run_mcmc(fixture$dat, "nodes", "outcome", fixture$tree$Depth,
                    epochs = 2000, burn_in = 100, cocktail_size = 1, beta = 2,
                    prob_type1 = 0.4, temperature = 0.8, max_score = 1.25,
                    seed = 207, store_trace = TRUE)
  expect_equal(nrow(result$trace), 1900)
  expect_equal(sum(result$score_distribution), 1900)
  expect_gt(result$statistics$rejected_moves, 0)
  expect_equal(result$statistics$accepted_moves + result$statistics$rejected_moves, 2000)
  expect_equal(sum(result$score_distribution_filtered), sum(result$trace$covered_patients > 2))
  expect_stream_matches_trace(result, 0.8, 2)
  expect_true(all(result$trace$score <= 1.25))
  expect_true(all(diff(result$uniform_reference$distribution$lower) > 0))
  expect_equal(tail(result$uniform_reference$distribution$lower, 1), 1.25)
  without_trace <- run_mcmc(fixture$dat, "nodes", "outcome", fixture$tree$Depth,
                           epochs = 10, cocktail_size = 1, seed = 1, n_results = 0,
                           beta = 100)
  expect_false("trace" %in% names(without_trace))
  expect_length(without_trace$top_scores, 0)
  expect_true(is.na(without_trace$uniform_reference_filtered$weight_ess))
})

test_that("burn-in discards a prefix without changing the trajectory", {
  fixture <- reference_fixture()
  args <- list(patient_data = fixture$dat, node_column = "nodes", target_column = "outcome",
               tree = fixture$tree, depth_column = "Depth", epochs = 200,
               cocktail_size = 1, seed = 32, store_trace = TRUE)
  all <- do.call(run_mcmc_df_tree, args)
  retained <- do.call(run_mcmc_df_tree, c(args, list(burn_in = 50)))
  expect_equal(retained$trace$score, all$trace$score[51:200])
  expect_equal(retained$trace$solution, all$trace$solution[51:200])
  expect_stream_matches_trace(retained, 1, 4)
})

test_that("both proposal types preserve fixed-size distinct-node subsets", {
  # All six pairs are observed, including ancestor-descendant pairs.
  dat <- data.frame(outcome = rep(0L, 4))
  dat$nodes <- rep(list(c(1L, 2L, 3L)), 4)
  tree <- data.frame(Depth = c(1L, 2L, 2L, 1L))
  result <- run_mcmc_df_tree(dat, "nodes", "outcome", tree, "Depth",
                            cocktail_size = 2, epochs = 30000, burn_in = 100,
                            prob_type1 = 0.2, seed = 452, store_trace = TRUE)
  combos <- result$trace$solution
  expect_true(all(lengths(combos) == 2L))
  expect_true(all(vapply(combos, function(x) length(unique(x)) == 2L, logical(1))))
  frequencies <- table(vapply(combos, paste, collapse = ",", FUN.VALUE = character(1))) / length(combos)
  expect_length(frequencies, 6)
  expect_lt(max(abs(as.numeric(frequencies) - 1/6)), 0.02)
  expect_equal(result$uniform_reference$weight_ess, 29900)
})

test_that("tempered sampling recovers exhaustive pair score probabilities", {
  dat <- data.frame(outcome = c(1L, 1L, 0L, 0L, 0L, 0L, 0L, 0L))
  dat$nodes <- list(c(0L, 1L), c(0L, 1L), c(0L, 2L), c(0L, 2L),
                    c(1L, 2L), c(1L, 2L), 0L, 1L)
  tree <- data.frame(Depth = rep(1L, 3))
  exact <- mcmc_size_2_true_score_distribution(dat, "nodes", "outcome", tree, "Depth", beta = 1)
  expect_equal(sum(exact$score_distribution), 3)
  expect_equal(exact$uniform_reference$weight_ess, 3)
  expect_equal(exact$statistics$total_iterations, 3)
  expect_true(is.na(exact$statistics$acceptance_rate))
  # Independently evaluate each supported pair using Fisher's upper tail.
  scores <- -phyper(c(2, 0, 0) - 1, m = 2, n = 6, k = 2, lower.tail = FALSE, log.p = TRUE)
  expected <- tabulate(floor(scores * 10) + 1, nbins = length(exact$score_distribution)) / 3
  expect_equal(exact$uniform_reference$distribution$probability, expected)
  result <- run_mcmc_df_tree(dat, "nodes", "outcome", tree, "Depth", epochs = 40000,
                            burn_in = 1000, prob_type1 = 1, temperature = 2,
                            cocktail_size = 2, beta = 1, seed = 537)
  estimate <- result$uniform_reference$distribution$probability
  expect_lt(max(abs(estimate - expected)), 0.02)
  # Raw tempered frequencies are deliberately different from the uniform law.
  expect_gt(max(abs(result$score_distribution / sum(result$score_distribution) - expected)), 0.25)
})

test_that("enumeration excludes absent pairs and never enumerates duplicate nodes", {
  dat <- data.frame(outcome = c(1L, 0L))
  dat$nodes <- list(1L, 2L)
  tree <- data.frame(Depth = c(1L, 2L, 1L))
  exact <- mcmc_size_2_true_score_distribution(dat, "nodes", "outcome", tree, "Depth")
  # {0,1} is covered by the first row through hierarchical inheritance.
  expect_equal(sum(exact$score_distribution), 1)
  expect_equal(exact$uniform_reference$n_samples, 1)
  result <- run_mcmc_df_tree(dat, "nodes", "outcome", tree, "Depth", epochs = 30,
                            cocktail_size = 2, seed = 21, store_trace = TRUE)
  expect_true(all(vapply(result$trace$solution, identical, logical(1), c(0L, 1L))))
  expect_equal(result$uniform_reference$weight_ess, 30)
})

test_that("MCMC rejects invalid controls and unsupported sizes without hanging", {
  fixture <- reference_fixture()
  base <- list(patient_data = fixture$dat, node_column = "nodes", target_column = "outcome",
               tree_depth = fixture$tree$Depth)
  for (name in c("temperature", "max_score")) {
    for (value in list(0, -1, Inf, NA_real_)) {
      expect_error(do.call(run_mcmc, c(base, list(epochs = 10), setNames(list(value), name))))
    }
  }
  for (name in c("epochs", "n_results", "cocktail_size", "beta", "burn_in")) {
    for (value in list(-1, 1.5, NA_real_, Inf)) {
      args <- c(base, list(epochs = 10)); args[[name]] <- value
      expect_error(do.call(run_mcmc, args), name)
    }
  }
  expect_error(do.call(run_mcmc, c(base, list(epochs = 10, burn_in = 10))), "burn_in")
  expect_error(do.call(run_mcmc, c(base, list(epochs = 10, prob_type1 = 0))), "prob_type1")
  expect_error(do.call(run_mcmc, c(base, list(epochs = 10, cocktail_size = 3))), "No observed combination")
  expect_error(do.call(run_mcmc, c(base, list(epochs = 10, cocktail_size = 4))), "cocktail_size")
})

test_that("opposite extreme finite scores are scaled before overflowing", {
  ref <- uniform_score_reference(c(-1e308, 1e308), temperature = 1e308)
  expect_equal(ref$weights, c(1, exp(-2)) / (1 + exp(-2)))
})

test_that("the continuous-outcome engine returns the same weighted reference", {
  dat <- data.frame(id = 1:30, outcome = as.numeric(1:30))
  dat$nodes <- c(rep(list(0L), 15), rep(list(1L), 15))
  result <- run_mcmc(dat, "nodes", "outcome", c(1L, 1L),
                    epochs = 200, cocktail_size = 1, score_type = "wilcoxon",
                    temperature = 3, max_score = 10, seed = 821, store_trace = TRUE, id_column = "id")
  expect_stream_matches_trace(result, 3, 4)
})
