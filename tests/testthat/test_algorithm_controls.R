make_algorithm_fixture <- function() {
  tree <- data.frame(
    Depth = rep(1L, 4L),
    UpperBound = 0:3,
    Name = LETTERS[1:4]
  )
  observations <- data.frame(outcome = c(1L, 0L, 1L, 0L))
  observations$nodes <- list(0L, 1L, 2L, 3L)
  list(tree = tree, observations = observations)
}

test_that("public GA uses prob_mutation_type1", {
  fixture <- make_algorithm_fixture()
  common <- list(
    patient_data = fixture$observations,
    node_column = "nodes",
    target_column = "outcome",
    tree = fixture$tree,
    depth_column = "Depth",
    upper_bound_column = "UpperBound",
    seed_population = rep(list(1L), 4L),
    population_size = 4L,
    epochs = 1L,
    mutation_rate = 1,
    crossover_rate = 0,
    elite_count = 0L,
    tournament_size = 2L,
    alpha = 1,
    seed = 41L
  )

  only_type1 <- do.call(
    run_genetic_algorithm_df_tree,
    c(common, list(prob_mutation_type1 = 1))
  )
  only_type2 <- do.call(
    run_genetic_algorithm_df_tree,
    c(common, list(prob_mutation_type1 = 0))
  )

  expect_true(all(lengths(only_type1$final_population) == 2L))
  expect_true(all(lengths(only_type2$final_population) == 1L))
})

test_that("generic GA and MCMC seeds reproduce their C++ streams", {
  fixture <- make_algorithm_fixture()

  ga_df_args <- list(
    patient_data = fixture$observations,
    node_column = "nodes",
    target_column = "outcome",
    tree = fixture$tree,
    depth_column = "Depth",
    upper_bound_column = "UpperBound",
    population_size = 6L,
    epochs = 4L,
    elite_count = 1L,
    tournament_size = 2L,
    seed = 73L
  )
  expect_identical(
    do.call(run_genetic_algorithm_df_tree, ga_df_args),
    do.call(run_genetic_algorithm_df_tree, ga_df_args)
  )

  ga_vector_args <- list(
    patient_data = fixture$observations,
    node_column = "nodes",
    target_column = "outcome",
    tree_depth = fixture$tree$Depth,
    population_size = 6L,
    epochs = 4L,
    elite_count = 1L,
    tournament_size = 2L,
    seed = 74L
  )
  expect_identical(
    do.call(run_genetic_algorithm, ga_vector_args),
    do.call(run_genetic_algorithm, ga_vector_args)
  )

  mcmc_vector_args <- list(
    patient_data = fixture$observations,
    node_column = "nodes",
    target_column = "outcome",
    tree_depth = fixture$tree$Depth,
    epochs = 30L,
    cocktail_size = 1L,
    max_score = 20,
    seed = 91L
  )
  expect_identical(
    do.call(run_mcmc, mcmc_vector_args),
    do.call(run_mcmc, mcmc_vector_args)
  )

  mcmc_df_args <- list(
    patient_data = fixture$observations,
    node_column = "nodes",
    target_column = "outcome",
    tree = fixture$tree,
    depth_column = "Depth",
    upper_bound_column = "UpperBound",
    epochs = 30L,
    cocktail_size = 1L,
    max_score = 20,
    seed = 92L
  )
  expect_identical(
    do.call(run_mcmc_df_tree, mcmc_df_args),
    do.call(run_mcmc_df_tree, mcmc_df_args)
  )
})

test_that("generic seeds and GA mutation probabilities are validated", {
  fixture <- make_algorithm_fixture()
  base_args <- list(
    patient_data = fixture$observations,
    node_column = "nodes",
    target_column = "outcome",
    tree = fixture$tree,
    depth_column = "Depth",
    population_size = 4L,
    epochs = 1L
  )

  for (bad_seed in list(-1L, NA_integer_, 1.5, c(1L, 2L))) {
    expect_error(
      do.call(run_genetic_algorithm_df_tree,
              c(base_args, list(seed = bad_seed))),
      "seed must be NULL or one non-negative integer"
    )
  }

  for (bad_probability in list(-0.1, 1.1, NA_real_, NaN, Inf)) {
    expect_error(
      do.call(run_genetic_algorithm_df_tree,
              c(base_args, list(prob_mutation_type1 = bad_probability))),
      "prob_mutation_type1 must be finite and between 0 and 1"
    )
  }
})
