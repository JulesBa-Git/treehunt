test_that("clustering input prefers the exact vec_cocktails column", {
  cocktails <- data.frame(label = c("first", "second"))
  cocktails$vec_cocktails <- list(c(1L, 2L), 3L)
  cocktails$vec_cocktail <- list(9L, 10L)

  expect_equal(.extract_cocktail_list(cocktails), cocktails$vec_cocktails)
})

test_that("clustering input requires one unambiguous list-column", {
  no_list <- data.frame(x = 1:2)
  expect_error(.extract_cocktail_list(no_list), "exactly one list-column")

  ambiguous <- data.frame(x = 1:2)
  ambiguous$a <- list(1L, 2L)
  ambiguous$b <- list(3L, 4L)
  expect_error(.extract_cocktail_list(ambiguous), "exactly one list-column")
})

test_that("dissimilarities handle roots, separate trees and empty populations", {
  dat <- data.frame(outcome = c(1L, 0L))
  dat$nodes <- list(1L, 2L)
  tree <- data.frame(Depth = c(1L, 2L, 1L))
  distance <- get_dissimilarity_of_list(list(1L, 2L, 3L), dat, "nodes", "outcome", tree, "Depth")
  expect_equal(dim(distance), c(3L, 3L))
  expect_true(all(is.finite(distance)))
  expect_equal(distance, t(distance))
  expect_equal(diag(distance), rep(0, 3))
  expect_true(all(distance >= 0))
  expect_equal(dim(get_dissimilarity_of_list(list(), dat, "nodes", "outcome", tree, "Depth")), c(0L, 0L))
  ga <- run_genetic_algorithm(dat, "nodes", "outcome", tree$Depth,
                             population_size = 4, epochs = 1, tournament_size = 2,
                             diversity = TRUE, seed = 14)
  expect_true(all(is.finite(ga$final_scores)))
})
