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
