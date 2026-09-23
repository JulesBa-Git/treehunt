test_that("tree creation works",{
  depth_vec <- c(1,2,3,4,5,5,5,5,5,5)
  df <- data.frame(
    id = 1:10,
    depth = depth_vec,
    name = rep("bla", 10)
  )
  
  expected_upper_bound <- c(9,9,9,9,4,5,6,7,8,9)
  
  tree <- make_test_tree_ctr1(df, "depth")
  info <- inspect_tree(tree)
  
  expect_equal(info$max_depth, 5)
  expect_false(info$has_name)
  expect_equal(info$upper_bound, expected_upper_bound)
  
  tree_vec <- make_test_tree_ctr2(depth_vec)
  info_vec <- inspect_tree(tree_vec)

  
  expect_equal(info$max_depth, info_vec$max_depth)
  expect_equal(info$upper_bound, info_vec$upper_bound)
  
  tree <- make_test_tree_ctr1(df, "depth", name = "name")
  info <- inspect_tree(tree)
  
  expect_true(info$has_name)
  expect_equal(info$name, rep("bla", 10))
  expect_equal(info$upper_bound, expected_upper_bound)
})

test_that("tree constructors agree with independently specified forest bounds", {
  depth <- c(1L, 2L, 3L, 2L, 1L, 2L, 3L, 3L)
  bounds <- c(3L, 2L, 2L, 3L, 7L, 7L, 6L, 7L)
  df <- data.frame(id = seq_along(depth), depth = depth)
  by_name <- inspect_tree(make_test_tree_ctr1(df, "depth"))
  by_index <- inspect_tree(make_test_tree_ctr1(df, 2))
  by_vector <- inspect_tree(make_test_tree_ctr2(depth))
  expect_equal(by_name$upper_bound, bounds)
  expect_equal(by_index$upper_bound, bounds)
  expect_equal(by_vector$upper_bound, bounds)
  expect_equal(by_name$max_depth, 3)
})

test_that("cannot create a tree from wrong depth or wrong index",{
  bad_depth <- c(1,2,3,5,1,3,5,5)
  df <- data.frame(
    id = 1:length(bad_depth),
    depth = bad_depth
  )
  
  expect_error(
    make_test_tree_ctr2(bad_depth),
    "son of a node must be in the next depth"
  )
  
  expect_error(
    make_test_tree_ctr1(df, "depth"),
    "son of a node must be in the next depth"
  )
  
  expect_error(
    make_test_tree_ctr1(df, "dep")
  )
  
})

test_that("empty tree should stop program",{
  depth_vec <- integer(0)
  
  expect_error(make_test_tree_ctr2(depth_vec),
               "No value in depth vector")
})
test_that("zero, negative and missing depths cannot index parent arrays", {
  for (depth in list(0L, -1L, c(1L, 0L), NA_integer_)) {
    expect_error(make_test_tree_ctr2(depth), "positive integers")
    expect_error(make_test_tree_ctr1(data.frame(depth = depth), "depth"), "positive integers")
  }
})
