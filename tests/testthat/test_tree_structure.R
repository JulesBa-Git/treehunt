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

test_that("tree creation with ATC_tree match",{
  library(emcAdr)
  modified_ATC <- ATC_Tree_UpperBound_2024
  modified_ATC$ATC_length[modified_ATC$ATC_length == 3] <- 2
  modified_ATC$ATC_length[modified_ATC$ATC_length == 4] <- 3
  modified_ATC$ATC_length[modified_ATC$ATC_length == 5] <- 4
  modified_ATC$ATC_length[modified_ATC$ATC_length == 7] <- 5
  depth_vec <- modified_ATC$ATC_length
  df <- data.frame(
    id = 1:length(modified_ATC$ATC_length),
    depth = depth_vec
  )
  
  tree <- make_test_tree_ctr1(df, "depth")
  info <- inspect_tree(tree)
  
  expect_equal(info$max_depth, 5)
  expect_false(info$has_name)
  expect_equal(info$upper_bound, modified_ATC$upperBound - 1)
  
  tree_vec <- make_test_tree_ctr2(depth_vec)
  info_vec <- inspect_tree(tree_vec)
  
  
  expect_equal(info$max_depth, info_vec$max_depth)
  expect_equal(info$upper_bound, info_vec$upper_bound)
  
  tree_2 <- make_test_tree_ctr1(df, 2)
  info_2 <- inspect_tree(tree_2)
  
  expect_equal(info$max_depth, info_2$max_depth)
  expect_equal(info$upper_bound, info_2$upper_bound)
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