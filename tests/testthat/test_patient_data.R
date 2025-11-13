test_that("PatientData creation with list format works", {
  tree_df <- data.frame(depth = c(1, 2, 3, 2, 3))
  tree <- make_test_tree_ctr2(tree_df$depth)

  # List format
  patient_df <- data.frame(
    patient_id = 1:3,
    outcome = c(1, 0, 1)
  )
  patient_df$drug_codes <- list(
    c(1, 2),    # Patient 1: nodes 1 and 2
    c(3),       # Patient 2: node 3
    c(2, 4)     # Patient 3: nodes 2 and 4
  )
  
  data <- .test_create_patient_data_int(patient_df, "drug_codes", "outcome", tree)
  
  expect_equal(.test_patient_data_size(data), 3)
  expect_equal(.test_patient_data_get_nodes(data, 0), c(1, 2))
  expect_equal(.test_patient_data_get_nodes(data, 1), 3)
  expect_equal(.test_patient_data_get_nodes(data, 2), c(2, 4))
  expect_equal(.test_patient_data_get_target(data, 0), 1)
  expect_equal(.test_patient_data_get_target(data, 1), 0)
})

test_that("PatientData creation with string format (comma-separated) works", {
  tree_df <- data.frame(depth = c(1, 2, 3, 2, 3))
  tree <- make_test_tree_ctr2(tree_df$depth)
  
  patient_df <- data.frame(
    patient_id = 1:3,
    drug_codes = c("1,2", "3", "2,4"),
    outcome = c(1, 0, 1),
    stringsAsFactors = FALSE
  )
  
  data <- .test_create_patient_data_int(patient_df, "drug_codes", "outcome", tree)
  
  expect_equal(.test_patient_data_size(data), 3)
  expect_equal(.test_patient_data_get_nodes(data, 0), c(1, 2))
  expect_equal(.test_patient_data_get_nodes(data, 1), 3)
  expect_equal(.test_patient_data_get_nodes(data, 2), c(2, 4))
})

test_that("PatientData handles strings with spaces", {
  tree_df <- data.frame(depth = c(1, 2, 2))
  tree <- make_test_tree_ctr2(tree_df$depth)
  
  patient_df <- data.frame(
    drug_codes = c("1, 2", " 1 ", " 2,1"),
    outcome = c(1, 0, 1),
    stringsAsFactors = FALSE
  )
  
  data <- .test_create_patient_data_int(patient_df, "drug_codes", "outcome", tree)
  
  expect_equal(.test_patient_data_get_nodes(data, 0), c(1, 2))
  expect_equal(.test_patient_data_get_nodes(data, 1), 1)
  expect_equal(.test_patient_data_get_nodes(data, 2), c(2, 1))
})

test_that("PatientData rejects plain integer vector", {
  tree_df <- data.frame(depth = c(1, 2, 2))
  tree <- make_test_tree_ctr2(tree_df$depth)
  
  patient_df <- data.frame(
    drug = c(1, 2, 3),
    outcome = c(1, 0, 1)
  )
  
  expect_error(
    .test_create_patient_data_int(patient_df, "drug", "outcome", tree),
    "Node column must be one of"
  )
})

test_that("PatientData supports double targets", {
  tree_df <- data.frame(depth = c(1, 2, 2))
  tree <- make_test_tree_ctr2(tree_df$depth)
  
  patient_df <- data.frame(
    outcome = c(0.5, 0.8, 0.3) 
  )
  patient_df$drugs <- list(c(1), c(2), c(1))
  
  data <- .test_create_patient_data_double(patient_df, "drugs", "outcome", tree)
  
  expect_equal(.test_patient_data_get_target_double(data, 0), 0.5)
  expect_equal(.test_patient_data_get_target_double(data, 1), 0.8)
  expect_equal(.test_patient_data_get_target_double(data, 2), 0.3)
})

test_that("PatientData hierarchical matching works", {
  # Tree: 0 (root) -> [1, 3], 1 -> [2], 3 -> [4]
  tree_df <- data.frame(depth = c(1, 2, 3, 2, 3))
  tree <- make_test_tree_ctr2(tree_df$depth)
  
  patient_df <- data.frame(patient = 1:3, outcome = c(1, 1, 0))
  patient_df$nodes <- list(
    c(2, 4),  # Patient 0: leaf nodes
    c(2),     # Patient 1: one leaf
    c(0)      # Patient 2: root
  )
  
  data <- .test_create_patient_data_int(patient_df, "nodes", "outcome", tree)
  
  # Solution [1, 3] covers [2, 4]
  expect_true(.test_patient_data_has_combination(data, 0, c(1, 3)))
  
  # Solution [1] only covers [2], not [4]
  expect_false(.test_patient_data_has_combination(data, 0, c(1)))
  
  # Solution [0] (root) covers everything
  expect_true(.test_patient_data_has_combination(data, 0, c(0)))
  expect_true(.test_patient_data_has_combination(data, 1, c(0)))
  expect_true(.test_patient_data_has_combination(data, 2, c(0)))
  
  # Patient 1 has [2], solution [3] doesn't cover it
  expect_false(.test_patient_data_has_combination(data, 1, c(3)))
  
  # Patient 1 has [2], solution [1] covers it
  expect_true(.test_patient_data_has_combination(data, 1, c(1)))
})

test_that("PatientData handles edge cases", {
  tree_df <- data.frame(depth = c(0, 1, 1))
  tree <- make_test_tree_ctr2(tree_df$depth)
  
  # Empty node list for a patient
  patient_df <- data.frame(outcome = c(1, 0))
  patient_df$nodes <- list(c(1), integer(0))  # Second patient has no nodes
  
  data <- .test_create_patient_data_int(patient_df, "nodes", "outcome", tree)
  
  # Patient with no nodes matches any solution (vacuous truth)
  expect_true(.test_patient_data_has_combination(data, 1, c(1)))
  expect_true(.test_patient_data_has_combination(data, 1, integer(0)))
  
  # Patient with nodes doesn't match empty solution
  expect_false(.test_patient_data_has_combination(data, 0, integer(0)))
})

test_that("PatientData column specification works by name and index", {
  tree_df <- data.frame(depth = c(0, 1, 1))
  tree <- make_test_tree_ctr2(tree_df$depth)
  
  patient_df <- data.frame(
    id = 1:3,
    outcome = c(1, 0, 1)
  )
  patient_df$drugs <- list(c(1), c(2), c(1))
  
  # By name
  data1 <- .test_create_patient_data_int(patient_df, "drugs", "outcome", tree)
  expect_equal(.test_patient_data_size(data1), 3)
  
  # By index (drugs is column 3, outcome is 2)
  data2 <- .test_create_patient_data_int(patient_df, 3L, 2L, tree)
  expect_equal(.test_patient_data_size(data2), 3)
  expect_equal(.test_patient_data_get_nodes(data2, 0), 1)
})

test_that("PatientData validates input properly", {
  tree_df <- data.frame(depth = c(0, 1, 1))
  tree <- make_test_tree_ctr2(tree_df$depth)
  
  # Empty DataFrame
  empty_df <- data.frame()
  expect_error(
    .test_create_patient_data_int(empty_df, "nodes", "outcome", tree),
    "empty"
  )
  
  # Invalid column name
  patient_df <- data.frame(outcome = c(1, 0))
  patient_df$nodes <- list(c(1), c(2))
  
  expect_error(
    .test_create_patient_data_int(patient_df, "nonexistent", "outcome", tree)
  )
})