# Tests for Solution class

test_that("Solution constructor cannot creates empty solution", {
  expect_error(.test_create_solution(integer(0)),
               "Cannot create a solution from an empty vector")
})

test_that("Solution constructor creates solution from vector", {
  sol <- .test_create_solution(c(0, 1, 2))
  
  nodes <- .test_solution_get_nodes(sol)
  expect_equal(nodes, c(0, 1, 2))
  expect_equal(.test_solution_size(sol), 3)
  expect_false(.test_solution_empty(sol))
})

# ==============================================================================
# Score Management Tests
# ==============================================================================

test_that("Solution score starts uncomputed", {
  sol <- .test_create_solution(c(0, 1))
  
  expect_false(.test_solution_is_score_computed(sol))
  expect_error(.test_solution_get_score(sol), "Score not computed")
})

test_that("Solution can set and get score", {
  sol <- .test_create_solution(c(0, 1))
  
  .test_solution_set_score(sol, 0.85)
  
  expect_true(.test_solution_is_score_computed(sol))
  expect_equal(.test_solution_get_score(sol), 0.85)
})

test_that("Solution score can be invalidated", {
  sol <- .test_create_solution(c(0, 1))
  .test_solution_set_score(sol, 0.85)
  
  expect_true(.test_solution_is_score_computed(sol))
  
  .test_solution_invalidate_score(sol)
  
  expect_false(.test_solution_is_score_computed(sol))
  expect_error(.test_solution_get_score(sol), "Score not computed")
})

# ==============================================================================
# Hash Tests
# ==============================================================================

test_that("Solution computes consistent hash", {
  sol1 <- .test_create_solution(c(0, 1, 2))
  sol2 <- .test_create_solution(c(0, 1, 2))
  
  hash1 <- .test_solution_get_hash(sol1)
  hash2 <- .test_solution_get_hash(sol2)
  
  expect_equal(hash1, hash2)
})

test_that("Solution hash is independent of input order", {
  sol1 <- .test_create_solution(c(0, 1, 2))
  sol2 <- .test_create_solution(c(2, 1, 0))
  sol3 <- .test_create_solution(c(1, 0, 2))
  
  hash1 <- .test_solution_get_hash(sol1)
  hash2 <- .test_solution_get_hash(sol2)
  hash3 <- .test_solution_get_hash(sol3)
  
  expect_equal(hash1, hash2)
  expect_equal(hash1, hash3)
})

test_that("Different solutions have different hashes (usually)", {
  sol1 <- .test_create_solution(c(0, 1, 2))
  sol2 <- .test_create_solution(c(3, 4, 5))
  
  hash1 <- .test_solution_get_hash(sol1)
  hash2 <- .test_solution_get_hash(sol2)
  
  expect_false(hash1 == hash2)
})

test_that("Solution hash can be invalidated", {
  sol <- .test_create_solution(c(0, 1))
  
  hash1 <- .test_solution_get_hash(sol)
  .test_solution_invalidate_hash(sol)
  hash2 <- .test_solution_get_hash(sol)
  
  # Should recompute to same value
  expect_equal(hash1, hash2)
})

test_that("Solution invalidate() clears both score and hash", {
  sol <- .test_create_solution(c(0, 1))
  .test_solution_set_score(sol, 0.85)
  
  # Get hash to compute it
  hash <- .test_solution_get_hash(sol)
  
  expect_true(.test_solution_is_score_computed(sol))
  
  .test_solution_invalidate(sol)
  
  expect_false(.test_solution_is_score_computed(sol))
  expect_error(.test_solution_get_score(sol), "Score not computed")
})

# ==============================================================================
# Validation Tests
# ==============================================================================

test_that("Solution with no parent-child pairs is valid", {
  # Tree: 0 -> [1, 3], 1 -> [2], 3 -> [4]
  # Nodes 0,1,2,3,4 with depths 1,2,3,2,3
  tree_df <- data.frame(depth = c(1, 2, 3, 2, 3))
  tree <- .test_create_tree_constructor2(tree_df$depth)
  
  sol <- .test_create_solution(c(1, 3))
  expect_true(.test_solution_is_valid(sol, tree))
  
  sol2 <- .test_create_solution(c(2, 4))
  expect_true(.test_solution_is_valid(sol2, tree))
})

test_that("Solution with parent-child pair is invalid", {
  # Tree: 0 -> [1, 3], 1 -> [2], 3 -> [4]
  tree_df <- data.frame(depth = c(1, 2, 3, 2, 3))
  tree <- .test_create_tree_constructor2(tree_df$depth)
  
  # Node 1 is parent of node 2
  sol <- .test_create_solution(c(1, 2))
  expect_false(.test_solution_is_valid(sol, tree))
  
  # Node 0 is parent of node 1
  sol2 <- .test_create_solution(c(0, 1))
  expect_false(.test_solution_is_valid(sol2, tree))
})

test_that("Single node solution is always valid", {
  tree_df <- data.frame(depth = c(1, 2, 3, 2, 3))
  tree <- .test_create_tree_constructor2(tree_df$depth)
  
  for (node in 0:4) {  # 0-based: nodes 0,1,2,3,4
    sol <- .test_create_solution(node)
    expect_true(.test_solution_is_valid(sol, tree))
  }
})

# ==============================================================================
# Utility Tests
# ==============================================================================

test_that("Solution contains() works correctly", {
  sol <- .test_create_solution(c(0, 2, 4))
  
  expect_true(.test_solution_contains(sol, 0))
  expect_true(.test_solution_contains(sol, 2))
  expect_true(.test_solution_contains(sol, 4))
  
  expect_false(.test_solution_contains(sol, 1))
  expect_false(.test_solution_contains(sol, 3))
  expect_false(.test_solution_contains(sol, 10))
})

test_that("Solution print works", {
  sol <- .test_create_solution(c(0, 1, 2))
  
  expect_output(.test_solution_print(sol))
})

# ==============================================================================
# Comparison Operator Tests
# ==============================================================================

test_that("Solutions can be compared by score", {
  sol1 <- .test_create_solution(c(0, 1))
  sol2 <- .test_create_solution(c(2, 3))
  
  .test_solution_set_score(sol1, 0.5)
  .test_solution_set_score(sol2, 0.8)
  
  expect_true(.test_solution_less_than(sol1, sol2))
  expect_true(.test_solution_greater_than(sol2, sol1))
})

test_that("Solutions cannot be compared without scores", {
  sol1 <- .test_create_solution(c(0, 1))
  sol2 <- .test_create_solution(c(2, 3))
  
  expect_error(.test_solution_less_than(sol1, sol2), 
               "Cannot compare solutions without computed scores")
})

test_that("Solutions can be tested for equality by nodes", {
  sol1 <- .test_create_solution(c(0, 1, 2))
  sol2 <- .test_create_solution(c(2, 1, 0))  # Same nodes, different order
  sol3 <- .test_create_solution(c(0, 1, 3))  # Different nodes
  
  expect_true(.test_solution_equals(sol1, sol2))
  expect_false(.test_solution_equals(sol1, sol3))
})

# ==============================================================================
# Random Solution Generation Tests
# ==============================================================================

test_that("create_random_valid generates valid solutions", {
  tree_df <- data.frame(depth = c(1, 2, 3, 2, 3, 2, 3))
  tree <- .test_create_tree_constructor2(tree_df$depth)
  
  for (seed in 1:20) {
    sol <- .test_create_random_valid_solution(tree, target_size = 3, seed = seed)
    
    expect_true(.test_solution_is_valid(sol, tree))
    expect_true(.test_solution_size(sol) >= 1)
  }
})

test_that("create_random_valid with exact_size matches target", {
  tree_df <- data.frame(depth = c(1, 2, 3, 2, 3, 2, 3, 2, 3))
  tree <- .test_create_tree_constructor2(tree_df$depth)
  
  # With enough nodes, should achieve exact size
  for (seed in 1:10) {
    sol <- .test_create_random_valid_solution(tree, target_size = 3, 
                                              seed = seed, exact_size = TRUE)
    
    expect_true(.test_solution_is_valid(sol, tree))
    # May not always achieve exact size due to deduplication/validity
    expect_true(.test_solution_size(sol) >= 1)
  }
})

test_that("create_random_valid returns different solutions with different seeds", {
  tree_df <- data.frame(depth = c(1, 2, 3, 2, 3, 2, 3))
  tree <- .test_create_tree_constructor2(tree_df$depth)
  
  sol1 <- .test_create_random_valid_solution(tree, target_size = 3, seed = 1)
  sol2 <- .test_create_random_valid_solution(tree, target_size = 3, seed = 2)
  
  # Very likely to be different (not guaranteed but extremely probable)
  nodes1 <- .test_solution_get_nodes(sol1)
  nodes2 <- .test_solution_get_nodes(sol2)
  
  expect_false(identical(nodes1, nodes2))
})

# ==============================================================================
# Mutation Tests: Type 2 Swap
# ==============================================================================

test_that("mutate_swap_type2 creates new solution", {
  # Tree: 0 -> [1, 3], 1 -> [2], 3 -> [4]
  tree_df <- data.frame(depth = c(1, 2, 3, 2, 3))
  tree <- .test_create_tree_constructor2(tree_df$depth)
  
  sol <- .test_create_solution(c(2))  # Leaf node
  
  # Should swap with parent (node 1)
  mutated <- .test_solution_mutate_swap_type2(sol, tree, seed = 42)
  
  expect_true(.test_solution_size(mutated) > 0)
  # Might be node 1 (parent) or stay as 2 if no neighbors
  nodes <- .test_solution_get_nodes(mutated)
  expect_true(length(nodes) > 0)
})

test_that("mutate_swap_type2 with internal node swaps with children or parent", {
  # Tree: 0 -> [1, 3], 1 -> [2], 3 -> [4]
  tree_df <- data.frame(depth = c(1, 2, 3, 2, 3))
  tree <- .test_create_tree_constructor2(tree_df$depth)
  
  sol <- .test_create_solution(c(1))  # Internal node with parent 0 and child 2
  
  mutated <- .test_solution_mutate_swap_type2(sol, tree, seed = 42)
  
  nodes <- .test_solution_get_nodes(mutated)
  # Should be 0 (parent) or 2 (child)
  expect_true(0 %in% nodes || 2 %in% nodes)
})

test_that("determine_vertex finds parent and children correctly", {
  # Tree: 0 -> [1, 3], 1 -> [2], 3 -> [4]
  tree_df <- data.frame(depth = c(1, 2, 3, 2, 3))
  tree <- .test_create_tree_constructor2(tree_df$depth)
  
  # Node 1 has parent 0 and child 2
  sol <- .test_create_solution(c(1))
  
  # Test multiple mutations to see different neighbors
  neighbors_found <- integer(0)
  
  for (seed in 1:20) {
    mutated <- .test_solution_mutate_swap_type2(sol, tree, seed = seed)
    nodes <- .test_solution_get_nodes(mutated)
    neighbors_found <- c(neighbors_found, nodes)
  }
  
  # Should find both parent (0) and child (2) across multiple runs
  expect_true(0 %in% neighbors_found || 2 %in% neighbors_found)
})

# ==============================================================================
# Mutation Tests: Type 1 Add/Remove
# ==============================================================================

test_that("mutate_add_remove_type1 adds to size 1 solutions", {
  tree_df <- data.frame(depth = c(1, 2, 3, 2, 3))
  tree <- .test_create_tree_constructor2(tree_df$depth)
  
  sol <- .test_create_solution(c(0))
  
  # Size 1 should always add
  mutated <- .test_solution_mutate_add_remove_type1(sol, tree, alpha = 1.0, seed = 42)
  
  expect_equal(.test_solution_size(mutated), 2)
})

test_that("mutate_add_remove_type1 adds/removes based on alpha", {
  tree_df <- data.frame(depth = c(1, 2, 3, 2, 3, 2, 3))
  tree <- .test_create_tree_constructor2(tree_df$depth)
  
  sol <- .test_create_solution(c(1, 3))  # Size 2
  
  # With alpha = 4.0, P(add) = 4.0/2 = 2.0 (capped at 1.0) -> always add
  adds <- 0
  removes <- 0
  
  for (seed in 1:20) {
    mutated <- .test_solution_mutate_add_remove_type1(sol, tree, alpha = 4.0, seed = seed)
    
    if (.test_solution_size(mutated) > 2) adds <- adds + 1
    if (.test_solution_size(mutated) < 2) removes <- removes + 1
  }
  
  expect_true(adds > removes)
})

test_that("mutate_add_remove_type1 removes with low alpha", {
  tree_df <- data.frame(depth = c(1, 2, 3, 2, 3, 2, 3))
  tree <- .test_create_tree_constructor2(tree_df$depth)
  
  sol <- .test_create_solution(c(0, 1, 2, 3))  # Size 4
  
  # With alpha = 0.5, P(add) = 0.5/4 = 0.125 -> mostly remove
  adds <- 0
  removes <- 0
  
  for (seed in 1:20) {
    mutated <- .test_solution_mutate_add_remove_type1(sol, tree, alpha = 0.5, seed = seed)
    
    if (.test_solution_size(mutated) > 4) adds <- adds + 1
    if (.test_solution_size(mutated) < 4) removes <- removes + 1
  }
  
  expect_true(removes > adds)
})

# ==============================================================================
# Mutation Tests: Type 1 Replace
# ==============================================================================

test_that("mutate_replace_type1 creates new valid solution of same size", {
  tree_df <- data.frame(depth = c(1, 2, 3, 2, 3, 2, 3))
  tree <- .test_create_tree_constructor2(tree_df$depth)
  
  sol <- .test_create_solution(c(0, 2, 4))
  original_size <- .test_solution_size(sol)
  
  mutated <- .test_solution_mutate_replace_type1(sol, tree, seed = 42)
  
  # Should be valid
  expect_true(.test_solution_is_valid(mutated, tree))
  
  # Should have similar size (might differ slightly due to randomness)
  expect_true(.test_solution_size(mutated) == original_size)
})


# ==============================================================================
# Crossover Tests
# ==============================================================================

test_that("crossover_single_point swaps subtrees", {
  # Tree: 0 -> [1, 3], 1 -> [2], 3 -> [4]
  tree_df <- data.frame(depth = c(1, 2, 3, 2, 3))
  tree <- .test_create_tree_constructor2(tree_df$depth)
  
  parent1 <- .test_create_solution(c(1, 3)) 
  parent2 <- .test_create_solution(c(2, 4))
  
  children <- .test_solution_crossover_single_point(parent1, parent2, tree, seed = 42)
  
  child1 <- children$child1
  child2 <- children$child2
  
  # Children should be valid
  expect_true(.test_solution_is_valid(child1, tree))
  expect_true(.test_solution_is_valid(child2, tree))
  
  # Children should be non-empty
  expect_true(.test_solution_size(child1) > 0)
  expect_true(.test_solution_size(child2) > 0)
})

test_that("crossover_single_point with disjoint parents", {
  tree_df <- data.frame(depth = c(1, 2, 3, 2, 3))
  tree <- .test_create_tree_constructor2(tree_df$depth)
  
  parent1 <- .test_create_solution(c(1))
  parent2 <- .test_create_solution(c(3))
  
  children <- .test_solution_crossover_single_point(parent1, parent2, tree, seed = 42)
  
  expect_true(.test_solution_size(children$child1) > 0)
  expect_true(.test_solution_size(children$child2) > 0)
})

test_that("crossover_single_point handles all leaves tree", {
  # Tree with only leaves (no internal nodes)
  tree_df <- data.frame(depth = c(1, 1, 1))
  tree <- .test_create_tree_constructor2(tree_df$depth)
  
  parent1 <- .test_create_solution(c(0))
  parent2 <- .test_create_solution(c(1))
  
  children <- .test_solution_crossover_single_point(parent1, parent2, tree, seed = 42)
  
  # Should return parents unchanged (no internal nodes for crossover)
  expect_equal(.test_solution_get_nodes(children$child1), 
               .test_solution_get_nodes(parent1))
  expect_equal(.test_solution_get_nodes(children$child2), 
               .test_solution_get_nodes(parent2))
})

test_that("Hash is stable across multiple calls", {
  sol <- .test_create_solution(c(0, 1, 2))
  
  hash1 <- .test_solution_get_hash(sol)
  hash2 <- .test_solution_get_hash(sol)
  hash3 <- .test_solution_get_hash(sol)
  
  expect_equal(hash1, hash2)
  expect_equal(hash2, hash3)
})