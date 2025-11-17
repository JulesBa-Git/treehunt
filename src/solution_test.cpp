#include <Rcpp.h>
#include "solution.h"
#include "tree_structure.h"

// ==============================================================================
// Basic Constructors and Getters
// ==============================================================================

//' @keywords internal
// [[Rcpp::export(.test_create_solution)]]
SEXP test_create_solution(Rcpp::IntegerVector nodes) {
 std::vector<int> node_vec = Rcpp::as<std::vector<int>>(nodes);
 Solution* sol = new Solution(node_vec);
 return Rcpp::XPtr<Solution>(sol, true);
}

//' @keywords internal
// [[Rcpp::export(.test_solution_get_nodes)]]
Rcpp::IntegerVector test_solution_get_nodes(SEXP sol_ptr) {
 Rcpp::XPtr<Solution> sol(sol_ptr);
 return Rcpp::wrap(sol->get_nodes());
}

//' @keywords internal
// [[Rcpp::export(.test_solution_size)]]
int test_solution_size(SEXP sol_ptr) {
 Rcpp::XPtr<Solution> sol(sol_ptr);
 return sol->size();
}

//' @keywords internal
// [[Rcpp::export(.test_solution_empty)]]
bool test_solution_empty(SEXP sol_ptr) {
 Rcpp::XPtr<Solution> sol(sol_ptr);
 return sol->empty();
}

// ==============================================================================
// Score Management
// ==============================================================================

//' @keywords internal
// [[Rcpp::export(.test_solution_set_score)]]
void test_solution_set_score(SEXP sol_ptr, double score) {
 Rcpp::XPtr<Solution> sol(sol_ptr);
 sol->set_score(score);
}

//' @keywords internal
// [[Rcpp::export(.test_solution_get_score)]]
double test_solution_get_score(SEXP sol_ptr) {
 Rcpp::XPtr<Solution> sol(sol_ptr);
 return sol->get_score();
}

//' @keywords internal
// [[Rcpp::export(.test_solution_is_score_computed)]]
bool test_solution_is_score_computed(SEXP sol_ptr) {
 Rcpp::XPtr<Solution> sol(sol_ptr);
 return sol->is_score_computed();
}

//' @keywords internal
// [[Rcpp::export(.test_solution_invalidate_score)]]
void test_solution_invalidate_score(SEXP sol_ptr) {
 Rcpp::XPtr<Solution> sol(sol_ptr);
 sol->invalidate_score();
}

// ==============================================================================
// Hash Management
// ==============================================================================

//' @keywords internal
// [[Rcpp::export(.test_solution_get_hash)]]
double test_solution_get_hash(SEXP sol_ptr) {
 Rcpp::XPtr<Solution> sol(sol_ptr);
 // Return as double because R doesn't have size_t
 return static_cast<double>(sol->get_hash());
}

//' @keywords internal
// [[Rcpp::export(.test_solution_invalidate_hash)]]
void test_solution_invalidate_hash(SEXP sol_ptr) {
 Rcpp::XPtr<Solution> sol(sol_ptr);
 sol->invalidate_hash();
}

//' @keywords internal
// [[Rcpp::export(.test_solution_invalidate)]]
void test_solution_invalidate(SEXP sol_ptr) {
 Rcpp::XPtr<Solution> sol(sol_ptr);
 sol->invalidate();
}

// ==============================================================================
// Validation
// ==============================================================================

//' @keywords internal
// [[Rcpp::export(.test_solution_is_valid)]]
bool test_solution_is_valid(SEXP sol_ptr, SEXP tree_ptr) {
 Rcpp::XPtr<Solution> sol(sol_ptr);
 Rcpp::XPtr<tree_structure> tree(tree_ptr);
 return sol->is_valid(*tree);
}

// ==============================================================================
// Utility Functions
// ==============================================================================

//' @keywords internal
// [[Rcpp::export(.test_solution_contains)]]
bool test_solution_contains(SEXP sol_ptr, int node) {
 Rcpp::XPtr<Solution> sol(sol_ptr);
 return sol->contains(node);
}

//' @keywords internal
// [[Rcpp::export(.test_solution_print)]]
void test_solution_print(SEXP sol_ptr) {
 Rcpp::XPtr<Solution> sol(sol_ptr);
 sol->print();
}

// ==============================================================================
// Comparison Operators
// ==============================================================================

//' @keywords internal
// [[Rcpp::export(.test_solution_less_than)]]
bool test_solution_less_than(SEXP sol1_ptr, SEXP sol2_ptr) {
 Rcpp::XPtr<Solution> sol1(sol1_ptr);
 Rcpp::XPtr<Solution> sol2(sol2_ptr);
 return *sol1 < *sol2;
}

//' @keywords internal
// [[Rcpp::export(.test_solution_greater_than)]]
bool test_solution_greater_than(SEXP sol1_ptr, SEXP sol2_ptr) {
 Rcpp::XPtr<Solution> sol1(sol1_ptr);
 Rcpp::XPtr<Solution> sol2(sol2_ptr);
 return *sol1 > *sol2;
}

//' @keywords internal
// [[Rcpp::export(.test_solution_equals)]]
bool test_solution_equals(SEXP sol1_ptr, SEXP sol2_ptr) {
 Rcpp::XPtr<Solution> sol1(sol1_ptr);
 Rcpp::XPtr<Solution> sol2(sol2_ptr);
 return *sol1 == *sol2;
}

// ==============================================================================
// Random Solution Generation
// ==============================================================================

//' @keywords internal
// [[Rcpp::export(.test_create_random_valid_solution)]]
SEXP test_create_random_valid_solution(SEXP tree_ptr, int target_size, 
                                      int seed, bool exact_size = true) {
 Rcpp::XPtr<tree_structure> tree(tree_ptr);
 std::mt19937 rng(seed);
 
 Solution sol = Solution::create_random_valid(*tree, rng, target_size, 1000, exact_size);
 Solution* result = new Solution(sol);
 return Rcpp::XPtr<Solution>(result, true);
}

// ==============================================================================
// Mutation Operators
// ==============================================================================

//' @keywords internal
// [[Rcpp::export(.test_solution_mutate_swap_type2)]]
SEXP test_solution_mutate_swap_type2(SEXP sol_ptr, SEXP tree_ptr, int seed) {
 Rcpp::XPtr<Solution> sol(sol_ptr);
 Rcpp::XPtr<tree_structure> tree(tree_ptr);
 
 std::mt19937 rng(seed);
 Solution mutated = sol->mutate_swap_type2(*tree, rng);
 Solution* result = new Solution(mutated);
 return Rcpp::XPtr<Solution>(result, true);
}

//' @keywords internal
// [[Rcpp::export(.test_solution_mutate_add_remove_type1)]]
SEXP test_solution_mutate_add_remove_type1(SEXP sol_ptr, SEXP tree_ptr, 
                                          double alpha, int seed) {
 Rcpp::XPtr<Solution> sol(sol_ptr);
 Rcpp::XPtr<tree_structure> tree(tree_ptr);
 
 std::mt19937 rng(seed);
 Solution mutated = sol->mutate_add_remove_type1(*tree, alpha, rng);
 Solution* result = new Solution(mutated);
 return Rcpp::XPtr<Solution>(result, true);
}

//' @keywords internal
// [[Rcpp::export(.test_solution_mutate_replace_type1)]]
SEXP test_solution_mutate_replace_type1(SEXP sol_ptr, SEXP tree_ptr, int seed) {
 Rcpp::XPtr<Solution> sol(sol_ptr);
 Rcpp::XPtr<tree_structure> tree(tree_ptr);
 
 std::mt19937 rng(seed);
 Solution mutated = sol->mutate_replace_type1(*tree, rng);
 Solution* result = new Solution(mutated);
 return Rcpp::XPtr<Solution>(result, true);
}

// ==============================================================================
// Crossover Operators
// ==============================================================================

//' @keywords internal
 // [[Rcpp::export(.test_solution_crossover_single_point)]]
 Rcpp::List test_solution_crossover_single_point(SEXP sol1_ptr, SEXP sol2_ptr, 
                                                 SEXP tree_ptr, int seed) {
   Rcpp::XPtr<Solution> sol1(sol1_ptr);
   Rcpp::XPtr<Solution> sol2(sol2_ptr);
   Rcpp::XPtr<tree_structure> tree(tree_ptr);
   
   std::mt19937 rng(seed);
   auto [child1, child2] = Solution::crossover_single_point(*sol1, *sol2, *tree, rng);
   
   Solution* c1 = new Solution(child1);
   Solution* c2 = new Solution(child2);
   
   return Rcpp::List::create(
     Rcpp::Named("child1") = Rcpp::XPtr<Solution>(c1, true),
     Rcpp::Named("child2") = Rcpp::XPtr<Solution>(c2, true)
   );
 }