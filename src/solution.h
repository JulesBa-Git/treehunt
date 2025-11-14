#ifndef SOLUTION_H
#define SOLUTION_H

#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <random>
#include <functional>
#include <utility>
#include "tree_structure.h"

class Solution {
private:
  std::vector<int> selected_nodes_;
  double score_;                  
  bool score_computed_;           
  mutable size_t hash_;
  mutable bool hash_computed_;    
  
  size_t compute_hash_internal() const;
  
public:

  Solution();
  explicit Solution(std::vector<int> nodes);
  Solution(std::initializer_list<int> nodes);
  
  Solution(const Solution&) = default;
  Solution(Solution&&) = default;
  Solution& operator=(const Solution&) = default;
  Solution& operator=(Solution&&) = default;
  
  // Create random VALID solution (no parent-child pairs)
  // Keeps sampling until valid or max_attempts reached
  static Solution create_random_valid(const tree_structure& tree,
                                      std::mt19937& rng,
                                      size_t target_size,
                                      int max_attempts = 1000,
                                      bool exact_size = true);
  
  
  const std::vector<int>& get_nodes() const { 
    return selected_nodes_; 
  }
  
  size_t size() const { 
    return selected_nodes_.size(); 
  }
  
  bool empty() const {
    return selected_nodes_.empty();
  }
  
  double get_score() const {
    if (!score_computed_) {
      Rcpp::stop("Score not computed.");
    }
    return score_;
  }
  
  bool is_score_computed() const {
    return score_computed_;
  }
  
  size_t get_hash() const {
    if (!hash_computed_) {
      hash_ = compute_hash_internal();
      hash_computed_ = true;
    }
    return hash_;
  }
  
  void set_score(double score) {
    score_ = score;
    score_computed_ = true;
  }
  
  void invalidate_score() {
    score_computed_ = false;
  }
  
  void invalidate_hash() {
    hash_computed_ = false;
  }
  
  void invalidate() {
    score_computed_ = false;
    hash_computed_ = false;
  }
  
  bool is_valid(const tree_structure& tree) const;
  
  std::vector<std::pair<int,int>> determine_vertex(const tree_structure& tree) const;

  // Mutation operators
  Solution mutate_swap_type2(const tree_structure& tree, std::mt19937& rng) const;
  Solution mutate_add_remove_type1(const tree_structure& tree,
                                   double alpha, std::mt19937& rng) const;
  Solution mutate_replace_type1(const tree_structure& tree, std::mt19937& rng) const;
  
  // Crossover operators (static) ! TODO !
  static std::pair<Solution, Solution> crossover_single_point(
      const Solution& parent1,
      const Solution& parent2,
      std::mt19937& rng);
  
  // Comparison operators (for sorting purposes (by score))
  bool operator<(const Solution& other) const {
    if (!score_computed_ || !other.score_computed_) {
      Rcpp::stop("Cannot compare solutions without computed scores");
    } 
    return score_ < other.score_;
  } 
  
  bool operator>(const Solution& other) const {
    return other < *this;
  }
  
  bool operator<=(const Solution& other) const {
    return !(other < *this);
  }
  
  bool operator>=(const Solution& other) const {
    return !(*this < other);
  }
  
  // Equality is base on nodes
  bool operator==(const Solution& other) const {
    return selected_nodes_ == other.selected_nodes_;
  }
  
  bool operator!=(const Solution& other) const {
    return !(*this == other);
  }

  void print() const;
  
  bool contains(int node) const {
    return std::find(selected_nodes_.begin(), selected_nodes_.end(), node)
    != selected_nodes_.end();
  }
};

#endif