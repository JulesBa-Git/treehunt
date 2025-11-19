#ifndef MCMC_ALGORITHM_H
#define MCMC_ALGORITHM_H

#include <Rcpp.h>
#include <vector>
#include <queue>
#include <random>
#include <algorithm>
#include <utility>
#include "solution.h"
#include "patient_data.h"
#include "tree_structure.h"
#include "score_functions.h"

// Parameters structure for MCMC
struct MCMCParams {
  size_t epochs;
  double temperature;
  size_t n_results;         
  size_t cocktail_size;     
  double prob_mutation_type1;
  size_t beta;              
  double max_score;         
  int seed;
  bool verbose;
  
  MCMCParams() : epochs(10000), temperature(1.0), n_results(10), 
  cocktail_size(3), prob_mutation_type1(0.5), 
  beta(5), max_score(100.0), seed(42), verbose(false) {}
};

// Results structure
struct MCMCResults {
  // Score distributions (0.1-wide bins)
  std::vector<size_t> score_distribution;
  std::vector<size_t> score_distribution_filtered;
  std::vector<double> outstanding_scores;  // Scores > max_score
  
  // Top solutions
  std::vector<std::vector<int>> top_solutions;
  std::vector<double> top_scores;
  
  // Top solutions with beta filter
  std::vector<std::vector<int>> top_solutions_filtered;
  std::vector<double> top_scores_filtered;
  
  // Statistics
  size_t total_iterations;
  size_t accepted_moves;
  size_t rejected_moves;
  size_t proposals_not_in_population;
  size_t type1_moves;
  size_t type2_moves;
  size_t type1_accepted;
  size_t type2_accepted;
  size_t type1_in_population;
  size_t type2_in_population;
  size_t cocktail_size;
  
  MCMCResults() : total_iterations(0), accepted_moves(0), rejected_moves(0),
  proposals_not_in_population(0), type1_moves(0), type2_moves(0),
  type1_accepted(0), type2_accepted(0), 
  type1_in_population(0), type2_in_population(0),
  cocktail_size(0) {}
};

template<typename TargetType>
class MCMCAlgorithm {
private:
  const PatientData<TargetType>& data_;
  const tree_structure& tree_;
  MCMCParams params_;
  std::mt19937 rng_;
  MCMCResults results_;
  
  // Min-heap comparator for tracking top solutions
  struct SolutionComparator {
    bool operator()(const std::pair<double, std::vector<int>>& a,
                 const std::pair<double, std::vector<int>>& b) const {
      return a.first > b.first;
    }
  };
  
  // Priority queues for top solutions
  std::priority_queue<
    std::pair<double, std::vector<int>>,
    std::vector<std::pair<double, std::vector<int>>>,
    SolutionComparator
  > top_heap_;
  
  std::priority_queue<
    std::pair<double, std::vector<int>>,
    std::vector<std::pair<double, std::vector<int>>>,
    SolutionComparator
  > top_heap_filtered_;
  
  double min_score_;           
  double min_score_filtered_; 
  
  bool is_in_population(const Solution& sol) const;
  Solution propose_type1_mutation();
  Solution propose_type2_mutation(const Solution& current);
  double compute_acceptance_probability_type1(double current_score, double proposed_score);
  double compute_acceptance_probability_type2(double current_score, double proposed_score,
                                             size_t n_vertex_current, size_t n_vertex_proposed);
  void update_score_distribution(double score, size_t covered_patients);
  void update_top_solutions(const Solution& sol, double score, int covered_patients);
  void finalize_results();
  
public:
  MCMCAlgorithm(const PatientData<TargetType>& data,
                const tree_structure& tree,
                const MCMCParams& params);
  
  MCMCResults run();
  
  const MCMCResults& get_results() const { return results_; }
};

#endif