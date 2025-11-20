#include "mcmc_algorithm.h"

template<typename TargetType>
MCMCAlgorithm<TargetType>::MCMCAlgorithm(const PatientData<TargetType>& data,
                                         const MCMCParams& params)
  : data_{data}, params_{params},
    min_score_(0.0), min_score_filtered_(0.0) {
  
  size_t n_bins = static_cast<size_t>(params_.max_score * 10) + 1;
  results_.score_distribution.resize(n_bins, 0);
  results_.score_distribution_filtered.resize(n_bins, 0);
  results_.cocktail_size = params_.cocktail_size;
  std::random_device rd;
  std::vector<unsigned int> seed_data;
  for (int i = 0; i < 4; ++i) {
    seed_data.push_back(rd());
  }
  
  std::seed_seq seq(seed_data.begin(), seed_data.end());
  
  rng_ = std::mt19937(seq);
}

template<typename TargetType>
bool MCMCAlgorithm<TargetType>::is_in_population(const Solution& sol) const {
  const auto& nodes = sol.get_nodes();
  
  for (size_t i = 0; i < data_.size(); ++i) {
    if (data_.patient_has_combination(i, nodes)) {
      return true;
    }
  }
  
  return false;
}

template<typename TargetType>
Solution MCMCAlgorithm<TargetType>::propose_type1_mutation() {
  // Type 1: Generate completely new random solution
  return Solution::create_random_valid(data_.get_tree(), rng_, params_.cocktail_size, 
                                       100, true);
}

template<typename TargetType>
Solution MCMCAlgorithm<TargetType>::propose_type2_mutation(const Solution& current) {
  // Type 2: Swap one node with its parent or child
  return current.mutate_swap_type2(data_.get_tree(), rng_);
}

template<typename TargetType>
double MCMCAlgorithm<TargetType>::compute_acceptance_probability_type1(
    double current_score, double proposed_score) {
  return std::exp((proposed_score - current_score) / params_.temperature);
}

template<typename TargetType>
double MCMCAlgorithm<TargetType>::compute_acceptance_probability_type2(
    double current_score, double proposed_score,
    size_t n_vertex_current, size_t n_vertex_proposed) {
  double score_ratio = std::exp((proposed_score - current_score) / params_.temperature);
  double vertex_ratio = static_cast<double>(n_vertex_current) / 
    static_cast<double>(n_vertex_proposed); 
  
  return score_ratio * vertex_ratio;
}

template<typename TargetType>
typename ScoreFunctions<TargetType>::ScoreData 
MCMCAlgorithm<TargetType>::compute_score(const Solution& sol) const{
  switch(params_.score_type_){
  case ScoreType::HYPERGEOMETRIC :
    return ScoreFunctions<TargetType>::compute_hypergeometric_with_data(
      data_, sol, true);
  
  case ScoreType::RELATIVE_RISK :
    return ScoreFunctions<TargetType>::compute_relative_risk_with_data(
      data_, sol);
  
  default:
    Rcpp::stop("Unknown score type");
  
  }  
}

template<typename TargetType>
void MCMCAlgorithm<TargetType>::update_score_distribution(double score, 
                                                          size_t covered_patients) {
  if (score < params_.max_score) {
    size_t bin_index = static_cast<size_t>(score * 10.0);
    ++results_.score_distribution[bin_index];
    
    if (covered_patients > params_.beta) {
      ++results_.score_distribution_filtered[bin_index];
    }
  } else {
    results_.outstanding_scores.push_back(score);
    ++results_.score_distribution.back(); 
    
    if (covered_patients > params_.beta) {
      ++results_.score_distribution_filtered.back();
    }
  }
}

template<typename TargetType>
void MCMCAlgorithm<TargetType>::update_top_solutions(const Solution& sol, 
                                                     double score, 
                                                     int covered_patients) {
  
  // Update regular top solutions
  if (top_heap_.size() < params_.n_results) {
    top_heap_.push(std::make_pair(score, sol.get_nodes()));
    min_score_ = top_heap_.top().first;
  } else if (score > min_score_) {
    top_heap_.pop();
    top_heap_.push(std::make_pair(score, sol.get_nodes()));
    min_score_ = top_heap_.top().first;
  }
  
  // Update filtered top solutions (beta threshold)
  if (covered_patients > params_.beta) {
    if (top_heap_filtered_.size() < params_.n_results) {
      top_heap_filtered_.push(std::make_pair(score, sol.get_nodes()));
      min_score_filtered_ = top_heap_filtered_.top().first;
    } else if (score > min_score_filtered_) {
      top_heap_filtered_.pop();
      top_heap_filtered_.push(std::make_pair(score, sol.get_nodes()));
      min_score_filtered_ = top_heap_filtered_.top().first;
    }
  }
}

template<typename TargetType>
void MCMCAlgorithm<TargetType>::finalize_results() {
  // Extract from heaps into result vectors
  
  // Regular top solutions
  std::vector<std::pair<double, std::vector<int>>> temp;
  temp.reserve(top_heap_.size());
  
  while (!top_heap_.empty()) {
    temp.push_back(top_heap_.top());
    top_heap_.pop();
  }
  std::reverse(temp.begin(), temp.end());
  
  results_.top_solutions.reserve(temp.size());
  results_.top_scores.reserve(temp.size());
  
  for (const auto& item : temp) {
    results_.top_scores.push_back(item.first);
    results_.top_solutions.push_back(item.second);
  }
  
  // Filtered top solutions
  temp.clear();
  temp.reserve(top_heap_filtered_.size());
  
  while (!top_heap_filtered_.empty()) {
    temp.push_back(top_heap_filtered_.top());
    top_heap_filtered_.pop();
  }
  
  std::reverse(temp.begin(), temp.end());
  
  results_.top_solutions_filtered.reserve(temp.size());
  results_.top_scores_filtered.reserve(temp.size());
  
  for (const auto& item : temp) {
    results_.top_scores_filtered.push_back(item.first);
    results_.top_solutions_filtered.push_back(item.second);
  }
}

template<typename TargetType>
MCMCResults MCMCAlgorithm<TargetType>::run() {
  Solution current = Solution::create_random_valid(data_.get_tree(), rng_, params_.cocktail_size,
                                                   1000, true);
  
  while (!is_in_population(current)) {
    current = Solution::create_random_valid(data_.get_tree(), rng_, params_.cocktail_size,
                                            1000, true);
  }
  
  auto current_score_data = compute_score(current);
  double current_score = std::min(current_score_data.score, params_.max_score);
  
  std::uniform_real_distribution<double> uniform(0.0, 1.0);
  
  
  for (size_t epoch = 0; epoch < params_.epochs; ++epoch) {
    bool is_type1 = uniform(rng_) < params_.prob_mutation_type1;
    
    if (is_type1) {
      Solution proposed = propose_type1_mutation();
      ++results_.type1_moves;
      
      if (!is_in_population(proposed)) {
        ++results_.proposals_not_in_population;
        ++results_.rejected_moves;
      } else {
        ++results_.type1_in_population;
        
        auto proposed_score_data = compute_score(proposed);
        double proposed_score = std::min(proposed_score_data.score, params_.max_score);
        
        double acceptance_prob = compute_acceptance_probability_type1(current_score, 
                                                                      proposed_score);
        
        if (uniform(rng_) < acceptance_prob) {
          current = proposed;
          current_score = proposed_score;
          current_score_data = proposed_score_data;
          ++results_.accepted_moves;
          ++results_.type1_accepted;
        } else {
          ++results_.rejected_moves;
        }
      }
      
    } else {
      // Type 2 mutation
      auto current_vertices = current.determine_vertex(data_.get_tree());
      
      if (current_vertices.empty()) {
        ++results_.type2_moves;
        ++results_.rejected_moves;
      } else {
        Solution proposed = propose_type2_mutation(current);
        ++results_.type2_moves;
        
        if (!is_in_population(proposed)) {
          ++results_.proposals_not_in_population;
          ++results_.rejected_moves;
        } else {
          ++results_.type2_in_population;
          
          auto proposed_score_data = compute_score(proposed);
          double proposed_score = std::min(proposed_score_data.score, params_.max_score);
          
          auto proposed_vertices = proposed.determine_vertex(data_.get_tree());
          
          double acceptance_prob = compute_acceptance_probability_type2(
            current_score, proposed_score,
            current_vertices.size(), proposed_vertices.size());
          
          if (uniform(rng_) < acceptance_prob) {
            // Accept
            current = proposed;
            current_score = proposed_score;
            current_score_data = proposed_score_data;
            ++results_.accepted_moves;
            ++results_.type2_accepted;
          } else {
            ++results_.rejected_moves;
          }
        }
      }
    }
    
    // Update distributions and top solutions with current state
    update_score_distribution(current_score, current_score_data.covered_patients);
    update_top_solutions(current, current_score, current_score_data.covered_patients);
  }
  
  results_.total_iterations = params_.epochs;
  finalize_results();
  
  if (params_.verbose) {
    Rcpp::Rcout << "=== MCMC Run Statistics ===\n";
    Rcpp::Rcout << "Total iterations: " << results_.total_iterations << "\n";
    Rcpp::Rcout << "Accepted moves: " << results_.accepted_moves 
                << " (" << (100.0 * results_.accepted_moves / results_.total_iterations) 
                << "%)\n";
    Rcpp::Rcout << "Type 1 moves: " << results_.type1_moves 
                << " (accepted: " << results_.type1_accepted << ")\n";
    Rcpp::Rcout << "Type 2 moves: " << results_.type2_moves 
                << " (accepted: " << results_.type2_accepted << ")\n";
    Rcpp::Rcout << "Proposals not in population: " << results_.proposals_not_in_population << "\n";
    
    if (results_.type1_in_population > 0) {
      Rcpp::Rcout << "Type 1 acceptance (in pop): " 
                  << (100.0 * results_.type1_accepted / results_.type1_in_population)
                  << "%\n";
    }
    if (results_.type2_in_population > 0) {
      Rcpp::Rcout << "Type 2 acceptance (in pop): "
                  << (100.0 * results_.type2_accepted / results_.type2_in_population)
                  << "%\n";
    }
  }
  
  return results_;
}

template class MCMCAlgorithm<int>;
template class MCMCAlgorithm<double>;