#include "mcmc_algorithm.h"
#include <set>
#include <climits>

template<typename TargetType>
MCMCAlgorithm<TargetType>::MCMCAlgorithm(const PatientData<TargetType>& data,
                                         const MCMCParams& params,
                                         const PWPScoreContext *pwp_context)
  : data_{data}, pwp_context_{pwp_context}, params_{params},
    min_score_(0.0), min_score_filtered_(0.0) {
  
  if (params_.epochs == 0 || params_.epochs > INT_MAX)
    Rcpp::stop("epochs must be a positive integer no greater than .Machine$integer.max.");
  if (!std::isfinite(params_.temperature) || params_.temperature <= 0)
    Rcpp::stop("temperature must be finite and positive.");
  if (!std::isfinite(params_.prob_mutation_type1) ||
      params_.prob_mutation_type1 <= 0 || params_.prob_mutation_type1 > 1)
    Rcpp::stop("prob_type1 must be finite and in (0, 1] to connect the state space.");
  if (!std::isfinite(params_.max_score) || params_.max_score <= 0 ||
      params_.max_score > (INT_MAX - 1) / 10.0)
    Rcpp::stop("max_score must be finite, positive, and small enough for histogram indexing.");
  if (params_.cocktail_size == 0 || params_.cocktail_size > data_.get_tree().size())
    Rcpp::stop("cocktail_size must lie between 1 and the number of tree nodes.");
  if (params_.burn_in >= params_.epochs)
    Rcpp::stop("burn_in must be smaller than epochs.");
  size_t n_bins = static_cast<size_t>(std::ceil(params_.max_score * 10)) + 1;
  results_.score_distribution.resize(n_bins, 0);
  results_.score_distribution_filtered.resize(n_bins, 0);
  results_.uniform_reference.resize(n_bins);
  results_.uniform_reference_filtered.resize(n_bins);
  results_.cocktail_size = params_.cocktail_size;
  if (params.seed >= 0) {
    rng_ = std::mt19937(static_cast<unsigned int>(params.seed));
  } else {
    std::random_device rd;
    std::vector<unsigned int> seed_data;
    for (int i = 0; i < 4; ++i) {
      seed_data.push_back(rd());
    }
    std::seed_seq seq(seed_data.begin(), seed_data.end());
    rng_ = std::mt19937(seq);
  }
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
  // Floyd sampling: every subset of k distinct nodes has the same probability.
  const size_t n = data_.get_tree().size();
  std::set<int> selected;
  for (size_t j = n - params_.cocktail_size; j < n; ++j) {
    const int candidate = std::uniform_int_distribution<int>(0, j)(rng_);
    if (!selected.insert(candidate).second) selected.insert(static_cast<int>(j));
  }
  return Solution(std::vector<int>(selected.begin(), selected.end()));
}

template<typename TargetType>
Solution MCMCAlgorithm<TargetType>::initial_solution() {
  // Build a supported state directly, avoiding unbounded rejection at startup.
  const auto& depth = data_.get_tree().get_depth();
  std::vector<int> parent(depth.size()), ancestors(data_.get_tree().max_depth() + 1, -1);
  for (size_t node = 0; node < depth.size(); ++node) {
    parent[node] = ancestors[depth[node] - 1];
    ancestors[depth[node]] = static_cast<int>(node);
  }
  for (size_t row = 0; row < data_.size(); ++row) {
    if (row % 1024 == 0) Rcpp::checkUserInterrupt();
    std::set<int> supported;
    for (int node : data_.get_patient_nodes(row)) {
      if (node < 0 || static_cast<size_t>(node) >= depth.size())
        Rcpp::stop("Observation nodes must be zero-based indices within the tree.");
      for (int ancestor = node; ancestor >= 0; ancestor = parent[ancestor])
        supported.insert(ancestor);
    }
    if (supported.size() >= params_.cocktail_size) {
      std::vector<int> nodes(supported.begin(), supported.end());
      std::shuffle(nodes.begin(), nodes.end(), rng_);
      nodes.resize(params_.cocktail_size);
      return Solution(nodes);
    }
  }
  Rcpp::stop("No observed combination has the requested cocktail_size.");
}

template<typename TargetType>
std::vector<std::pair<int, int>> MCMCAlgorithm<TargetType>::local_moves(
    const Solution& current) const {
  auto moves = current.determine_vertex(data_.get_tree());
  const auto& nodes = current.get_nodes();
  moves.erase(std::remove_if(moves.begin(), moves.end(), [&](const auto& move) {
    return std::binary_search(nodes.begin(), nodes.end(), move.second);
  }), moves.end());
  return moves;
}

template<typename TargetType>
double MCMCAlgorithm<TargetType>::capped_score(double score) const {
  // Positive infinity is supported by the existing finite score cap.
  if (std::isnan(score) || score < 0)
    Rcpp::stop("MCMC requires non-negative scores; use positive_only = TRUE for PWP.");
  return std::min(score, params_.max_score);
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
    
  case ScoreType::WILCOXON :
    return ScoreFunctions<TargetType>::compute_wilcoxon_risk_with_data(
      data_, sol);
  
  case ScoreType::COMPOSITE :
    return ScoreFunctions<TargetType>::compute_multifactor_risk_QT_with_data(
      data_, sol);
    
  case ScoreType::RESIDUALS :
    return ScoreFunctions<TargetType>::compute_residuals_risk_with_stats(
      data_, sol).first;

  case ScoreType::PWP_RAO : {
    if (pwp_context_ == nullptr)
      Rcpp::stop("PWP Rao scoring requires a PWP score context.");
    const auto pwp = pwp_context_->compute(sol);
    typename ScoreFunctions<TargetType>::ScoreData out;
    out.score = pwp.fitness;
    out.covered_patients = pwp.covered_patients;
    out.covered_nonzero_target = pwp.covered_events;
    return out;
  }
    
  default:
    Rcpp::stop("Unknown score type");
  
  }  
}

template<typename TargetType>
void MCMCAlgorithm<TargetType>::update_score_distribution(
    double score, size_t covered_patients, bool exhaustive) {
  const size_t bin = score < params_.max_score ?
    static_cast<size_t>(score * 10.0) : results_.score_distribution.size() - 1;
  ++results_.score_distribution[bin];
  if (score >= params_.max_score) results_.outstanding_scores.push_back(score);
  // Exhaustive enumeration already follows the uniform reference.
  results_.uniform_reference.add(bin, exhaustive ? 0.0 : score, params_.temperature);
  if (covered_patients > params_.beta) {
    ++results_.score_distribution_filtered[bin];
    results_.uniform_reference_filtered.add(bin, exhaustive ? 0.0 : score,
                                             params_.temperature);
  }
}

template<typename TargetType>
void MCMCAlgorithm<TargetType>::update_top_solutions(const Solution& sol, 
                                                     double score, 
                                                     int covered_patients) {
  
  if (params_.n_results == 0) return;
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
  Solution current = initial_solution();

  auto current_score_data = compute_score(current);
  double current_score = capped_score(current_score_data.score);
  
  std::uniform_real_distribution<double> uniform(0.0, 1.0);
  
  
  for (size_t epoch = 0; epoch < params_.epochs; ++epoch) {
    if (epoch % 1024 == 0) Rcpp::checkUserInterrupt();
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
        double proposed_score = capped_score(proposed_score_data.score);
        
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
      auto current_vertices = local_moves(current);
      
      if (current_vertices.empty()) {
        ++results_.type2_moves;
        ++results_.rejected_moves;
      } else {
        const auto move = current_vertices[std::uniform_int_distribution<size_t>(
          0, current_vertices.size() - 1)(rng_)];
        auto nodes = current.get_nodes();
        *std::find(nodes.begin(), nodes.end(), move.first) = move.second;
        Solution proposed(nodes);
        ++results_.type2_moves;
        
        if (!is_in_population(proposed)) {
          ++results_.proposals_not_in_population;
          ++results_.rejected_moves;
        } else {
          ++results_.type2_in_population;
          
          auto proposed_score_data = compute_score(proposed);
          double proposed_score = capped_score(proposed_score_data.score);
          
          auto proposed_vertices = local_moves(proposed);
          
          double acceptance_prob = compute_acceptance_probability_type2(
            current_score, proposed_score,
            current_vertices.size(), proposed_vertices.size());
          
          if (uniform(rng_) < acceptance_prob) {
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
    
    if (epoch >= params_.burn_in) {
      update_score_distribution(current_score, current_score_data.covered_patients);
      update_top_solutions(current, current_score, current_score_data.covered_patients);
      if (params_.store_trace) {
        results_.trace_scores.push_back(current_score);
        results_.trace_support.push_back(current_score_data.covered_patients);
        results_.trace_solutions.push_back(current.get_nodes());
      }
    }
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

template<typename TargetType>
MCMCResults MCMCAlgorithm<TargetType>::true_size2_distribution() {
  
  const size_t tree_size = data_.get_tree().size();
  results_.exhaustive = true;
  for (size_t i = 0; i + 1 < tree_size; ++i) {
    Rcpp::checkUserInterrupt();
    for (size_t j = i + 1; j < tree_size; ++j) {
      Solution explorer({static_cast<int>(i), static_cast<int>(j)});
      if (!is_in_population(explorer)) continue;
      const auto scored = compute_score(explorer);
      update_score_distribution(capped_score(scored.score), scored.covered_patients, true);
      ++results_.total_iterations;
    }
  }
  finalize_results();
  return results_;
}


template class MCMCAlgorithm<int>;
template class MCMCAlgorithm<double>;
