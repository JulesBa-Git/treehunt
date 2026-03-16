#include "genetic_algorithm.h"

template<typename TargetType>
GeneticAlgorithm<TargetType>::GeneticAlgorithm(const PatientData<TargetType>& data,
                                               const GAParams& params,
                                               const Rcpp::Nullable<Rcpp::List>& seed_population)
  : population_{}, data_{data}, params_{params}, cache_hits_{0} {
    
    population_.reserve(params.population_size);
    std::random_device rd;
    std::vector<unsigned int> seed_data;
    for (int i = 0; i < 4; ++i) {
      seed_data.push_back(rd());
    }
    
    std::seed_seq seq(seed_data.begin(), seed_data.end());
    
    rng_ = std::mt19937(seq);
    
    if(seed_population.isNotNull()){
      Rcpp::List seed(seed_population);
      if(seed.size() > params.population_size)
        Rcpp::stop("Seed population size (%i) exceeds the total population size (%i)." 
                     ,seed.size(), params.population_size);
      
      for(int i = 0; i < seed.size(); ++i){
        std::vector<int> nodes = Rcpp::as<std::vector<int>>(seed[i]);
        
        // subtract 1 from every element (1-based to 0-based)
        std::transform(nodes.begin(), nodes.end(), nodes.begin(), 
                       [](int x) { return x - 1; });
        
        Solution tmp_sol(std::move(nodes));
        
        if(tmp_sol.is_valid(data.get_tree()))
          population_.emplace_back(std::move(tmp_sol));
        else
          Rcpp::Rcout<< "Warning: Seed individual at index " << i << " is invalid and was skipped.\n";
      }
    }
}

template<typename TargetType>
void GeneticAlgorithm<TargetType>::initialize(){
  std::uniform_int_distribution<size_t> node(0, data_.get_tree().get_depth().size()-1);
  std::poisson_distribution<int> solution_size(data_.mean_patient_nodes());
  
  for(size_t i = population_.size(); i < params_.population_size; ++i){
    int solution_length = solution_size(rng_);
    if (solution_length == 0)
      solution_length = 1;
    population_.emplace_back(Solution::create_random_valid(data_.get_tree(), 
                                                 rng_,
                                                 solution_length));
  }
  population_.shrink_to_fit();
  assert(population_.size() == params_.population_size);
}

template<typename TargetType>
double GeneticAlgorithm<TargetType>::compute_score(const Solution& sol) const{
  switch(params_.score_type){
  case ScoreType::HYPERGEOMETRIC :
    return ScoreFunctions<TargetType>::compute_hypergeometric(data_, sol);
  
  case ScoreType::RELATIVE_RISK :
    return ScoreFunctions<TargetType>::compute_relative_risk(data_, sol);
    
  case ScoreType::WILCOXON : 
    return ScoreFunctions<TargetType>::compute_wilcoxon_risk(data_, sol);
  
  default :
    Rcpp::stop("Unknown score type");
  }
}

//try to put a caching system. Queue might become too large, set a limit?
template<typename TargetType>
void GeneticAlgorithm<TargetType>::evaluate(){
  
  for(Solution& solution : population_){
    auto search = score_hash_map_.find(solution.get_hash());
    if(search != score_hash_map_.end()){
      ++cache_hits_;
      solution.set_score(search->second);
    }
    else{
      solution.set_score(compute_score(solution));
      score_hash_map_.insert({solution.get_hash(), solution.get_score()});
    }
  }
}

template<typename TargetType>
void GeneticAlgorithm<TargetType>::penalize(){
  auto [M, index] = dissimilarity_input_computation();
  
  std::vector<std::vector<double>> 
    similarity(population_.size(), std::vector<double>(population_.size(), -1.0));

  for(size_t i = 0; i < population_.size() -1 ; ++i){
    similarity[i][i] = 1;
    for(size_t j = i+1; j < population_.size(); ++j){
      
      if(similarity[i][j] < 0){
        double i_j_similarity = 1.0 -
          normalized_distance(i, j, M, index, data_.get_tree().max_depth());
        
        similarity[i][j] = i_j_similarity;
        similarity[j][i] = i_j_similarity;
      }
    }
  }
  similarity[population_.size()-1][population_.size()-1] = 1;
  
  for(size_t i = 0 ; i < population_.size(); ++i){
    //Accumulation starts at -1 because the similarity with himself is 1
    double penalized_score = 
      population_[i].get_score() / std::accumulate(similarity[i].begin(),
                                                 similarity[i].end(), -1.0);
    population_[i].set_score(penalized_score);
  }
}

template<typename TargetType>
void GeneticAlgorithm<TargetType>::keep_elite(std::vector<Solution>& offspring) const{
  std::vector<std::reference_wrapper<const Solution>> sorted_pop(
      population_.begin(), population_.end());
  
  std::partial_sort(sorted_pop.begin(), 
                    sorted_pop.begin() + params_.elite_count,
                    sorted_pop.end(),
                    [](const Solution& a, const Solution& b) {
                      return a > b; 
                    });
  
  for(size_t i = 0; i < params_.elite_count; ++i){
    offspring.emplace_back(sorted_pop[i].get());
  }
}

template<typename TargetType>
void GeneticAlgorithm<TargetType>::tournament_selection(
    std::vector<Solution>& offspring) const{
  const size_t solutions_to_select = population_.size() - params_.elite_count;
  
  std::uniform_int_distribution<size_t> index_dist(0, population_.size()-1);
  
  for(size_t i = 0; i < solutions_to_select; ++i){
    
    size_t best_index = index_dist(rng_);
    for(int j = 1 ; j < params_.tournament_size; ++j){
      const size_t current_index = index_dist(rng_);
      
      if(population_[current_index] > population_[best_index]){
        best_index = current_index;
      }
    } 
    offspring.emplace_back(population_[best_index]);
  }
}

template<typename TargetType>
void GeneticAlgorithm<TargetType>::mutate(std::vector<Solution>& offspring) const{
  std::uniform_real_distribution<double> prob_dist(0.0, 1.0);
  for(size_t i = params_.elite_count; i < offspring.size(); ++i){
    double mutation_probability = prob_dist(rng_);
    
    if(mutation_probability < params_.mutation_rate){
      Solution tmp = offspring[i].mutate_genetic_algorithm(
        data_.get_tree(), 
        params_.alpha,
        rng_
      );
      if(tmp.is_valid(data_.get_tree())){
        offspring[i]= std::move(tmp);
      }
    }
  }
}

template<typename TargetType>
void GeneticAlgorithm<TargetType>::crossover(std::vector<Solution>& offspring) const{
  std::uniform_real_distribution<double> prob_dist(0.0, 1.0);
  
  for(size_t i = params_.elite_count; i +1 < offspring.size(); i+=2){
    if(prob_dist(rng_) < params_.crossover_rate){
      auto [c1,c2] = Solution::crossover_single_point(offspring[i],
                                                      offspring[i+1],
                                                      data_.get_tree(),
                                                      rng_);
      offspring[i] = std::move(c1);
      offspring[i+1] = std::move(c2);
    }
  }
}

template<typename TargetType>
void GeneticAlgorithm<TargetType>::replace_population(std::vector<Solution>& offspring){
  population_ = std::move(offspring);
  offspring.clear();
  offspring.reserve(params_.population_size);
}

template<typename TargetType>
GAResults GeneticAlgorithm<TargetType>::extract_results() const {
  GAResults results;
  results.cache_hits = cache_hits_;
  results.total_generations = params_.epochs;
  
  // Final population
  results.final_population.reserve(population_.size());
  results.final_scores.reserve(population_.size());
  
  for(const auto& sol : population_){
    results.final_population.push_back(sol.get_nodes());
    results.final_scores.push_back(sol.get_score());
  }
  
  return results;
}

template<typename TargetType>
GAResults GeneticAlgorithm<TargetType>::run(){
  initialize();

  std::vector<Solution> offspring;
  offspring.reserve(params_.population_size);

  for(size_t i = 0 ; i < params_.epochs; ++i){
    evaluate();
    
    if(params_.diversity)
      penalize();
      
    
    if(params_.verbose && (i % 10 == 0 || i == params_.epochs - 1)){
      double best_score = 0.0;
      for(const auto& sol : population_){
        if(sol.get_score() > best_score){
          best_score = sol.get_score();
        }
      } 
      Rcpp::Rcout << "Generation " << i << "/" << params_.epochs 
                  << " - Best score: " << best_score 
                  << " - Cache size: " << score_hash_map_.size()
                  << " - Cache hits: " << cache_hits_ << "\n";
    } 
    
    if(params_.elite_count > 0){
      keep_elite(offspring);
    }
    
    tournament_selection(offspring);
    mutate(offspring);
    crossover(offspring);
    
    replace_population(offspring);
  }
  
  // Final evaluation
  evaluate();
  
  GAResults results = extract_results();
  
  if(params_.verbose){
    Rcpp::Rcout << "\n=== GA Run Complete ===\n";
    Rcpp::Rcout << "Total generations: " << results.total_generations << "\n";
    Rcpp::Rcout << "Cache hits: " << results.cache_hits << "\n";
  }
  
  return results;
  
}

template<typename TargetType>
std::pair<std::vector<std::vector<int>>, std::vector<int>> 
  GeneticAlgorithm<TargetType>::dissimilarity_input_computation() const {
  std::vector<std::vector<int>> M;
  std::vector<int> index;
  
  const auto& tree = data_.get_tree();
  const auto& father = tree.get_father();
  const auto& depth = tree.get_depth();
  
  int node_count = 0;
  
  for(const auto& nodes : population_){
    node_count+= nodes.size();
  }
  
  index.reserve(population_.size());
  M.resize(node_count);
  
  for(size_t i = 0; i < M.size(); ++i)
    M[i].reserve(tree.max_depth() + 1);
  
  size_t row = 0;
  size_t current_patient_node_index = 0;
  for(const Solution& sol : population_){
    const auto& nodes_vector = sol.get_nodes();
    for(const auto& node : nodes_vector){
      //if we are not on max_depth, we should push the node on the lower
      //depth 
      for(int i = tree.max_depth(); i >= depth[node]; --i){
        //add the node
        M[row].push_back(node);
      }
      int current_father = father[node];
      //add the father of this node
      while(depth[current_father] >= 1){
        M[row].push_back(current_father);
        current_father = father[current_father];
      }// push the father of depth 1
      M[row++].push_back(current_father);
    }
    index.push_back(current_patient_node_index);
    current_patient_node_index += nodes_vector.size();
  }
  
  return std::make_pair(std::move(M), std::move(index));
}

template<typename TargetType>
 Rcpp::NumericMatrix GeneticAlgorithm<TargetType>::population_dissimilarity() const{
   auto [M, index] = dissimilarity_input_computation();
   size_t n = population_.size();
   Rcpp::NumericMatrix dissimilarity(n, n);
   dissimilarity.fill(-1.0);
   
   for(size_t i = 0; i < n -1 ; ++i){ 
     dissimilarity(i, i) = 0.0;
     for(size_t j = i+1; j < n ; ++j){
       
       if(dissimilarity(i, j) < 0){
         double i_j_dissimilarity = normalized_distance(i, j, M, index, data_.get_tree().max_depth());
         
         dissimilarity(i,j) = i_j_dissimilarity;
         dissimilarity(j, i) = i_j_dissimilarity;
       }
     }
   }
   dissimilarity(n-1, n-1) = 0.0;
   
   return dissimilarity;
}

template class GeneticAlgorithm<int>;
template class GeneticAlgorithm<double>;