#include "genetic_algorithm.h"

template<typename TargetType>
GeneticAlgorithm<TargetType>::GeneticAlgorithm(const PatientData<TargetType>& data,
                                               const GAParams& params)
  : data_{data}, params_{params} {
    
    population_.reserve(params.population_size);
    std::random_device rd;
    std::vector<unsigned int> seed_data;
    for (int i = 0; i < 4; ++i) {
      seed_data.push_back(rd());
    }
    
    std::seed_seq seq(seed_data.begin(), seed_data.end());
    
    rng_ = std::mt19937(seq);
}

template<typename TargetType>
void GeneticAlgorithm<TargetType>::initialize(){
  std::uniform_int_distribution<size_t> node(0, data_.get_tree().get_depth().size()-1);
  std::poisson_distribution<int> solution_size(data_.mean_patient_nodes());
  
  for(size_t i = 0; i < params_.population_size; ++i){
    population_.emplace_back(Solution::create_random_valid(data_.get_tree(), 
                                                 rng_,
                                                 solution_size(rng_)));
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
  
  default :
    Rcpp::stop("Unknown score type");
  }
}

//try to put a caching system. Queue might become too large, set a limit?
template<typename TargetType>
void GeneticAlgorithm<TargetType>::evaluate(){
  
  for(Solution& solution : population_){
    auto search = score_hash_map.find(solution.get_hash());
    if(search != score_hash_map.end()){
      solution.set_score(search->second);
    }
    else{
      solution.set_score(compute_score(solution));
      score_hash_map.insert({solution.get_hash(), solution.get_score()});
    }
  }
}

template<typename TargetType>
void GeneticAlgorithm<TargetType>::penalize(){
  auto [M, index] = data_.dissimilarity_input_computation();
  std::vector<std::vector<double>> 
    similarity(population_.size(), std::vector<double>(population_.size(), -1.0));
  
  for(size_t i = 0; i < index.size() -1 ; ++i){
    similarity[i][i] = 1;
    for(size_t j = i+1; j < index.size(); ++j){
      
      if(similarity[i][j] < 0){
        double i_j_similarity = 1.0 -
          normalized_distance(i, j, M, index, data_.get_tree().max_depth());
        
        similarity[i][j] = i_j_similarity;
        similarity[j][i] = i_j_similarity;
      }
    }
  }
  
  
}

template<typename TargetType>
void GeneticAlgorithm<TargetType>::keep_elite(std::vector<Solution>& offspring) const{
  std::priority_queue<Solution> elite_queue(population_.begin(),
                                            population_.end());
  
  for(int i = 0; i < params_.elite_count; ++i){
    offspring.emplace_back(elite_queue.top());
    elite_queue.pop();
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
  
  for(size_t i = params_.elite_count; i < offspring.size(); ++i){
    double mutation_probability = R::runif(0.0,1.0);
    
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
  
  for(size_t i = 0; i < offspring.size() - 2; i+=2){
    auto [c1,c2] = Solution::crossover_single_point(offspring[i],
                                                    offspring[i+1],
                                                    data_.get_tree(),
                                                    rng_);
    offspring[i] = c1;
    offspring[i+1] = c2;
  }
  
  //for even size population, we have to do one more crossover, if not we 
  //do no crossover on the last solution.
  if(offspring.size() % 2 == 0){
    auto [c1,c2] = Solution::crossover_single_point(offspring[offspring.size()-2],
                                                    offspring[offspring.size()-1],
                                                             data_.get_tree(),
                                                             rng_);
    offspring[offspring.size()-2] = c1;
    offspring[offspring.size()-1] = c2;
  }
}

template<typename TargetType>
void GeneticAlgorithm<TargetType>::replace_population(std::vector<Solution>& offspring){
  population_ = std::move(offspring);
  offspring.clear();
  offspring.reserve(population_.size());
}

template<typename TargetType>
void GeneticAlgorithm<TargetType>::run(){
  initialize();
  
  std::vector<Solution> offspring;
  offspring.reserve(population_.size());
  
  for(size_t i = 0 ; i < params_.epochs; ++i){
    evaluate();
    if(params_.diversity)
      //diversity here
    
    if(params_.elite_count > 0){
      keep_elite(offspring);
    }
    
    tournament_selection(offspring);
    mutate(offspring);
    crossover(offspring);
    
    replace_population(offspring);
  }
}

template class GeneticAlgorithm<int>;
template class GeneticAlgorithm<double>;