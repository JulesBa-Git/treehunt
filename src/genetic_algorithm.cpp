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

template class GeneticAlgorithm<int>;
template class GeneticAlgorithm<double>;