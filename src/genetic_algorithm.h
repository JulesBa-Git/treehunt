#ifndef GENETIC_ALGORITHM_H
#define GENETIC_ALGORITHM_H

#include "score_functions.h"
#include "solution.h"
#include <random>


struct GAParams {
  size_t population_size;
  size_t epochs;   
  double mutation_rate;
  double prob_mutation_type1;
  double crossover_rate;  
  size_t elite_count;
  size_t tournament_size;
  double alpha;   
  ScoreType score_type;     
  bool penalize;
  bool verbose;             
  
  GAParams() : population_size(100), epochs(1000), mutation_rate(0.1),
  prob_mutation_type1(0.5), crossover_rate(0.8), elite_count(10),
  tournament_size(3), alpha(1.0), score_type(ScoreType::HYPERGEOMETRIC), 
  penalize(false), verbose(false) {}
};

template <typename TargetType>
class GeneticAlgorithm{
private:
  std::vector<Solution> population_;
  std::unordered_map<size_t, double> score_hash_map;
  const PatientData<TargetType>& data_;
  GAParams params_;
  std::mt19937 rng_;
   
public :
  GeneticAlgorithm() = delete;
  GeneticAlgorithm(const PatientData<TargetType>& data,
                   const GAParams& params);
  
  void initialize();
  
  double compute_score(const Solution& sol) const;
  
  void evaluate();
  
  
};

#endif