#ifndef GENETIC_ALGORITHM_H
#define GENETIC_ALGORITHM_H

#include "helpers.h"
#include "score_functions.h"
#include "solution.h"
#include <random>
#include <memory>


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
  bool diversity;
  bool verbose;             
  
  GAParams() : population_size(100), epochs(1000), mutation_rate(0.1),
  prob_mutation_type1(0.5), crossover_rate(0.8), elite_count(10),
  tournament_size(3), alpha(1.0), score_type(ScoreType::HYPERGEOMETRIC), 
  diversity(false), verbose(false) {}
};

struct GAResults {
  std::vector<std::vector<int>> final_population;
  std::vector<double> final_scores;
  size_t total_generations;
  size_t cache_hits;
  
  GAResults() : total_generations(0) {}
};

template <typename TargetType>
class GeneticAlgorithm{
private:
  std::vector<Solution> population_;
  std::unordered_map<size_t, double> score_hash_map_;
  const PatientData<TargetType>& data_;
  GAParams params_;
  mutable std::mt19937 rng_; //concurrency not really possible here
  //tracking of cache_hits
  size_t cache_hits_;
   
public :
  GeneticAlgorithm() = delete;
  GeneticAlgorithm(const PatientData<TargetType>& data,
                   const GAParams& params);
  
  void initialize();
  
  double compute_score(const Solution& sol) const;
  void evaluate();
  void penalize();
  void keep_elite(std::vector<Solution>& offspring) const;
  void tournament_selection(std::vector<Solution>& offspring) const;
  void mutate(std::vector<Solution>& offspring) const;
  void crossover(std::vector<Solution>& offspring) const;
  void replace_population(std::vector<Solution>& offspring);
  GAResults run();
  
  GAResults extract_results() const;
  const std::vector<Solution>& get_population() const { return population_; }
  
};

#endif