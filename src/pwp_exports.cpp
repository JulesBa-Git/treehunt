#include <Rcpp.h>

#include "genetic_algorithm.h"
#include "mcmc_algorithm.h"
#include "patient_data.h"
#include "pwp_score.h"
#include "solution.h"
#include "tree_structure.h"

namespace {

std::vector<int> convert_nodes(const Rcpp::IntegerVector& nodes,
                               int index_base) {
  if (index_base != 0 && index_base != 1) {
    Rcpp::stop("index_base must be 0 or 1.");
  }
  std::vector<int> out = Rcpp::as<std::vector<int>>(nodes);
  if (index_base == 1) {
    std::transform(out.begin(), out.end(), out.begin(),
                   [](int x) { return x - 1; });
  }
  return out;
}

Rcpp::List format_ga_results(const GAResults& results,
                             std::size_t population_size,
                             std::size_t epochs) {
  Rcpp::List population(results.final_population.size());
  for (std::size_t i = 0; i < results.final_population.size(); ++i) {
    population[i] = Rcpp::wrap(results.final_population[i]);
  }
  return Rcpp::List::create(
      Rcpp::Named("final_population") = population,
      Rcpp::Named("final_scores") = Rcpp::wrap(results.final_scores),
      Rcpp::Named("parameters") = Rcpp::List::create(
          Rcpp::Named("population_size") = population_size,
          Rcpp::Named("epochs") = epochs,
          Rcpp::Named("score_type") = "pwp_rao"),
      Rcpp::Named("statistics") = Rcpp::List::create(
          Rcpp::Named("total_generations") = results.total_generations,
          Rcpp::Named("cache_hits") = results.cache_hits));
}

Rcpp::List format_mcmc_results(const MCMCResults& results) {
  Rcpp::List top(results.top_solutions.size());
  for (std::size_t i = 0; i < results.top_solutions.size(); ++i) {
    top[i] = Rcpp::wrap(results.top_solutions[i]);
  }
  Rcpp::List filtered(results.top_solutions_filtered.size());
  for (std::size_t i = 0; i < results.top_solutions_filtered.size(); ++i) {
    filtered[i] = Rcpp::wrap(results.top_solutions_filtered[i]);
  }
  const double acceptance_rate = results.total_iterations == 0 ? 0.0 :
      static_cast<double>(results.accepted_moves) /
      static_cast<double>(results.total_iterations);
  return Rcpp::List::create(
      Rcpp::Named("top_solutions") = top,
      Rcpp::Named("top_scores") = Rcpp::wrap(results.top_scores),
      Rcpp::Named("top_solutions_filtered") = filtered,
      Rcpp::Named("top_scores_filtered") =
          Rcpp::wrap(results.top_scores_filtered),
      Rcpp::Named("score_distribution") =
          Rcpp::wrap(results.score_distribution),
      Rcpp::Named("score_distribution_filtered") =
          Rcpp::wrap(results.score_distribution_filtered),
      Rcpp::Named("outstanding_scores") =
          Rcpp::wrap(results.outstanding_scores),
      Rcpp::Named("statistics") = Rcpp::List::create(
          Rcpp::Named("total_iterations") = results.total_iterations,
          Rcpp::Named("accepted_moves") = results.accepted_moves,
          Rcpp::Named("rejected_moves") = results.rejected_moves,
          Rcpp::Named("acceptance_rate") = acceptance_rate,
          Rcpp::Named("proposals_not_in_population") =
              results.proposals_not_in_population,
          Rcpp::Named("type1_moves") = results.type1_moves,
          Rcpp::Named("type2_moves") = results.type2_moves,
          Rcpp::Named("type1_accepted") = results.type1_accepted,
          Rcpp::Named("type2_accepted") = results.type2_accepted,
          Rcpp::Named("cocktail_size") = results.cocktail_size));
}

} // namespace

//' Internal C++ engine for the PWP Rao genetic algorithm
//'
//' Use \code{run_pwp_genetic_algorithm()} from R.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List run_genetic_algorithm_pwp_cpp(
    Rcpp::DataFrame patient_data,
    SEXP node_column,
    SEXP target_column,
    SEXP id_column,
    Rcpp::DataFrame tree,
    SEXP depth_column,
    SEXP upper_bound_column,
    SEXP name_column,
    Rcpp::List score_context,
    Rcpp::Nullable<Rcpp::List> seed_population = R_NilValue,
    std::size_t population_size = 100,
    std::size_t epochs = 1000,
    double mutation_rate = 0.1,
    double prob_mutation_type1 = 0.2,
    double crossover_rate = 0.8,
    std::size_t elite_count = 2,
    std::size_t tournament_size = 3,
    double alpha = 1.0,
    bool diversity = false,
    int seed = 1,
    bool verbose = false) {
  if (patient_data.nrows() == 0 || tree.nrows() == 0) {
    Rcpp::stop("patient_data and tree must be non-empty.");
  }
  if (population_size < 2 || epochs == 0 || elite_count >= population_size ||
      tournament_size < 1 || tournament_size > population_size) {
    Rcpp::stop("Invalid genetic-algorithm size parameter.");
  }
  if (mutation_rate < 0.0 || mutation_rate > 1.0 ||
      prob_mutation_type1 < 0.0 || prob_mutation_type1 > 1.0 ||
      crossover_rate < 0.0 || crossover_rate > 1.0) {
    Rcpp::stop("Genetic-algorithm probabilities must lie in [0, 1].");
  }

  tree_structure cpp_tree(tree, depth_column, upper_bound_column, name_column);
  PatientData<int> data(patient_data, node_column, target_column, cpp_tree,
                        id_column);
  PWPScoreContext pwp_context(data, score_context);

  GAParams params;
  params.population_size = population_size;
  params.epochs = epochs;
  params.mutation_rate = mutation_rate;
  params.prob_mutation_type1 = prob_mutation_type1;
  params.crossover_rate = crossover_rate;
  params.elite_count = elite_count;
  params.tournament_size = tournament_size;
  params.alpha = alpha;
  params.score_type = ScoreType::PWP_RAO;
  params.diversity = diversity;
  params.seed = seed;
  params.verbose = verbose;

  GeneticAlgorithm<int> algorithm(data, params, seed_population, &pwp_context);
  Rcpp::List out = format_ga_results(algorithm.run(), population_size, epochs);
  out["pwp_context_statistics"] = Rcpp::List::create(
      Rcpp::Named("event_time_groups") = pwp_context.number_event_groups(),
      Rcpp::Named("nuisance_covariates") =
          pwp_context.number_covariates());
  return out;
}

//' Internal C++ engine for PWP Rao MCMC
//'
//' Use \code{run_pwp_mcmc()} from R.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List run_mcmc_pwp_cpp(
    Rcpp::DataFrame patient_data,
    SEXP node_column,
    SEXP target_column,
    SEXP id_column,
    Rcpp::DataFrame tree,
    SEXP depth_column,
    SEXP upper_bound_column,
    SEXP name_column,
    Rcpp::List score_context,
    std::size_t epochs = 100000,
    double temperature = 1.0,
    std::size_t n_results = 20,
    std::size_t cocktail_size = 2,
    double prob_type1 = 0.01,
    std::size_t beta = 20,
    double max_score = 50.0,
    int seed = 1,
    bool verbose = false) {
  if (epochs == 0 || temperature <= 0.0 || prob_type1 < 0.0 ||
      prob_type1 > 1.0) {
    Rcpp::stop("Invalid MCMC parameter.");
  }
  tree_structure cpp_tree(tree, depth_column, upper_bound_column, name_column);
  PatientData<int> data(patient_data, node_column, target_column, cpp_tree,
                        id_column);
  PWPScoreContext pwp_context(data, score_context);

  MCMCParams params;
  params.epochs = epochs;
  params.temperature = temperature;
  params.n_results = n_results;
  params.cocktail_size = cocktail_size;
  params.prob_mutation_type1 = prob_type1;
  params.beta = beta;
  params.max_score = max_score;
  params.score_type_ = ScoreType::PWP_RAO;
  params.seed = seed;
  params.verbose = verbose;

  MCMCAlgorithm<int> algorithm(data, params, &pwp_context);
  Rcpp::List out = format_mcmc_results(algorithm.run());
  out["pwp_context_statistics"] = Rcpp::List::create(
      Rcpp::Named("event_time_groups") = pwp_context.number_event_groups(),
      Rcpp::Named("nuisance_covariates") =
          pwp_context.number_covariates());
  return out;
}

//' Evaluate PWP Rao statistics for specified combinations
//'
//' @keywords internal
// [[Rcpp::export]]
Rcpp::DataFrame score_pwp_combinations_cpp(
    Rcpp::List combinations,
    int index_base,
    Rcpp::DataFrame patient_data,
    SEXP node_column,
    SEXP target_column,
    SEXP id_column,
    Rcpp::DataFrame tree,
    SEXP depth_column,
    SEXP upper_bound_column,
    SEXP name_column,
    Rcpp::List score_context) {
  tree_structure cpp_tree(tree, depth_column, upper_bound_column, name_column);
  PatientData<int> data(patient_data, node_column, target_column, cpp_tree,
                        id_column);
  PWPScoreContext pwp_context(data, score_context);

  const int n = combinations.size();
  Rcpp::NumericVector fitness(n), signed_z(n), statistic(n), p_value(n),
      score(n), information(n);
  Rcpp::LogicalVector eligible(n);
  Rcpp::IntegerVector covered_rows(n), covered_patients(n), covered_events(n);
  for (int i = 0; i < n; ++i) {
    const auto nodes = convert_nodes(Rcpp::as<Rcpp::IntegerVector>(
                                         combinations[i]), index_base);
    const auto details = pwp_context.compute(Solution(nodes));
    eligible[i] = details.eligible;
    fitness[i] = details.fitness;
    signed_z[i] = details.signed_z;
    statistic[i] = details.score_statistic;
    p_value[i] = details.p_value;
    score[i] = details.score;
    information[i] = details.efficient_information;
    covered_rows[i] = static_cast<int>(details.covered_rows);
    covered_patients[i] = static_cast<int>(details.covered_patients);
    covered_events[i] = static_cast<int>(details.covered_events);
  }
  return Rcpp::DataFrame::create(
      Rcpp::Named("eligible") = eligible,
      Rcpp::Named("fitness") = fitness,
      Rcpp::Named("signed_z") = signed_z,
      Rcpp::Named("rao_chisq") = statistic,
      Rcpp::Named("p_value_model_based") = p_value,
      Rcpp::Named("score_u") = score,
      Rcpp::Named("efficient_information") = information,
      Rcpp::Named("covered_intervals") = covered_rows,
      Rcpp::Named("covered_patients") = covered_patients,
      Rcpp::Named("covered_events") = covered_events);
}

//' Construct one hierarchical combination indicator
//'
//' @keywords internal
// [[Rcpp::export]]
Rcpp::LogicalVector pwp_combination_flag_cpp(
    Rcpp::IntegerVector combination,
    int index_base,
    Rcpp::DataFrame patient_data,
    SEXP node_column,
    SEXP target_column,
    SEXP id_column,
    Rcpp::DataFrame tree,
    SEXP depth_column,
    SEXP upper_bound_column,
    SEXP name_column) {
  tree_structure cpp_tree(tree, depth_column, upper_bound_column, name_column);
  PatientData<int> data(patient_data, node_column, target_column, cpp_tree,
                        id_column);
  const auto nodes = convert_nodes(combination, index_base);
  Rcpp::LogicalVector flag(data.size());
  for (std::size_t i = 0; i < data.size(); ++i) {
    flag[i] = data.patient_has_combination(i, nodes);
  }
  return flag;
}
