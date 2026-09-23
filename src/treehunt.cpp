#include <Rcpp.h>
#include "mcmc_algorithm.h"
#include "mcmc_output.h"
#include "genetic_algorithm.h"
#include "patient_data.h"
#include "tree_structure.h"
#include <cmath>
#include <limits>


// Helper: Detect target type from R vector

enum class TargetTypeDetected {
  BINARY,
  CONTINUOUS
};

TargetTypeDetected detect_target_type(const Rcpp::DataFrame& df, SEXP target_column) {
  SEXP col;
  if (TYPEOF(target_column) == STRSXP) {
    col = df[Rcpp::as<Rcpp::String>(target_column)];
  } else {
    col = df[Rcpp::as<int>(target_column) - 1];
  }
  
  if (TYPEOF(col) == INTSXP || TYPEOF(col) == LGLSXP) {
    return TargetTypeDetected::BINARY;
  } else if (TYPEOF(col) == REALSXP) {
    // Check if it's actually binary (0/1 doubles)
    Rcpp::NumericVector nums(col);
    bool is_binary = true;
    for (int i = 0; i < nums.size() && is_binary; ++i) {
      if (nums[i] != 0.0 && nums[i] != 1.0) {
        is_binary = false;
      }
    }
    return is_binary ? TargetTypeDetected::BINARY : TargetTypeDetected::CONTINUOUS;
  }
  
  return TargetTypeDetected::BINARY;
}

// Helper: Convert ScoreType string to enum

ScoreType parse_score_type(const std::string& score_type_str) {
  if (score_type_str == "hypergeometric" || score_type_str == "HYPERGEOMETRIC") {
    return ScoreType::HYPERGEOMETRIC;
  } else if (score_type_str == "relative_risk" || score_type_str == "RELATIVE_RISK" ||
    score_type_str == "rr" || score_type_str == "RR") {
    return ScoreType::RELATIVE_RISK;
  } else if (score_type_str == "Wilcoxon" || score_type_str == "WILCOXON" ||
    score_type_str == "wilcoxon" || score_type_str == "wilc") {
    return ScoreType::WILCOXON;
  }else if (score_type_str == "composite" || score_type_str == "Composite" ||
    score_type_str == "QT") {
    return ScoreType::COMPOSITE;
  } else if(score_type_str == "residuals" || score_type_str == "RESIDUALS" ||
    score_type_str == "LMM"){
    return ScoreType::RESIDUALS;
  } else if(score_type_str == "bootstrap" || score_type_str == "BOOTSTRAP" ||
    score_type_str == "Boot"){
    return ScoreType::BOOTSTRAP;
  } else if(score_type_str == "pwp_rao" || score_type_str == "PWP_RAO" ||
    score_type_str == "rao"){
    return ScoreType::PWP_RAO;
  }
  else {
    Rcpp::stop("Unknown score type: '" + score_type_str + 
      "'. Use one of 'hypergeometric', 'relative_risk', 'wilcoxon'.");
  }
}

// A NULL seed preserves the historical non-deterministic initialization.
// Otherwise require a scalar non-negative integer that can be stored in the
// engine parameter structs.
int parse_optional_seed(SEXP seed) {
  if (Rf_isNull(seed))
    return -1;

  if (Rf_length(seed) != 1 ||
      (TYPEOF(seed) != INTSXP && TYPEOF(seed) != REALSXP)) {
    Rcpp::stop("seed must be NULL or one non-negative integer.");
  }

  const double value = Rcpp::as<double>(seed);
  if (!std::isfinite(value) || value < 0.0 || std::floor(value) != value ||
      value > static_cast<double>(std::numeric_limits<int>::max())) {
    Rcpp::stop("seed must be NULL or one non-negative integer.");
  }

  return static_cast<int>(value);
}

// MCMC Algorithm Interface


//' Run MCMC Algorithm for Estimation of Score Distribution Among Nodes of The
//' Tree
//'
//' Performs a Modified Metropolis-Hastings MCMC sampling to estimate the score
//' distribution of nodes combination of a given \emph{cocktail_size}. The
//' algorithm explores the space of tree combinations using a proposal law 
//' composed of two mutation types.
//'
//' @param patient_data A data.frame containing patient information with at least
//'   a node column and a target column.
//' @param node_column Either a string (column name) or integer (column index, 1-based)
//'   specifying the column containing drug codes. This column should be either:
//'   \itemize{
//'     \item A list of integer vectors: \code{list(c(1,2), c(3), c(4,5))}
//'     \item A character vector with comma-separated values: \code{c("1,2", "3", "4,5")}
//'   }
//' @param id_column Optional observation-unit identifier column, given by name
//'   or one-based position. Required for Wilcoxon and residual patient-level scores.
//' @param target_column Either a string (column name) or integer (column index, 1-based)
//'   specifying the target/outcome column. Integer values are treated as binary for now,
//'   numeric values with non-0/1 entries are treated as continuous.
//' @param tree_depth An integer vector specifying the depth of each node in the
//'   tree structure. Must start at depth 1 and children must be at depth+1 of
//'   their parent.
//' @param epochs Number of MCMC iterations to run.
//' @param temperature Temperature parameter for the Metropolis-Hastings acceptance
//'   probability. Higher values lead to an easiest acceptance of lower score. Default: 1.0.
//' @param n_results Number of top solutions to track and return. Default: 10.
//' @param cocktail_size Target size of drug combinations to search for. Default: 2.
//' @param prob_type1 Probability of using Type 1 mutation (random generation) vs
//'   Type 2 mutation (local swap). Default: 0.01.
//' @param beta Strict support threshold: filtered results require more than
//'   beta covered observations (or distinct patients for patient-level scores). Default: 4.
//' @param max_score Finite positive cap applied to scores before Metropolis-Hastings
//'   acceptance, importance weighting and binning. Default: 200.0.
//' @param score_type Scoring function to use. Either "hypergeometric" for the
//'   hypergeometric test, "relative_risk" for relative risk calculation, or "wilcoxon"
//'   for the wilcoxon test with continuous output.
//'   Default: "hypergeometric".
//' @param verbose If TRUE, prints progress and statistics during the run.
//'   Default: FALSE.
//' @param seed Optional non-negative integer seed for the C++ random-number
//'   generator. If \code{NULL}, the generator is initialized non-deterministically.
//'
//' @param burn_in Number of initial iterations excluded from all distributions,
//'   top solutions and the optional trace. Must be smaller than epochs. Default: 0.
//' @param store_trace If TRUE, return the retained capped score, support and
//'   zero-based node vector for each iteration. Default FALSE avoids this storage.
//'
//' @return A list containing:
//'   \describe{
//'     \item{top_solutions}{List of zero-based node vectors for the top scoring solutions}
//'     \item{top_scores}{Numeric vector of scores for the top solutions}
//'     \item{top_solutions_filtered}{Top solutions meeting the beta threshold}
//'     \item{top_scores_filtered}{Scores for the filtered solutions}
//'     \item{score_distribution}{Histogram of retained capped scores (0.1-wide bins and a cap bin)}
//'     \item{score_distribution_filtered}{Histogram for solutions meeting beta threshold}
//'     \item{outstanding_scores}{Retained capped scores equal to max_score}
//'     \item{uniform_reference}{Uniform-reference distribution and weight ESS}
//'     \item{uniform_reference_filtered}{Support-filtered uniform reference and ESS}
//'     \item{reference_parameters}{Temperature, score cap, support and burn-in settings}
//'     \item{trace}{Retained trajectory, only when store_trace is TRUE}
//'     \item{statistics}{List of run statistics including acceptance rates}
//'   }
//'
//' @details
//' The stationary target on observed subsets of \code{cocktail_size} distinct
//' nodes is proportional to \eqn{\exp(S(C)/T)}, where \eqn{S} is capped at
//' \code{max_score}. Ancestor-descendant pairs are included in this MCMC space;
//' genetic-algorithm validity filtering is a separate operation.
//'
//' Type 1 proposes a uniform subset of distinct nodes. Type 2 replaces a node
//' by a parent or child not already in the combination. Its Hastings correction
//' is \eqn{|V(C)|/|V(C')|}, counting these possible local moves before rejecting
//' proposals absent from the data. Each acceptance probability is the minimum
//' of 1 and \eqn{\exp((S(C')-S(C))/T)} times the applicable proposal correction.
//' A positive \code{prob_type1} connects the observed state space.
//'
//' At every retained iteration, including rejected moves, weights proportional
//' to \eqn{\exp(-S(C)/T)} recover the uniform combination reference. Weights are
//' computed from exact capped scores before grouping into bins. They are
//' rescaled internally to avoid numerical underflow. The filtered reference
//' conditions on support strictly greater than \code{beta}, with its own ESS.
//' This estimates a uniform law on combinations, not on distinct score values.
//'
//' \code{uniform_reference} and \code{uniform_reference_filtered} each contain
//' \code{distribution}, \code{n_samples}, and \code{weight_ess}. The distribution
//' columns are \code{lower}, \code{upper}, \code{probability}, \code{cdf} and
//' \code{upper_tail}. Bins are [lower, upper), except the final singleton bin
//' at the cap. CDF includes the current bin; upper_tail includes it and all
//' subsequent bins. Empty references have NA probabilities and ESS.
//' ESS is \eqn{(\sum w)^2/\sum w^2}; it measures weight concentration and does
//' not account for MCMC autocorrelation. Check mixing separately and choose
//' an appropriate \code{burn_in}. See \code{\link{uniform_score_reference}}
//' for post-processing an optionally stored trace without binning.
//'
//' @examples
//' dat <- data.frame(outcome = c(1L, 1L, 0L, 0L))
//' dat$nodes <- list(c(0L, 1L), 0L, 1L, 1L)
//' tree <- data.frame(Depth = c(1L, 1L))
//' result <- run_mcmc(dat, "nodes", "outcome", c(1L, 1L),
//'   epochs = 100, cocktail_size = 1, seed = 42)
//' result$uniform_reference$weight_ess
//'
//' @export
// [[Rcpp::export]]
Rcpp::List run_mcmc(
   Rcpp::DataFrame patient_data,
   SEXP node_column,
   SEXP target_column,
   Rcpp::IntegerVector tree_depth,
   double epochs,
   double temperature = 1.0,
   double n_results = 10,
   double cocktail_size = 2,
   double prob_type1 = 0.01,
   double beta = 4,
   double max_score = 200.0,
   std::string score_type = "hypergeometric",
   bool verbose = false,
   SEXP seed = R_NilValue,
   double burn_in = 0,
   bool store_trace = false,
   SEXP id_column = R_NilValue) {
 
 // Validate inputs
 if (patient_data.nrows() == 0) {
   Rcpp::stop("patient_data cannot be empty");
 }
 if (tree_depth.size() == 0) {
   Rcpp::stop("tree_depth cannot be empty");
 }
 if (epochs == 0) {
   Rcpp::stop("epochs must be positive");
 }
 if (temperature <= 0) {
   Rcpp::stop("temperature must be positive");
 }
 if (prob_type1 < 0 || prob_type1 > 1) {
   Rcpp::stop("prob_type1 must be between 0 and 1");
 }
 
 // Build tree structure
 tree_structure tree(tree_depth);
 
 // Setup MCMC parameters
 MCMCParams params;
 params.epochs = mcmc_count(epochs, "epochs");
 params.temperature = temperature;
 params.n_results = mcmc_count(n_results, "n_results");
 params.cocktail_size = mcmc_count(cocktail_size, "cocktail_size");
 params.prob_mutation_type1 = prob_type1;
 params.beta = mcmc_count(beta, "beta");
 params.max_score = max_score;
 params.score_type_ = parse_score_type(score_type);
 params.verbose = verbose;
 params.seed = parse_optional_seed(seed);
 params.burn_in = mcmc_count(burn_in, "burn_in");
 params.store_trace = store_trace;
 
 // Detect target type and run appropriate template
 TargetTypeDetected target_type = detect_target_type(patient_data, target_column);
 MCMCResults results;
 
 if (target_type == TargetTypeDetected::BINARY) {
   PatientData<int> data(patient_data, node_column, target_column, tree, id_column);
   MCMCAlgorithm<int> algorithm(data, params);
   results = algorithm.run();
 } else {
   PatientData<double> data(patient_data, node_column, target_column, tree, id_column);
   MCMCAlgorithm<double> algorithm(data, params);
   results = algorithm.run();
 }
 
 return mcmc_output(results, params);
}

// Genetic Algorithm Interface

//' Run Genetic Algorithm for High Score Nodes Combination Search
//'
//' Performs a genetic algorithm search to find optimal node combinations that
//' maximize a specified score function. The algorithm evolves a population of
//' solutions through selection, crossover, and mutation operations.
//'
//' @param patient_data A data.frame containing patient information with at least
//'   a node column and a target column.
//' @param node_column Either a string (column name) or integer (column index, 1-based)
//'   specifying the column containing node indexes. This column should be either:
//'   \itemize{
//'     \item A list of integer vectors: \code{list(c(1,2), c(3), c(4,5))} (0 indexed Note : Might be wise to change this)
//'     \item A character vector with comma-separated values: \code{c("1,2", "3", "4,5")}
//'   }
//' @param target_column Either a string (column name) or integer (column index, 1-based)
//'   specifying the target/outcome column. Integer values are treated as binary,
//'   numeric values with non-0/1 entries are treated as continuous.
//' @param tree_depth An integer vector specifying the depth of each node in the
//'   tree structure. Must start at depth 1 and children must be at depth+1 of
//'   their parent.
//' @param seed_population A \code{list} of integer vectors representing an initial 
//'   population (vector of 1-based tree index). If provided, these individuals 
//'   will be included in the first generation . If the list contains fewer 
//'   individuals than \code{population_size}, the remainder will be initialized 
//'   randomly. If \code{NULL} (the default), the entire initial population 
//'   is generated randomly.
//' @param population_size Number of solutions in the population. Default: 100.
//' @param epochs Number of generations to evolve. Default: 1000.
//' @param mutation_rate Probability of mutating each offspring. Default: 0.1.
//' @param prob_mutation_type1 When mutation occurs, probability of using Type 1
//'   (add/remove) vs Type 2 (swap) mutation. Default: 0.2.
//' @param crossover_rate Probability of applying crossover to selected parents.
//'   Default: 0.8.
//' @param elite_count Number of top solutions to preserve unchanged each generation.
//'   Default: 0.
//' @param tournament_size Number of solutions competing in tournament selection.
//'   Default: 2.
//' @param alpha Parameter controlling the add/remove mutation bias. Higher values
//'   favor adding nodes. Default: 1.0.
//' @param score_type Scoring function to use. Either "hypergeometric" for the
//'   hypergeometric test, "relative_risk" for relative risk calculation, or "wilcoxon"
//'   for the wilcoxon test with continuous output.
//' @param diversity If TRUE, applies a diversity penalty to encourage exploration
//'   of different solutions. Default: FALSE.
//' @param verbose If TRUE, prints progress during the run. Default: FALSE.
//' @param seed Optional non-negative integer seed for the C++ random-number
//'   generator. If \code{NULL}, the generator is initialized non-deterministically.
//'
//' @return A list containing:
//'   \describe{
//'     \item{final_population}{The complete final population of solutions}
//'     \item{final_scores}{Scores for all solutions in the final population}
//'     \item{parameters}{List of parameters used for the run}
//'     \item{statistics}{Additional information about cache hits etc.}
//'   }
//'
//' @details
//' The genetic algorithm uses the following operators:
//' \itemize{
//'   \item \strong{Selection}: Tournament selection with configurable size
//'   \item \strong{Crossover}: Single-point crossover on tree subtrees
//'   \item \strong{Mutation Type 1}: Add or remove a node (probability controlled by alpha)
//'   \item \strong{Mutation Type 2}: Swap a node with its parent or child
//'   \item \strong{Elitism}: Top solutions preserved unchanged
//' }
//'
//' The algorithm maintains a score cache to avoid redundant computations for
//' solutions that have been seen before (identified by hash).
//'
//' When \code{diversity = TRUE}, solutions are penalized based on their similarity
//' to other solutions in the population, encouraging exploration of diverse regions.
//'
//' @examples
//' \dontrun{
//' # Create example data
//' patient_df <- data.frame(
//'   patient_id = 1:100,
//'   outcome = rbinom(100, 1, 0.3)
//' )
//' patient_df$drugs <- lapply(1:100, function(i) sample(1:20, sample(1:5, 1)))
//'
//' # Define tree structure
//' tree_depth <- c(1, rep(2, 5), rep(3, 15))
//'
//' # Run GA
//' results <- run_genetic_algorithm(
//'   patient_data = patient_df,
//'   node_column = "drugs",
//'   target_column = "outcome",
//'   tree_depth = tree_depth,
//'   population_size = 50,
//'   epochs = 500,
//'   score_type = "hypergeometric",
//'   verbose = TRUE
//' )
//'
//' # View top results
//' print(results$top_scores)
//' print(results$top_solutions)
//' }
//'
//' @export
//' @seealso \code{\link{run_mcmc}} for an MCMC-based optimization approach
// [[Rcpp::export]]
Rcpp::List run_genetic_algorithm(
   Rcpp::DataFrame patient_data,
   SEXP node_column,
   SEXP target_column,
   Rcpp::IntegerVector tree_depth,
   Rcpp::Nullable<Rcpp::List> seed_population = R_NilValue,
   size_t population_size = 100,
   size_t epochs = 1000,
   double mutation_rate = 0.1,
   double prob_mutation_type1 = 0.2,
   double crossover_rate = 0.8,
   size_t elite_count = 0,
   size_t tournament_size = 3,
   double alpha = 1.0,
   std::string score_type = "hypergeometric",
   bool diversity = false,
   bool verbose = false,
   SEXP seed = R_NilValue) {
 
 // Validate inputs
 if (patient_data.nrows() == 0) {
   Rcpp::stop("patient_data cannot be empty");
 }
 if (tree_depth.size() == 0) {
   Rcpp::stop("tree_depth cannot be empty");
 }
 if (population_size < 2) {
   Rcpp::stop("population_size must be at least 2");
 }
 if (epochs <= 0) {
   Rcpp::stop("epochs must be positive");
 }
 if (mutation_rate < 0 || mutation_rate > 1) {
   Rcpp::stop("mutation_rate must be between 0 and 1");
 }
 if (!std::isfinite(prob_mutation_type1) || prob_mutation_type1 < 0 ||
     prob_mutation_type1 > 1) {
   Rcpp::stop("prob_mutation_type1 must be finite and between 0 and 1");
 }
 if (crossover_rate < 0 || crossover_rate > 1) {
   Rcpp::stop("crossover_rate must be between 0 and 1");
 }
 if (elite_count >= population_size || elite_count < 0 ) {
   Rcpp::stop("elite_count must be positive and less than population_size");
 }
 if (tournament_size < 1 || tournament_size > population_size) {
   Rcpp::stop("tournament_size must be between 1 and population_size");
 }
 
 // Build tree structure
 tree_structure tree(tree_depth);
 
 // Setup GA parameters
 GAParams params;
 params.population_size = population_size;
 params.epochs = epochs;
 params.mutation_rate = mutation_rate;
 params.prob_mutation_type1 = prob_mutation_type1;
 params.crossover_rate = crossover_rate;
 params.elite_count = elite_count;
 params.tournament_size = tournament_size;
 params.alpha = alpha;
 params.score_type = parse_score_type(score_type);
 params.diversity = diversity;
 params.verbose = verbose;
 params.seed = parse_optional_seed(seed);
 
 // Detect target type and run appropriate template
 TargetTypeDetected target_type = detect_target_type(patient_data, target_column);
 GAResults results;
 
 if (target_type == TargetTypeDetected::BINARY) {
   PatientData<int> data(patient_data, node_column, target_column, tree);
   GeneticAlgorithm<int> algorithm(data, params, seed_population);
   results = algorithm.run();
 } else {
   PatientData<double> data(patient_data, node_column, target_column, tree);
   GeneticAlgorithm<double> algorithm(data, params, seed_population);
   results = algorithm.run();
 }
 
 // Convert final_population to R list
 Rcpp::List final_population_list(results.final_population.size());
 for (size_t i = 0; i < results.final_population.size(); ++i) {
   final_population_list[i] = Rcpp::wrap(results.final_population[i]);
 }
 
 Rcpp::List parameters = Rcpp::List::create(
   Rcpp::Named("population_size") = population_size,
   Rcpp::Named("epochs") = epochs,
   Rcpp::Named("mutation_rate") = mutation_rate,
   Rcpp::Named("prob_mutation_type1") = prob_mutation_type1,
   Rcpp::Named("crossover_rate") = crossover_rate,
   Rcpp::Named("elite_count") = elite_count,
   Rcpp::Named("tournament_size") = tournament_size,
   Rcpp::Named("alpha") = alpha,
   Rcpp::Named("score_type") = score_type,
   Rcpp::Named("diversity") = diversity,
   Rcpp::Named("seed") =
     (params.seed >= 0 ? Rcpp::wrap(params.seed) : R_NilValue)
 );
 
 Rcpp::List statistics = Rcpp::List::create(
   Rcpp::Named("total_generations") = results.total_generations,
   Rcpp::Named("cache_hits") = results.cache_hits
 );
 
 return Rcpp::List::create(
   Rcpp::Named("final_population") = final_population_list,
   Rcpp::Named("final_scores") = Rcpp::wrap(results.final_scores),
   Rcpp::Named("parameters") = parameters,
   Rcpp::Named("statistics") = statistics
 );
}

//' Run MCMC Algorithm for Estimation of Score Distribution Among Nodes of The
//' Tree
//'
//' Performs a Modified Metropolis-Hastings MCMC sampling to estimate the score
//' distribution of nodes combination of a given \emph{cocktail_size}. The
//' algorithm explores the space of tree combinations using a proposal law 
//' composed of two mutation types.
//'
//' @param patient_data A data.frame containing patient information with at least
//'   a node column and a target column.
//' @param node_column Either a string (column name) or integer (column index, 1-based)
//'   specifying the column containing drug codes. This column should be either:
//'   \itemize{
//'     \item A list of integer vectors: \code{list(c(1,2), c(3), c(4,5))}
//'     \item A character vector with comma-separated values: \code{c("1,2", "3", "4,5")}
//'   }
//' @param id_column Optional observation-unit identifier column, given by name
//'   or one-based column position; needed for patient-level scoring.
//' @param target_column Either a string (column name) or integer (column index, 1-based)
//'   specifying the target/outcome column. Integer values are treated as binary for now,
//'   numeric values with non-0/1 entries are treated as continuous.
//' @param tree A data.frame containing the structural definition of the tree.
//' @param depth_column Either a string or integer specifying the column in 
//'   \code{tree} that contains the node depth levels.
//' @param upper_bound_column (Optional) Either a string or integer (1-based index) 
//'  specifying the column in \code{tree_depth} that contains upper 
//'  bound of nodes. Defaults to \code{NULL}.
//' @param name_column (Optional) Either a string or integer (1-based index) 
//'  specifying the column in \code{tree} that contains the corresponding name
//'  of nodes. Defaults to \code{NULL}.
//' @param epochs Number of MCMC iterations to run.
//' @param temperature Temperature parameter for the Metropolis-Hastings acceptance
//'   probability. Higher values lead to an easiest acceptance of lower score. Default: 1.0.
//' @param n_results Number of top solutions to track and return. Default: 10.
//' @param cocktail_size Target size of drug combinations to search for. Default: 2.
//' @param prob_type1 Probability of using Type 1 mutation (random generation) vs
//'   Type 2 mutation (local swap). Default: 0.01.
//' @param beta Strict support threshold: filtered results require more than
//'   beta covered observations (or distinct patients for patient-level scores). Default: 4.
//' @param max_score Finite positive cap applied to scores before Metropolis-Hastings
//'   acceptance, importance weighting and binning. Default: 200.0.
//' @param score_type Scoring function to use. Either "hypergeometric" for the
//'   hypergeometric test, "relative_risk" for relative risk calculation, or "wilcoxon"
//'   for the wilcoxon test with continuous output.
//'   Default: "hypergeometric".
//' @param verbose If TRUE, prints progress and statistics during the run.
//'   Default: FALSE.
//' @param seed Optional non-negative integer seed for the C++ random-number
//'   generator. If \code{NULL}, the generator is initialized non-deterministically.
//'
//' @param burn_in Number of initial iterations excluded from all distributions,
//'   top solutions and the optional trace. Must be smaller than epochs. Default: 0.
//' @param store_trace If TRUE, return the retained capped score, support and
//'   zero-based node vector for each iteration. Default FALSE avoids this storage.
//'
//' @return A list containing:
//'   \describe{
//'     \item{top_solutions}{List of zero-based node vectors for the top scoring solutions}
//'     \item{top_scores}{Numeric vector of scores for the top solutions}
//'     \item{top_solutions_filtered}{Top solutions meeting the beta threshold}
//'     \item{top_scores_filtered}{Scores for the filtered solutions}
//'     \item{score_distribution}{Histogram of retained capped scores (0.1-wide bins and a cap bin)}
//'     \item{score_distribution_filtered}{Histogram for solutions meeting beta threshold}
//'     \item{outstanding_scores}{Retained capped scores equal to max_score}
//'     \item{uniform_reference}{Uniform-reference distribution and weight ESS}
//'     \item{uniform_reference_filtered}{Support-filtered uniform reference and ESS}
//'     \item{reference_parameters}{Temperature, score cap, support and burn-in settings}
//'     \item{trace}{Retained trajectory, only when store_trace is TRUE}
//'     \item{statistics}{List of run statistics including acceptance rates}
//'   }
//'
//' @details
//' The stationary target on observed subsets of \code{cocktail_size} distinct
//' nodes is proportional to \eqn{\exp(S(C)/T)}, where \eqn{S} is capped at
//' \code{max_score}. Ancestor-descendant pairs are included in this MCMC space;
//' genetic-algorithm validity filtering is a separate operation.
//'
//' Type 1 proposes a uniform subset of distinct nodes. Type 2 replaces a node
//' by a parent or child not already in the combination. Its Hastings correction
//' is \eqn{|V(C)|/|V(C')|}, counting these possible local moves before rejecting
//' proposals absent from the data. Each acceptance probability is the minimum
//' of 1 and \eqn{\exp((S(C')-S(C))/T)} times the applicable proposal correction.
//' A positive \code{prob_type1} connects the observed state space.
//'
//' At every retained iteration, including rejected moves, weights proportional
//' to \eqn{\exp(-S(C)/T)} recover the uniform combination reference. Weights are
//' computed from exact capped scores before grouping into bins. They are
//' rescaled internally to avoid numerical underflow. The filtered reference
//' conditions on support strictly greater than \code{beta}, with its own ESS.
//' This estimates a uniform law on combinations, not on distinct score values.
//'
//' \code{uniform_reference} and \code{uniform_reference_filtered} each contain
//' \code{distribution}, \code{n_samples}, and \code{weight_ess}. The distribution
//' columns are \code{lower}, \code{upper}, \code{probability}, \code{cdf} and
//' \code{upper_tail}. Bins are [lower, upper), except the final singleton bin
//' at the cap. CDF includes the current bin; upper_tail includes it and all
//' subsequent bins. Empty references have NA probabilities and ESS.
//' ESS is \eqn{(\sum w)^2/\sum w^2}; it measures weight concentration and does
//' not account for MCMC autocorrelation. Check mixing separately and choose
//' an appropriate \code{burn_in}. See \code{\link{uniform_score_reference}}
//' for post-processing an optionally stored trace without binning.
//'
//' @examples
//' dat <- data.frame(outcome = c(1L, 1L, 0L, 0L))
//' dat$nodes <- list(c(0L, 1L), 0L, 1L, 1L)
//' tree <- data.frame(Depth = c(1L, 1L))
//' result <- run_mcmc_df_tree(dat, "nodes", "outcome", tree, "Depth",
//'   epochs = 100, cocktail_size = 1, seed = 42)
//' result$uniform_reference$weight_ess
//'
//' @export
// [[Rcpp::export]]
Rcpp::List run_mcmc_df_tree(
   Rcpp::DataFrame patient_data,
   SEXP node_column,
   SEXP target_column,
   Rcpp::DataFrame tree,
   SEXP depth_column,
   SEXP upper_bound_column = R_NilValue,
   SEXP name_column = R_NilValue,
   SEXP id_column = R_NilValue,
   double epochs = 1e6,
   double temperature = 1.0,
   double n_results = 10,
   double cocktail_size = 2,
   double prob_type1 = 0.01,
   double beta = 4,
   double max_score = 200.0,
   std::string score_type = "hypergeometric",
   bool verbose = false,
   SEXP seed = R_NilValue,
   double burn_in = 0,
   bool store_trace = false) {
 
 // Validate inputs
 if (patient_data.nrows() == 0) {
   Rcpp::stop("patient_data cannot be empty");
 }
 if (tree.nrow() == 0) {
   Rcpp::stop("tree_depth cannot be empty");
 }
 if (epochs == 0) {
   Rcpp::stop("epochs must be positive");
 }
 if (temperature <= 0) {
   Rcpp::stop("temperature must be positive");
 }
 if (prob_type1 < 0 || prob_type1 > 1) {
   Rcpp::stop("prob_type1 must be between 0 and 1");
 }
 
 // Build tree structure
 tree_structure cppTree(tree, depth_column, upper_bound_column, name_column);
 
 // Setup MCMC parameters
 MCMCParams params;
 params.epochs = mcmc_count(epochs, "epochs");
 params.temperature = temperature;
 params.n_results = mcmc_count(n_results, "n_results");
 params.cocktail_size = mcmc_count(cocktail_size, "cocktail_size");
 params.prob_mutation_type1 = prob_type1;
 params.beta = mcmc_count(beta, "beta");
 params.max_score = max_score;
 params.score_type_ = parse_score_type(score_type);
 params.verbose = verbose;
 params.seed = parse_optional_seed(seed);
 params.burn_in = mcmc_count(burn_in, "burn_in");
 params.store_trace = store_trace;
 
 // Detect target type and run appropriate template
 TargetTypeDetected target_type = detect_target_type(patient_data, target_column);
 MCMCResults results;
 
 if (target_type == TargetTypeDetected::BINARY) {
   PatientData<int> data(patient_data, node_column, target_column, cppTree,
                         id_column);
   MCMCAlgorithm<int> algorithm(data, params);
   results = algorithm.run();
 } else {
   PatientData<double> data(patient_data, node_column, target_column, cppTree,
                            id_column);
   MCMCAlgorithm<double> algorithm(data, params);
   results = algorithm.run();
 }
 
 return mcmc_output(results, params);
}



//' Run Genetic Algorithm for High Score Nodes Combination Search
//'
//' Performs a genetic algorithm search to find optimal node combinations that
//' maximize a specified score function. This version builds the internal tree 
//' structure from a provided data frame mapping node depths and bounds.
//'
//' @param patient_data A data.frame containing patient information.
//' @param node_column Either a string (column name) or integer (1-based index)
//'   specifying the column containing node indexes (list of vectors or comma-separated strings).
//' @param id_column Optional observation-unit identifier column, given by name
//'   or one-based column position; needed for patient-level scoring.
//' @param target_column Either a string (column name) or integer (1-based index)
//'   specifying the target/outcome column.
//' @param tree A data.frame containing the structural definition of the tree.
//' @param depth_column Either a string or integer specifying the column in 
//'   \code{tree} that contains the node depth levels.
//' @param upper_bound_column (Optional) Either a string or integer (1-based index) 
//'  specifying the column in \code{tree_depth} that contains upper 
//'  bound of nodes. Defaults to \code{NULL}.
//' @param name_column (Optional) Either a string or integer (1-based index) 
//'  specifying the column in \code{tree} that contains the corresponding name
//'  of nodes. Defaults to \code{NULL}.
//' @param seed_population A \code{list} of integer vectors representing an initial 
//'   population (vector of 1-based tree index). If provided, these individuals 
//'   will be included in the first generation . If the list contains fewer 
//'   individuals than \code{population_size}, the remainder will be initialized 
//'   randomly. If \code{NULL} (the default), the entire initial population 
//'   is generated randomly.
//' @param population_size Number of solutions in the population. Default: 100.
//' @param epochs Number of generations to evolve. Default: 1000.
//' @param mutation_rate Probability of mutating each offspring. Default: 0.1.
//' @param prob_mutation_type1 Probability of using Type 1 (add/remove) vs Type 2 (swap) mutation. Default: 0.2.
//' @param crossover_rate Probability of applying crossover to selected parents. Default: 0.8.
//' @param elite_count Number of top solutions to preserve unchanged each generation. Default: 0.
//' @param tournament_size Number of solutions competing in tournament selection. Default: 3.
//' @param alpha Parameter controlling add/remove mutation bias. Higher values favor adding nodes. Default: 1.0.
//' @param score_type Scoring function: "hypergeometric", "relative_risk", or "wilcoxon".
//' @param diversity If TRUE, applies a diversity penalty to encourage exploration. Default: FALSE.
//' @param verbose If TRUE, prints progress during the run. Default: FALSE.
//' @param seed Optional non-negative integer seed for the C++ random-number
//'   generator. If \code{NULL}, the generator is initialized non-deterministically.
//'
//' @return A list containing:
//'   \describe{
//'     \item{final_population}{The complete final population of solutions, 0-based
//'     tree index. TODO -> return 1 based (R-like) index for more coherence}
//'     \item{final_scores}{Scores for all solutions in the final population}
//'     \item{parameters}{List of parameters used for the run}
//'     \item{statistics}{Additional information about cache hits, etc.}
//'   }
//'
//' @details
//' Unlike the vector-based version, this function extracts tree hierarchy from 
//' the \code{tree_depth} data frame. It uses \code{depth_column} and optionally 
//' \code{upper_bound_column} and \code{name_column} to define the tree.
//'
//' @examples
//' \dontrun{
//' # Define tree structure via data frame
//' tree_df <- data.frame(
//'   node_id = 1:21,
//'   depth_level = c(1, rep(2, 5), rep(3, 15)), # User may add, upper bound 
//'   # and name column
//' )
//'
//' results <- run_genetic_algorithm_df_tree(
//'   patient_data = patient_df,
//'   node_column = "drugs",
//'   target_column = "outcome",
//'   tree = tree_df,
//'   depth_column = "depth_level", # or 2
//'   population_size = 100,
//'   score_type = "hypergeometric"
//' )
//' }
//'
//' @export
// [[Rcpp::export]]
Rcpp::List run_genetic_algorithm_df_tree(
   Rcpp::DataFrame patient_data,
   SEXP node_column,
   SEXP target_column,
   Rcpp::DataFrame tree,
   SEXP depth_column,
   SEXP upper_bound_column = R_NilValue,
   SEXP name_column = R_NilValue,
   SEXP id_column = R_NilValue,
   Rcpp::Nullable<Rcpp::List> seed_population = R_NilValue,
   size_t population_size = 100,
   size_t epochs = 1000,
   double mutation_rate = 0.1,
   double prob_mutation_type1 = 0.2,
   double crossover_rate = 0.8,
   size_t elite_count = 0,
   size_t tournament_size = 3,
   double alpha = 1.0,
   std::string score_type = "hypergeometric",
   bool diversity = false,
   bool verbose = false,
   SEXP seed = R_NilValue) {
 
 // Validate inputs
 if (patient_data.nrows() == 0) {
   Rcpp::stop("patient_data cannot be empty");
 }
 if (tree.nrow() == 0) {
   Rcpp::stop("tree cannot be empty");
 }
 if (population_size < 2) {
   Rcpp::stop("population_size must be at least 2");
 }
 if (epochs <= 0) {
   Rcpp::stop("epochs must be positive");
 }
 if (mutation_rate < 0 || mutation_rate > 1) {
   Rcpp::stop("mutation_rate must be between 0 and 1");
 }
 if (!std::isfinite(prob_mutation_type1) || prob_mutation_type1 < 0 ||
     prob_mutation_type1 > 1) {
   Rcpp::stop("prob_mutation_type1 must be finite and between 0 and 1");
 }
 if (crossover_rate < 0 || crossover_rate > 1) {
   Rcpp::stop("crossover_rate must be between 0 and 1");
 }
 if (elite_count >= population_size || elite_count < 0 ) {
   Rcpp::stop("elite_count must be positive and less than population_size");
 }
 if (tournament_size < 1 || tournament_size > population_size) {
   Rcpp::stop("tournament_size must be between 1 and population_size");
 }
 
 // Build tree structure
 tree_structure cppTree(tree, depth_column, upper_bound_column, name_column);
 
 // Setup GA parameters
 GAParams params;
 params.population_size = population_size;
 params.epochs = epochs;
 params.mutation_rate = mutation_rate;
 params.prob_mutation_type1 = prob_mutation_type1;
 params.crossover_rate = crossover_rate;
 params.elite_count = elite_count;
 params.tournament_size = tournament_size;
 params.alpha = alpha;
 params.score_type = parse_score_type(score_type);
 params.diversity = diversity;
 params.verbose = verbose;
 params.seed = parse_optional_seed(seed);

 // Detect target type and run appropriate template
 TargetTypeDetected target_type = detect_target_type(patient_data, target_column);
 GAResults results;

 if (target_type == TargetTypeDetected::BINARY) {
   PatientData<int> data(patient_data, node_column, target_column, cppTree, 
                            id_column);
   GeneticAlgorithm<int> algorithm(data, params, seed_population);
   results = algorithm.run();
 } else {
   PatientData<double> data(patient_data, node_column, target_column, cppTree, 
                            id_column);
   GeneticAlgorithm<double> algorithm(data, params, seed_population);
   results = algorithm.run();
 }

 // Convert final_population to R list
 Rcpp::List final_population_list(results.final_population.size());
 for (size_t i = 0; i < results.final_population.size(); ++i) {
   final_population_list[i] = Rcpp::wrap(results.final_population[i]);
 }
 
 Rcpp::List parameters = Rcpp::List::create(
   Rcpp::Named("population_size") = population_size,
   Rcpp::Named("epochs") = epochs,
   Rcpp::Named("mutation_rate") = mutation_rate,
   Rcpp::Named("prob_mutation_type1") = prob_mutation_type1,
   Rcpp::Named("crossover_rate") = crossover_rate,
   Rcpp::Named("elite_count") = elite_count,
   Rcpp::Named("tournament_size") = tournament_size,
   Rcpp::Named("alpha") = alpha,
   Rcpp::Named("score_type") = score_type,
   Rcpp::Named("diversity") = diversity,
   Rcpp::Named("seed") =
     (params.seed >= 0 ? Rcpp::wrap(params.seed) : R_NilValue)
 );
 
 Rcpp::List statistics = Rcpp::List::create(
   Rcpp::Named("total_generations") = results.total_generations,
   Rcpp::Named("cache_hits") = results.cache_hits
 );
 
 return Rcpp::List::create(
   Rcpp::Named("final_population") = final_population_list,
   Rcpp::Named("final_scores") = Rcpp::wrap(results.final_scores),
   Rcpp::Named("parameters") = parameters,
   Rcpp::Named("statistics") = statistics
 );
} 


//' Run MCMC Algorithm to Compute the Score Distribution Among Node Combinations 
//' of Size 2
//'
//' Performs a Modified Metropolis-Hastings MCMC sampling to compute the score
//' distribution of nodes combination of size \emph{cocktail_size}.
//'
//' @param patient_data A data.frame containing patient information with at least
//'   a node column and a target column.
//' @param node_column Either a string (column name) or integer (column index, 1-based)
//'   specifying the column containing drug codes. This column should be either:
//'   \itemize{
//'     \item A list of integer vectors: \code{list(c(1,2), c(3), c(4,5))}
//'     \item A character vector with comma-separated values: \code{c("1,2", "3", "4,5")}
//'   }
//' @param id_column Optional observation-unit identifier column, given by name
//'   or one-based column position; needed for patient-level scoring.
//' @param target_column Either a string (column name) or integer (column index, 1-based)
//'   specifying the target/outcome column. Integer values are treated as binary for now,
//'   numeric values with non-0/1 entries are treated as continuous.
//' @param tree A data.frame containing the structural definition of the tree.
//' @param depth_column Either a string or integer specifying the column in 
//'   \code{tree} that contains the node depth levels.
//' @param upper_bound_column (Optional) Either a string or integer (1-based index) 
//'  specifying the column in \code{tree_depth} that contains upper 
//'  bound of nodes. Defaults to \code{NULL}.
//' @param name_column (Optional) Either a string or integer (1-based index) 
//'  specifying the column in \code{tree} that contains the corresponding name
//'  of nodes. Defaults to \code{NULL}.
//' @param beta Strict support threshold: filtered results require more than
//'   beta covered observations (or distinct patients for patient-level scores). Default: 4.
//' @param max_score Finite positive cap applied to scores before Metropolis-Hastings
//'   acceptance, importance weighting and binning. Default: 200.0.
//' @param score_type Scoring function to use. Either "hypergeometric" for the
//'   hypergeometric test, "relative_risk" for relative risk calculation, or "wilcoxon"
//'   for the wilcoxon test with continuous output.
//'   Default: "hypergeometric".
//'
//' @return A list containing:
//'   \describe{
//'     \item{top_solutions}{List of zero-based node vectors for the top scoring solutions}
//'     \item{top_scores}{Numeric vector of scores for the top solutions}
//'     \item{top_solutions_filtered}{Top solutions meeting the beta threshold}
//'     \item{top_scores_filtered}{Scores for the filtered solutions}
//'     \item{score_distribution}{Histogram of retained capped scores (0.1-wide bins and a cap bin)}
//'     \item{score_distribution_filtered}{Histogram for solutions meeting beta threshold}
//'     \item{outstanding_scores}{Retained capped scores equal to max_score}
//'     \item{uniform_reference}{Uniform-reference distribution and weight ESS}
//'     \item{uniform_reference_filtered}{Support-filtered uniform reference and ESS}
//'     \item{reference_parameters}{Temperature, score cap, support and burn-in settings}
//'     \item{trace}{Retained trajectory, only when store_trace is TRUE}
//'     \item{statistics}{List of run statistics including acceptance rates}
//'   }
//'
//' @details
//' Enumerates each unordered pair of distinct nodes exactly once and retains
//' pairs present in at least one observation, including ancestor-descendant
//' pairs. It uses the same state space, cap and strict support filter as
//' \code{run_mcmc_df_tree(cocktail_size = 2)}. No MCMC or importance weighting
//' is performed: each enumerated pair has unit weight. The returned
//' \code{uniform_reference} is exact; its \code{weight_ess} equals the number
//' of enumerated pairs. Acceptance statistics are not applicable.
//'
//' @examples
//' dat <- data.frame(outcome = c(1L, 1L, 0L, 0L))
//' dat$nodes <- list(c(0L, 1L), 0L, 1L, 1L)
//' tree <- data.frame(Depth = c(1L, 1L))
//' result <- mcmc_size_2_true_score_distribution(dat, "nodes", "outcome", tree, "Depth")
//' result$uniform_reference$weight_ess
//'
//' @export
// [[Rcpp::export]]
Rcpp::List mcmc_size_2_true_score_distribution(
   Rcpp::DataFrame patient_data,
   SEXP node_column,
   SEXP target_column,
   Rcpp::DataFrame tree,
   SEXP depth_column,
   SEXP upper_bound_column = R_NilValue,
   SEXP name_column = R_NilValue,
   SEXP id_column = R_NilValue,
   double beta = 4,
   double max_score = 200.0,
   std::string score_type = "hypergeometric") {
 
 // Validate inputs
 if (patient_data.nrows() == 0) {
   Rcpp::stop("patient_data cannot be empty");
 }
 if (tree.nrow() == 0) {
   Rcpp::stop("tree_depth cannot be empty");
 }
 
 
 // Build tree structure
 tree_structure cppTree(tree, depth_column, upper_bound_column, name_column);
 
 // Setup MCMC parameters
 MCMCParams params;
 params.cocktail_size = 2;
 params.beta = mcmc_count(beta, "beta");
 params.max_score = max_score;
 params.score_type_ = parse_score_type(score_type);
 
 // Detect target type and run appropriate template
 TargetTypeDetected target_type = detect_target_type(patient_data, target_column);
 MCMCResults results;
 
 if (target_type == TargetTypeDetected::BINARY) {
   PatientData<int> data(patient_data, node_column, target_column, cppTree, id_column);
   MCMCAlgorithm<int> algorithm(data, params);
   results = algorithm.true_size2_distribution();
 } else {
   PatientData<double> data(patient_data, node_column, target_column, cppTree,
                            id_column);
   MCMCAlgorithm<double> algorithm(data, params);
   results = algorithm.true_size2_distribution();
 }
 
 return mcmc_output(results, params);
}


//' Compute scores for supplied node combinations
//'
//' @param cocktail_list List of integer vectors identifying tree rows with
//'   one-based R indices. The returned \code{solutions} retain these indices.
//' @param patient_data Observation data frame. Its node list-column uses the
//'   package's zero-based internal tree indices.
//' @param node_column Name or one-based position of the node list-column.
//' @param target_column Name or one-based position of the outcome column.
//' @param tree Tree data frame in depth-first order.
//' @param depth_column Name or one-based position of the tree-depth column.
//' @param id_column Optional name or one-based position of the observation-unit
//'   identifier. It is required by patient-level continuous-outcome scores.
//' @param upper_bound_column Optional name or one-based position of the
//'   zero-based inclusive subtree upper-bound column.
//' @param name_column Optional name or one-based position of the node-label column.
//' @param score_type Registered score implementation to evaluate.
//' @return A list containing the supplied combinations, their scores, coverage
//'   counts, and score-specific summaries.
//' @details Candidate vectors deliberately use ordinary one-based R row
//'   positions at this user-facing boundary. They are checked against the tree
//'   before being converted once to the zero-based representation used by the
//'   C++ scoring engine.
//' @export
// [[Rcpp::export]]
Rcpp::List compute_score(
   Rcpp::List cocktail_list,
   Rcpp::DataFrame patient_data,
   SEXP node_column,
   SEXP target_column,
   Rcpp::DataFrame tree,
   SEXP depth_column,
   SEXP id_column = R_NilValue,
   SEXP upper_bound_column = R_NilValue,
   SEXP name_column = R_NilValue,
   std::string score_type = "hypergeometric"
){
 tree_structure cppTree(tree, depth_column, upper_bound_column, name_column);
 TargetTypeDetected target_type = detect_target_type(patient_data, target_column);
 ScoreType Cpp_score_type = parse_score_type(score_type);
 if (target_type == TargetTypeDetected::CONTINUOUS && Rf_isNull(id_column)) {
   Rcpp::stop("id_column is required for continuous-outcome scores.");
 }
 std::vector<Solution> sols;
 sols.reserve(cocktail_list.size());
 
 std::vector<double> cocktail_scores;
 cocktail_scores.reserve(cocktail_list.size());
 std::vector<double> cocktail_taker;
 cocktail_taker.reserve(cocktail_list.size());
 std::vector<std::vector<double>> diff_QT_distribution;
 diff_QT_distribution.reserve(cocktail_list.size());
 
 for (R_xlen_t cocktail_index = 0; cocktail_index < cocktail_list.size();
      ++cocktail_index) {
   std::vector<int> nodes =
     Rcpp::as<std::vector<int>>(cocktail_list[cocktail_index]);

   for (int node : nodes) {
     if (node == NA_INTEGER || node < 1 ||
         static_cast<size_t>(node) > cppTree.size()) {
       Rcpp::stop(
         "cocktail_list[[%i]] contains an invalid tree row; indices must lie "
         "between 1 and %i.",
         static_cast<int>(cocktail_index + 1),
         static_cast<int>(cppTree.size())
       );
     }
   }
   
   // subtract 1 from every element (1-based to 0-based)
   std::transform(nodes.begin(), nodes.end(), nodes.begin(), 
                  [](int x) { return x - 1; });
   
   sols.emplace_back(std::move(nodes));
 }
 
 if (target_type == TargetTypeDetected::BINARY) {
   PatientData<int> data(patient_data, node_column, target_column, cppTree);
   
   for(const auto& solution : sols){
     ScoreFunctions<int>::ScoreData score_data;
     switch(Cpp_score_type){
     case ScoreType::HYPERGEOMETRIC :
       score_data = ScoreFunctions<int>::compute_hypergeometric_with_data(data, solution);
       break;
     case ScoreType::RELATIVE_RISK :
       score_data = ScoreFunctions<int>::compute_relative_risk_with_data(data, solution);
       break;
     default :
       Rcpp::stop("Unknown score type");
     }
     cocktail_scores.push_back(score_data.score);
     cocktail_taker.push_back(score_data.covered_patients);
   }
   
 } else {
   PatientData<double> data(patient_data, node_column, target_column, cppTree, 
                            id_column);
   
   for(const auto& solution : sols){
     std::pair<ScoreFunctions<double>::ScoreData,std::vector<double>> score_data;
     switch(Cpp_score_type){
     case ScoreType::WILCOXON :
       score_data = ScoreFunctions<double>::compute_wilcoxon_risk_with_stats(data, solution);
       break;
     case ScoreType::COMPOSITE :
       score_data = ScoreFunctions<double>::compute_multifactor_risk_QT_with_stats(
         data, solution
       );
       break;
     case ScoreType::RESIDUALS :
       score_data = ScoreFunctions<double>::compute_residuals_risk_with_stats(
         data, solution
       );
       break;
     case ScoreType::BOOTSTRAP : 
       score_data = ScoreFunctions<double>::compute_residuals_risk_bootstrap_with_stats(
         data, solution
       );
       break;
       
     default :
       Rcpp::stop("Wrong score type");
     }
     cocktail_scores.push_back(score_data.first.score);
     cocktail_taker.push_back(score_data.first.covered_patients);
     diff_QT_distribution.emplace_back(score_data.second);
   }
 }
 return Rcpp::List::create(
   Rcpp::Named("solutions") = cocktail_list,
   Rcpp::Named("scores") = Rcpp::wrap(cocktail_scores),
   Rcpp::Named("number of takers") = Rcpp::wrap(cocktail_taker),
   Rcpp::Named("QT_diff_distribution") = Rcpp::wrap(diff_QT_distribution)
 );
}


//' Compute the dissimilarity matrix on a list of cocktails
//'
//' @inheritParams compute_score
//' @return A symmetric numeric matrix of normalized greedy tree-edit dissimilarities.
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix get_dissimilarity_of_list(
    Rcpp::List cocktail_list,
    Rcpp::DataFrame patient_data,
    SEXP node_column,
    SEXP target_column,
    Rcpp::DataFrame tree,
    SEXP depth_column,
    SEXP upper_bound_column = R_NilValue,
    SEXP name_column = R_NilValue,
    std::string score_type = "wilcoxon"
){
  
  tree_structure cppTree(tree, depth_column, upper_bound_column, name_column);
  
  // Setup GA parameters
  GAParams params;
  params.population_size = cocktail_list.size();
  
  // Detect target type and run appropriate template
  TargetTypeDetected target_type = detect_target_type(patient_data, target_column);
  
  Rcpp::NumericMatrix results;
  
  if (target_type == TargetTypeDetected::BINARY) {
    PatientData<int> data(patient_data, node_column, target_column, cppTree);
    GeneticAlgorithm<int> algorithm(data, params, cocktail_list);
    results = algorithm.population_dissimilarity();
  } else {
    PatientData<double> data(patient_data, node_column, target_column, cppTree);
    GeneticAlgorithm<double> algorithm(data, params, cocktail_list);
    results = algorithm.population_dissimilarity();
  }
  return results;
}

//' Identify observations covered by hierarchical combinations
//'
//' @inheritParams compute_score
//' @param cocktail_list List of zero-based node-index vectors.
//' @param id_column Name of the integer patient identifier column.
//' @param hadm_column Name of the integer admission identifier column.
//' @return A list with the supplied combinations (`cocktail_list`), distinct
//'   patient identifiers (`Takers_set`), admission identifiers (`hadm_set`),
//'   and one-based observation rows (`idx_set`) for each combination.
//' @export
// [[Rcpp::export]]
Rcpp::List get_taker(
    Rcpp::List cocktail_list,
    Rcpp::DataFrame patient_data,
    SEXP node_column,
    SEXP target_column,
    SEXP id_column,
    SEXP hadm_column,
    Rcpp::DataFrame tree,
    SEXP depth_column,
    SEXP upper_bound_column = R_NilValue,
    SEXP name_column = R_NilValue
){
  
  tree_structure cppTree(tree, depth_column, upper_bound_column, name_column);
  TargetTypeDetected target_type = detect_target_type(patient_data, target_column);
  
  std::vector<int> ids = Rcpp::as<std::vector<int>>(
    patient_data[Rcpp::as<Rcpp::String>(id_column)]
  );
  std::vector<int> hadms = Rcpp::as<std::vector<int>>(
    patient_data[Rcpp::as<Rcpp::String>(hadm_column)]
  );
  
  
  std::vector<std::set<int>> ids_of_takers;
  std::vector<std::set<int>> hadm_id_of_takers;
  std::vector<std::set<int>> idx_of_takers;
  
  auto process_data = [&](auto& data) {
    for(int i = 0; i < cocktail_list.size(); ++i){
      std::set<int> current_taker_set;
      std::set<int> current_hadm_set;
      std::set<int> current_idx_set;
      std::vector<int> cocktail_crt = Rcpp::as<std::vector<int>>(cocktail_list[i]);
      
      for(int j = 0; j < data.size(); ++j){
        if(data.patient_has_combination(j, cocktail_crt)){
          current_taker_set.insert(ids[j]);
          current_hadm_set.insert(hadms[j]);
          current_idx_set.insert(j + 1);
        }
      }
      ids_of_takers.push_back(std::move(current_taker_set));
      hadm_id_of_takers.push_back(std::move(current_hadm_set));
      idx_of_takers.push_back(std::move(current_idx_set));
    }
  };
  
  if (target_type == TargetTypeDetected::BINARY) {
    PatientData<int> data(patient_data, node_column, target_column, cppTree);
    process_data(data);
  } else {
    PatientData<double> data(patient_data, node_column, target_column, cppTree);
    process_data(data);
  }
  
  return Rcpp::List::create(
    Rcpp::Named("cocktail_list") = cocktail_list,
    Rcpp::Named("Takers_set") = Rcpp::wrap(ids_of_takers),
    Rcpp::Named("hadm_set") = Rcpp::wrap(hadm_id_of_takers),
    Rcpp::Named("idx_set") = Rcpp::wrap(idx_of_takers)
  );
}
