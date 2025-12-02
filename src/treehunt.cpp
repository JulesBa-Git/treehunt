#include <Rcpp.h>
#include "mcmc_algorithm.h"
#include "genetic_algorithm.h"
#include "patient_data.h"
#include "tree_structure.h"


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
  } else {
    Rcpp::stop("Unknown score type: '" + score_type_str + 
      "'. Use 'hypergeometric' or 'relative_risk'.");
  }
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
//' @param beta Minimum number of patients that must be covered for a solution to
//'   be included in the filtered results. Default: 4.
//' @param max_score Maximum score value for binning in the score distribution.
//'   Scores above this are tracked separately. Default: 200.0.
//' @param score_type Scoring function to use. Either "hypergeometric" for the
//'   hypergeometric test or "relative_risk" for relative risk calculation.
//'   Default: "hypergeometric".
//' @param verbose If TRUE, prints progress and statistics during the run.
//'   Default: FALSE.
//'
//' @return A list containing:
//'   \describe{
//'     \item{top_solutions}{List of node vectors for the top scoring solutions}
//'     \item{top_scores}{Numeric vector of scores for the top solutions}
//'     \item{top_solutions_filtered}{Top solutions meeting the beta threshold}
//'     \item{top_scores_filtered}{Scores for the filtered solutions}
//'     \item{score_distribution}{Histogram of scores (0.1-wide bins)}
//'     \item{score_distribution_filtered}{Histogram for solutions meeting beta threshold}
//'     \item{outstanding_scores}{Scores that exceeded max_score}
//'     \item{statistics}{List of run statistics including acceptance rates}
//'   }
//'
//' @details
//' The MCMC algorithm uses a Modified Metropolis-Hastings approach with two
//' proposal types:
//' \itemize{
//'   \item \strong{Type 1}: Generates a completely new random valid solution
//'   \item \strong{Type 2}: Swaps one node with its parent or child in the tree
//' }
//'
//' The acceptance probability for Type 1 proposal is:
//' \deqn{\alpha = \exp((S_{proposed} - S_{current}) / T)}
//'
//' For Type 2 proposals, a proposal ratio correction is applied since the ratio
//' of \mathbb{P}(current | proposed) \noteq \mathbb{P}(proposed | current):
//' \deqn{\alpha = \exp((S_{proposed} - S_{current}) / T) \times \frac{|V_{current}|}{|V_{proposed}|}}
//'
//' where \eqn{|V|} is the number of possible swap vertices for a solution.
//'
//' Solutions are only accepted if they appear in at least one patient's data
//' ("modified" constraint).
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
//' # Define tree structure (simple 3-level tree)
//' tree_depth <- c(1, rep(2, 5), rep(3, 15))
//'
//' # Run MCMC
//' results <- run_mcmc(
//'   patient_data = patient_df,
//'   node_column = "drugs",
//'   target_column = "outcome",
//'   tree_depth = tree_depth,
//'   epochs = 5000,
//'   cocktail_size = 2,
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
// [[Rcpp::export]]
Rcpp::List run_mcmc(
   Rcpp::DataFrame patient_data,
   SEXP node_column,
   SEXP target_column,
   Rcpp::IntegerVector tree_depth,
   size_t epochs,
   double temperature = 1.0,
   size_t n_results = 10,
   size_t cocktail_size = 2,
   double prob_type1 = 0.01,
   size_t beta = 4,
   double max_score = 200.0,
   std::string score_type = "hypergeometric",
   bool verbose = false) {
 
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
 params.epochs = epochs;
 params.temperature = temperature;
 params.n_results = n_results;
 params.cocktail_size = cocktail_size;
 params.prob_mutation_type1 = prob_type1;
 params.beta = beta;
 params.max_score = max_score;
 params.score_type_ = parse_score_type(score_type);
 params.verbose = verbose;
 
 // Detect target type and run appropriate template
 TargetTypeDetected target_type = detect_target_type(patient_data, target_column);
 MCMCResults results;
 
 if (target_type == TargetTypeDetected::BINARY) {
   PatientData<int> data(patient_data, node_column, target_column, tree);
   MCMCAlgorithm<int> algorithm(data, params);
   results = algorithm.run();
 } else {
   PatientData<double> data(patient_data, node_column, target_column, tree);
   MCMCAlgorithm<double> algorithm(data, params);
   results = algorithm.run();
 }
 
 // Convert top_solutions to R list
 Rcpp::List top_solutions_list(results.top_solutions.size());
 for (size_t i = 0; i < results.top_solutions.size(); ++i) {
   top_solutions_list[i] = Rcpp::wrap(results.top_solutions[i]);
 }
 
 Rcpp::List top_solutions_filtered_list(results.top_solutions_filtered.size());
 for (size_t i = 0; i < results.top_solutions_filtered.size(); ++i) {
   top_solutions_filtered_list[i] = Rcpp::wrap(results.top_solutions_filtered[i]);
 }
 
 // Build statistics list
 Rcpp::List statistics = Rcpp::List::create(
   Rcpp::Named("total_iterations") = results.total_iterations,
   Rcpp::Named("accepted_moves") = results.accepted_moves,
   Rcpp::Named("rejected_moves") = results.rejected_moves,
   Rcpp::Named("acceptance_rate") = static_cast<double>(results.accepted_moves) / 
     static_cast<double>(results.total_iterations),
     Rcpp::Named("proposals_not_in_population") = results.proposals_not_in_population,
     Rcpp::Named("type1_moves") = results.type1_moves,
     Rcpp::Named("type2_moves") = results.type2_moves,
     Rcpp::Named("type1_accepted") = results.type1_accepted,
     Rcpp::Named("type2_accepted") = results.type2_accepted,
     Rcpp::Named("type1_in_population") = results.type1_in_population,
     Rcpp::Named("type2_in_population") = results.type2_in_population,
     Rcpp::Named("cocktail_size") = results.cocktail_size
 );
 
 return Rcpp::List::create(
   Rcpp::Named("top_solutions") = top_solutions_list,
   Rcpp::Named("top_scores") = Rcpp::wrap(results.top_scores),
   Rcpp::Named("top_solutions_filtered") = top_solutions_filtered_list,
   Rcpp::Named("top_scores_filtered") = Rcpp::wrap(results.top_scores_filtered),
   Rcpp::Named("score_distribution") = Rcpp::wrap(results.score_distribution),
   Rcpp::Named("score_distribution_filtered") = Rcpp::wrap(results.score_distribution_filtered),
   Rcpp::Named("outstanding_scores") = Rcpp::wrap(results.outstanding_scores),
   Rcpp::Named("statistics") = statistics
 );
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
//'     \item A list of integer vectors: \code{list(c(1,2), c(3), c(4,5))}
//'     \item A character vector with comma-separated values: \code{c("1,2", "3", "4,5")}
//'   }
//' @param target_column Either a string (column name) or integer (column index, 1-based)
//'   specifying the target/outcome column. Integer values are treated as binary,
//'   numeric values with non-0/1 entries are treated as continuous.
//' @param tree_depth An integer vector specifying the depth of each node in the
//'   tree structure. Must start at depth 1 and children must be at depth+1 of
//'   their parent.
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
//'   hypergeometric test or "relative_risk" for relative risk calculation.
//'   Default: "hypergeometric".
//' @param diversity If TRUE, applies a diversity penalty to encourage exploration
//'   of different solutions. Default: FALSE.
//' @param verbose If TRUE, prints progress during the run. Default: FALSE.
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
//' \dontrun
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
   bool verbose = false) {
 
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
 
 // Detect target type and run appropriate template
 TargetTypeDetected target_type = detect_target_type(patient_data, target_column);
 GAResults results;
 
 if (target_type == TargetTypeDetected::BINARY) {
   PatientData<int> data(patient_data, node_column, target_column, tree);
   GeneticAlgorithm<int> algorithm(data, params);
   results = algorithm.run();
 } else {
   PatientData<double> data(patient_data, node_column, target_column, tree);
   GeneticAlgorithm<double> algorithm(data, params);
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
   Rcpp::Named("diversity") = diversity
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