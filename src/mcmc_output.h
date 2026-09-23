#ifndef TREEHUNT_MCMC_OUTPUT_H
#define TREEHUNT_MCMC_OUTPUT_H
#include "mcmc_algorithm.h"
#include <climits>

inline size_t mcmc_count(double value, const char* name) {
  if (!std::isfinite(value) || value < 0 || value > INT_MAX || value != std::floor(value))
    Rcpp::stop(std::string(name) + " must be a non-negative integer no greater than .Machine$integer.max.");
  return static_cast<size_t>(value);
}

inline Rcpp::List mcmc_output(const MCMCResults& r, const MCMCParams& p) {
  Rcpp::List out = Rcpp::List::create(
    Rcpp::Named("top_solutions") = Rcpp::wrap(r.top_solutions),
    Rcpp::Named("top_scores") = Rcpp::wrap(r.top_scores),
    Rcpp::Named("top_solutions_filtered") = Rcpp::wrap(r.top_solutions_filtered),
    Rcpp::Named("top_scores_filtered") = Rcpp::wrap(r.top_scores_filtered),
    Rcpp::Named("score_distribution") = Rcpp::wrap(r.score_distribution),
    Rcpp::Named("score_distribution_filtered") = Rcpp::wrap(r.score_distribution_filtered),
    Rcpp::Named("outstanding_scores") = Rcpp::wrap(r.outstanding_scores),
    Rcpp::Named("uniform_reference") = r.uniform_reference.as_list(p.max_score),
    Rcpp::Named("uniform_reference_filtered") = r.uniform_reference_filtered.as_list(p.max_score),
    Rcpp::Named("statistics") = Rcpp::List::create(
      Rcpp::Named("total_iterations") = r.total_iterations,
      Rcpp::Named("accepted_moves") = r.accepted_moves,
      Rcpp::Named("rejected_moves") = r.rejected_moves,
      Rcpp::Named("acceptance_rate") = r.exhaustive || r.total_iterations == 0 ? NA_REAL :
        static_cast<double>(r.accepted_moves) / r.total_iterations,
      Rcpp::Named("proposals_not_in_population") = r.proposals_not_in_population,
      Rcpp::Named("type1_moves") = r.type1_moves,
      Rcpp::Named("type2_moves") = r.type2_moves,
      Rcpp::Named("type1_accepted") = r.type1_accepted,
      Rcpp::Named("type2_accepted") = r.type2_accepted,
      Rcpp::Named("type1_in_population") = r.type1_in_population,
      Rcpp::Named("type2_in_population") = r.type2_in_population,
      Rcpp::Named("cocktail_size") = r.cocktail_size,
      Rcpp::Named("seed") = p.seed >= 0 ? Rcpp::wrap(p.seed) : R_NilValue),
    Rcpp::Named("reference_parameters") = Rcpp::List::create(
      Rcpp::Named("method") = r.exhaustive ? "enumeration" : "self_normalized_importance_sampling",
      Rcpp::Named("temperature") = p.temperature,
      Rcpp::Named("max_score") = p.max_score,
      Rcpp::Named("beta") = p.beta,
      Rcpp::Named("burn_in") = p.burn_in,
      Rcpp::Named("support") = "observed subsets of distinct nodes of fixed cocktail_size; ancestor-descendant pairs included"));
  if (p.store_trace) {
    Rcpp::List trace = Rcpp::List::create(Rcpp::Named("score") = r.trace_scores,
      Rcpp::Named("covered_patients") = r.trace_support);
    trace["solution"] = Rcpp::wrap(r.trace_solutions);
    trace.attr("class") = "data.frame";
    trace.attr("row.names") = Rcpp::IntegerVector::create(NA_INTEGER,
      -static_cast<int>(r.trace_scores.size()));
    out["trace"] = trace;
  }
  return out;
}
#endif
