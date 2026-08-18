#ifndef PWP_SCORE_H
#define PWP_SCORE_H

#include <Rcpp.h>
#include <cstddef>
#include <cstdint>
#include <unordered_map>
#include <vector>

#include "patient_data.h"
#include "solution.h"

struct PWPScoreResult {
  bool eligible;
  double fitness;
  double signed_z;
  double score_statistic;
  double p_value;
  double score;
  double efficient_information;
  std::size_t covered_rows;
  std::size_t covered_patients;
  std::size_t covered_events;
};

// Precomputed context for the score test of one time-varying binary covariate
// added to a reduced, stratified PWP gap-time Cox model. The risk-set constants
// do not depend on the candidate combination and are built once.
class PWPScoreContext {
private:
  const PatientData<int>& data_;
  std::size_t n_rows_;
  std::size_t n_covariates_;
  bool positive_only_;
  std::size_t min_covered_patients_;
  std::size_t min_covered_events_;

  std::vector<double> risk_;
  std::vector<double> x_; // row-major: row * p + covariate
  std::vector<double> nuisance_variance_; // inverse nuisance information
  std::vector<int> status_;

  // For every tree node, sorted row indexes covered by that node/subtree.
  std::vector<std::vector<int>> node_rows_;

  // Event-time group representation. Groups are contiguous within strata.
  std::vector<std::size_t> stratum_starts_;
  std::vector<std::size_t> stratum_ends_;
  std::vector<int> row_risk_end_group_;
  std::vector<int> row_event_group_;

  // Efron coefficients, one value per event-time group.
  std::vector<double> coeff_a_;
  std::vector<double> coeff_b_;
  std::vector<double> coeff_c_;
  std::vector<double> coeff_d_;
  std::vector<double> coeff_e_;
  std::vector<double> coeff_p0_; // group * p + covariate
  std::vector<double> coeff_p1_; // group * p + covariate

  // Patient markers provide exact distinct-patient coverage without allocating
  // a hash set for every candidate.
  std::vector<std::size_t> row_patient_index_;
  mutable std::vector<std::uint32_t> patient_seen_;
  mutable std::uint32_t patient_marker_;

  std::vector<int> covered_rows(const Solution& solution) const;
  std::size_t count_covered_patients(const std::vector<int>& rows) const;

public:
  PWPScoreContext(const PatientData<int>& data, const Rcpp::List& context);

  PWPScoreResult compute(const Solution& solution) const;
  std::size_t number_event_groups() const { return coeff_a_.size(); }
  std::size_t number_covariates() const { return n_covariates_; }
};

#endif
