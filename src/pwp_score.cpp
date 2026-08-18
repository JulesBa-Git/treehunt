#include "pwp_score.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <numeric>
#include <stdexcept>
#include <unordered_map>

namespace {

template <typename T>
void require_length(const std::vector<T>& x, std::size_t n,
                    const std::string& name) {
  if (x.size() != n) {
    Rcpp::stop("PWP context field '%s' has length %i; expected %i.",
               name, static_cast<int>(x.size()), static_cast<int>(n));
  }
}

std::vector<int> intersect_sorted(const std::vector<int>& left,
                                  const std::vector<int>& right) {
  std::vector<int> out;
  out.reserve(std::min(left.size(), right.size()));
  std::set_intersection(left.begin(), left.end(), right.begin(), right.end(),
                        std::back_inserter(out));
  return out;
}

} // namespace

PWPScoreContext::PWPScoreContext(const PatientData<int>& data,
                                 const Rcpp::List& context)
    : data_(data), n_rows_(data.size()), n_covariates_(0),
      positive_only_(true), min_covered_patients_(1),
      min_covered_events_(1), patient_marker_(0) {
  if (!context.containsElementNamed("time") ||
      !context.containsElementNamed("status") ||
      !context.containsElementNamed("strata") ||
      !context.containsElementNamed("risk") ||
      !context.containsElementNamed("x") ||
      !context.containsElementNamed("nuisance_variance")) {
    Rcpp::stop("Incomplete PWP score context.");
  }

  const std::vector<double> time =
      Rcpp::as<std::vector<double>>(context["time"]);
  status_ = Rcpp::as<std::vector<int>>(context["status"]);
  const std::vector<int> strata =
      Rcpp::as<std::vector<int>>(context["strata"]);
  risk_ = Rcpp::as<std::vector<double>>(context["risk"]);
  Rcpp::NumericMatrix x_matrix(context["x"]);
  Rcpp::NumericMatrix nuisance_variance(context["nuisance_variance"]);

  if (context.containsElementNamed("positive_only")) {
    positive_only_ = Rcpp::as<bool>(context["positive_only"]);
  }
  if (context.containsElementNamed("min_covered_patients")) {
    min_covered_patients_ =
        Rcpp::as<std::size_t>(context["min_covered_patients"]);
  }
  if (context.containsElementNamed("min_covered_events")) {
    min_covered_events_ =
        Rcpp::as<std::size_t>(context["min_covered_events"]);
  }

  require_length(time, n_rows_, "time");
  require_length(status_, n_rows_, "status");
  require_length(strata, n_rows_, "strata");
  require_length(risk_, n_rows_, "risk");

  if (static_cast<std::size_t>(x_matrix.nrow()) != n_rows_) {
    Rcpp::stop("PWP context matrix 'x' must have one row per analysis row.");
  }
  n_covariates_ = static_cast<std::size_t>(x_matrix.ncol());
  if (nuisance_variance.nrow() != nuisance_variance.ncol() ||
      static_cast<std::size_t>(nuisance_variance.nrow()) != n_covariates_) {
    Rcpp::stop("'nuisance_variance' must be a p by p matrix matching 'x'.");
  }

  x_.resize(n_rows_ * n_covariates_);
  for (std::size_t i = 0; i < n_rows_; ++i) {
    for (std::size_t j = 0; j < n_covariates_; ++j) {
      x_[i * n_covariates_ + j] = x_matrix(i, j);
    }
  }
  nuisance_variance_.resize(n_covariates_ * n_covariates_);
  for (std::size_t i = 0; i < n_covariates_; ++i) {
    for (std::size_t j = 0; j < n_covariates_; ++j) {
      nuisance_variance_[i * n_covariates_ + j] =
          nuisance_variance(i, j);
    }
  }

  for (std::size_t i = 0; i < n_rows_; ++i) {
    if (!std::isfinite(time[i]) || time[i] <= 0.0 ||
        !std::isfinite(risk_[i]) || risk_[i] <= 0.0) {
      Rcpp::stop("PWP times and risk weights must be finite and positive.");
    }
    if (status_[i] != 0 && status_[i] != 1) {
      Rcpp::stop("PWP status must contain only 0 and 1.");
    }
    for (std::size_t j = 0; j < n_covariates_; ++j) {
      if (!std::isfinite(x_[i * n_covariates_ + j])) {
        Rcpp::stop("PWP nuisance design matrix contains non-finite values.");
      }
    }
  }

  // Build exact node/subtree coverage once. Patient rows contain 0-based node
  // indexes; each observed node covers itself and all its ancestors.
  const auto& tree = data_.get_tree();
  const auto& fathers = tree.get_father();
  node_rows_.resize(tree.size());
  for (std::size_t row = 0; row < n_rows_; ++row) {
    std::vector<int> ancestors;
    for (const int node : data_.get_patient_nodes(row)) {
      if (node < 0 || static_cast<std::size_t>(node) >= tree.size()) {
        Rcpp::stop("Patient node index %i is outside [0, %i].", node,
                   static_cast<int>(tree.size() - 1));
      }
      int current = node;
      while (current >= 0) {
        ancestors.push_back(current);
        current = fathers[current];
      }
    }
    std::sort(ancestors.begin(), ancestors.end());
    ancestors.erase(std::unique(ancestors.begin(), ancestors.end()),
                    ancestors.end());
    for (const int node : ancestors) {
      node_rows_[node].push_back(static_cast<int>(row));
    }
  }

  // Compress patient identifiers for exact coverage counts.
  std::unordered_map<int, std::size_t> patient_index;
  patient_index.reserve(n_rows_ / 2 + 1);
  row_patient_index_.resize(n_rows_);
  for (std::size_t i = 0; i < n_rows_; ++i) {
    const int id = data_.get_id(i);
    const auto inserted = patient_index.emplace(id, patient_index.size());
    row_patient_index_[i] = inserted.first->second;
  }
  patient_seen_.assign(patient_index.size(), 0);

  // Collect unique event times by stratum. Censored rows remain in the risk
  // sets but never create an event-time group.
  std::map<int, std::vector<double>> event_times_by_stratum;
  for (std::size_t i = 0; i < n_rows_; ++i) {
    if (status_[i] == 1) {
      event_times_by_stratum[strata[i]].push_back(time[i]);
    }
  }
  if (event_times_by_stratum.empty()) {
    Rcpp::stop("The PWP context contains no event.");
  }
  for (auto& item : event_times_by_stratum) {
    auto& times = item.second;
    std::sort(times.begin(), times.end());
    times.erase(std::unique(times.begin(), times.end()), times.end());
  }

  std::unordered_map<int, std::size_t> stratum_lookup;
  std::vector<std::vector<double>> event_times;
  event_times.reserve(event_times_by_stratum.size());
  std::size_t total_groups = 0;
  for (const auto& item : event_times_by_stratum) {
    stratum_lookup[item.first] = event_times.size();
    event_times.push_back(item.second);
    stratum_starts_.push_back(total_groups);
    total_groups += item.second.size();
    stratum_ends_.push_back(total_groups);
  }

  row_risk_end_group_.assign(n_rows_, -1);
  row_event_group_.assign(n_rows_, -1);
  std::vector<double> s0_end(total_groups, 0.0);
  std::vector<double> d0(total_groups, 0.0);
  std::vector<double> sx_end(total_groups * n_covariates_, 0.0);
  std::vector<double> dx(total_groups * n_covariates_, 0.0);
  std::vector<int> tied_events(total_groups, 0);

  for (std::size_t i = 0; i < n_rows_; ++i) {
    const auto stratum_it = stratum_lookup.find(strata[i]);
    if (stratum_it == stratum_lookup.end()) {
      continue; // a stratum without events contributes nothing
    }
    const std::size_t s = stratum_it->second;
    const auto& times = event_times[s];
    const auto upper = std::upper_bound(times.begin(), times.end(), time[i]);
    if (upper != times.begin()) {
      const std::size_t local_end =
          static_cast<std::size_t>(std::distance(times.begin(), upper) - 1);
      const std::size_t group = stratum_starts_[s] + local_end;
      row_risk_end_group_[i] = static_cast<int>(group);
      s0_end[group] += risk_[i];
      for (std::size_t j = 0; j < n_covariates_; ++j) {
        sx_end[group * n_covariates_ + j] +=
            risk_[i] * x_[i * n_covariates_ + j];
      }
    }
    if (status_[i] == 1) {
      const auto exact = std::lower_bound(times.begin(), times.end(), time[i]);
      if (exact == times.end() || *exact != time[i]) {
        Rcpp::stop("Internal error while matching an event time.");
      }
      const std::size_t group = stratum_starts_[s] +
          static_cast<std::size_t>(std::distance(times.begin(), exact));
      row_event_group_[i] = static_cast<int>(group);
      d0[group] += risk_[i];
      ++tied_events[group];
      for (std::size_t j = 0; j < n_covariates_; ++j) {
        dx[group * n_covariates_ + j] +=
            risk_[i] * x_[i * n_covariates_ + j];
      }
    }
  }

  std::vector<double> s0(total_groups, 0.0);
  std::vector<double> sx(total_groups * n_covariates_, 0.0);
  for (std::size_t s = 0; s < stratum_starts_.size(); ++s) {
    double cumulative0 = 0.0;
    std::vector<double> cumulative_x(n_covariates_, 0.0);
    for (std::size_t g = stratum_ends_[s]; g-- > stratum_starts_[s];) {
      cumulative0 += s0_end[g];
      s0[g] = cumulative0;
      for (std::size_t j = 0; j < n_covariates_; ++j) {
        cumulative_x[j] += sx_end[g * n_covariates_ + j];
        sx[g * n_covariates_ + j] = cumulative_x[j];
      }
    }
  }

  coeff_a_.assign(total_groups, 0.0);
  coeff_b_.assign(total_groups, 0.0);
  coeff_c_.assign(total_groups, 0.0);
  coeff_d_.assign(total_groups, 0.0);
  coeff_e_.assign(total_groups, 0.0);
  coeff_p0_.assign(total_groups * n_covariates_, 0.0);
  coeff_p1_.assign(total_groups * n_covariates_, 0.0);

  for (std::size_t g = 0; g < total_groups; ++g) {
    const int m = tied_events[g];
    if (m <= 0) {
      Rcpp::stop("Internal PWP event group without events.");
    }
    for (int l = 0; l < m; ++l) {
      const double fraction = static_cast<double>(l) / static_cast<double>(m);
      const double denominator = s0[g] - fraction * d0[g];
      if (!std::isfinite(denominator) || denominator <= 0.0) {
        Rcpp::stop("Non-positive Efron denominator in PWP context.");
      }
      const double inv = 1.0 / denominator;
      const double inv2 = inv * inv;
      coeff_a_[g] += inv;
      coeff_b_[g] += fraction * inv;
      coeff_c_[g] += inv2;
      coeff_d_[g] += fraction * inv2;
      coeff_e_[g] += fraction * fraction * inv2;
      for (std::size_t j = 0; j < n_covariates_; ++j) {
        const double nuisance_numerator =
            sx[g * n_covariates_ + j] -
            fraction * dx[g * n_covariates_ + j];
        coeff_p0_[g * n_covariates_ + j] +=
            nuisance_numerator * inv2;
        coeff_p1_[g * n_covariates_ + j] +=
            fraction * nuisance_numerator * inv2;
      }
    }
  }
}

std::vector<int>
PWPScoreContext::covered_rows(const Solution& solution) const {
  const auto& nodes = solution.get_nodes();
  if (nodes.empty()) {
    return {};
  }
  for (const int node : nodes) {
    if (node < 0 || static_cast<std::size_t>(node) >= node_rows_.size()) {
      Rcpp::stop("Solution node index is outside the tree.");
    }
    if (node_rows_[node].empty()) {
      return {};
    }
  }

  auto smallest = std::min_element(
      nodes.begin(), nodes.end(), [this](int left, int right) {
        return node_rows_[left].size() < node_rows_[right].size();
      });
  std::vector<int> rows = node_rows_[*smallest];
  for (const int node : nodes) {
    if (node == *smallest) {
      continue;
    }
    rows = intersect_sorted(rows, node_rows_[node]);
    if (rows.empty()) {
      break;
    }
  }
  return rows;
}

std::size_t PWPScoreContext::count_covered_patients(
    const std::vector<int>& rows) const {
  if (++patient_marker_ == 0) {
    std::fill(patient_seen_.begin(), patient_seen_.end(), 0);
    patient_marker_ = 1;
  }
  std::size_t count = 0;
  for (const int row : rows) {
    const std::size_t patient = row_patient_index_[row];
    if (patient_seen_[patient] != patient_marker_) {
      patient_seen_[patient] = patient_marker_;
      ++count;
    }
  }
  return count;
}

PWPScoreResult PWPScoreContext::compute(const Solution& solution) const {
  PWPScoreResult result{};
  const std::vector<int> rows = covered_rows(solution);
  result.covered_rows = rows.size();
  if (rows.empty()) {
    result.p_value = 1.0;
    return result;
  }
  result.covered_patients = count_covered_patients(rows);

  const std::size_t groups = coeff_a_.size();
  std::vector<double> sz_end(groups, 0.0);
  std::vector<double> dz(groups, 0.0);
  std::vector<double> sxz_end(groups * n_covariates_, 0.0);
  std::vector<double> dxz(groups * n_covariates_, 0.0);
  std::vector<int> event_z(groups, 0);

  for (const int row_int : rows) {
    const std::size_t row = static_cast<std::size_t>(row_int);
    const int risk_group = row_risk_end_group_[row];
    if (risk_group >= 0) {
      const std::size_t g = static_cast<std::size_t>(risk_group);
      sz_end[g] += risk_[row];
      for (std::size_t j = 0; j < n_covariates_; ++j) {
        sxz_end[g * n_covariates_ + j] +=
            risk_[row] * x_[row * n_covariates_ + j];
      }
    }
    const int event_group = row_event_group_[row];
    if (event_group >= 0) {
      const std::size_t g = static_cast<std::size_t>(event_group);
      ++event_z[g];
      ++result.covered_events;
      dz[g] += risk_[row];
      for (std::size_t j = 0; j < n_covariates_; ++j) {
        dxz[g * n_covariates_ + j] +=
            risk_[row] * x_[row * n_covariates_ + j];
      }
    }
  }

  std::vector<double> sz(groups, 0.0);
  std::vector<double> sxz(groups * n_covariates_, 0.0);
  for (std::size_t s = 0; s < stratum_starts_.size(); ++s) {
    double cumulative_z = 0.0;
    std::vector<double> cumulative_xz(n_covariates_, 0.0);
    for (std::size_t g = stratum_ends_[s]; g-- > stratum_starts_[s];) {
      cumulative_z += sz_end[g];
      sz[g] = cumulative_z;
      for (std::size_t j = 0; j < n_covariates_; ++j) {
        cumulative_xz[j] += sxz_end[g * n_covariates_ + j];
        sxz[g * n_covariates_ + j] = cumulative_xz[j];
      }
    }
  }

  double u = 0.0;
  double izz = 0.0;
  std::vector<double> ixz(n_covariates_, 0.0);
  for (std::size_t g = 0; g < groups; ++g) {
    const double expected_z = coeff_a_[g] * sz[g] - coeff_b_[g] * dz[g];
    u += static_cast<double>(event_z[g]) - expected_z;
    izz += expected_z -
        (coeff_c_[g] * sz[g] * sz[g] -
         2.0 * coeff_d_[g] * sz[g] * dz[g] +
         coeff_e_[g] * dz[g] * dz[g]);
    for (std::size_t j = 0; j < n_covariates_; ++j) {
      ixz[j] += coeff_a_[g] * sxz[g * n_covariates_ + j] -
          coeff_b_[g] * dxz[g * n_covariates_ + j] -
          sz[g] * coeff_p0_[g * n_covariates_ + j] +
          dz[g] * coeff_p1_[g * n_covariates_ + j];
    }
  }

  double nuisance_adjustment = 0.0;
  for (std::size_t i = 0; i < n_covariates_; ++i) {
    for (std::size_t j = 0; j < n_covariates_; ++j) {
      nuisance_adjustment += ixz[i] *
          nuisance_variance_[i * n_covariates_ + j] * ixz[j];
    }
  }
  const double efficient_information = izz - nuisance_adjustment;
  result.score = u;
  result.efficient_information = efficient_information;
  if (!std::isfinite(efficient_information) || efficient_information <= 1e-12) {
    result.p_value = 1.0;
    return result;
  }

  result.signed_z = u / std::sqrt(efficient_information);
  result.score_statistic = result.signed_z * result.signed_z;
  result.p_value = R::pchisq(result.score_statistic, 1.0, false, false);
  result.eligible =
      result.covered_patients >= min_covered_patients_ &&
      result.covered_events >= min_covered_events_;
  if (result.eligible) {
    result.fitness = positive_only_ ? std::max(0.0, result.signed_z)
                                    : result.signed_z;
  }
  return result;
}
