#ifndef TREEHUNT_UNIFORM_REFERENCE_H
#define TREEHUNT_UNIFORM_REFERENCE_H

#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>

// Self-normalized importance weights, rescaled whenever a smaller score appears.
// Each retained iteration contributes, including repeated states after rejection.
class UniformReference {
  std::vector<double> bins_;
  double minimum_ = 0.0, sum_ = 0.0, sum_squared_ = 0.0;
  size_t count_ = 0;
public:
  void resize(size_t n) { bins_.assign(n, 0.0); }
  void add(size_t bin, double score, double temperature) {
    if (count_ == 0) minimum_ = score;
    if (score < minimum_) {
      const double scale = std::exp(-(minimum_ - score) / temperature);
      for (double& value : bins_) value *= scale;
      sum_ *= scale;
      sum_squared_ *= scale * scale;
      minimum_ = score;
    }
    const double weight = std::exp(-(score - minimum_) / temperature);
    bins_[bin] += weight;
    sum_ += weight;
    sum_squared_ += weight * weight;
    ++count_;
  }
  Rcpp::List as_list(double cap) const {
    const size_t n = bins_.size();
    Rcpp::NumericVector lower(n), upper(n), probability(n, NA_REAL),
      cdf(n, NA_REAL), tail(n, NA_REAL);
    for (size_t i = 0; i < n; ++i) {
      lower[i] = i + 1 == n ? cap : i / 10.0;
      upper[i] = std::min((i + 1) / 10.0, cap);
      if (count_) probability[i] = bins_[i] / sum_;
    }
    if (count_) {
      double cumulative = 0.0, survival = 0.0;
      for (size_t i = 0; i < n; ++i) {
        cumulative += probability[i];
        cdf[i] = std::min(1.0, cumulative);
        survival += probability[n - 1 - i];
        tail[n - 1 - i] = std::min(1.0, survival);
      }
      cdf[n - 1] = 1.0;
      tail[0] = 1.0;
    }
    const double ess = count_ ? std::min(static_cast<double>(count_),
      std::max(1.0, sum_ * sum_ / sum_squared_)) : NA_REAL;
    return Rcpp::List::create(
      Rcpp::Named("distribution") = Rcpp::DataFrame::create(
        Rcpp::Named("lower") = lower, Rcpp::Named("upper") = upper,
        Rcpp::Named("probability") = probability,
        Rcpp::Named("cdf") = cdf, Rcpp::Named("upper_tail") = tail),
      Rcpp::Named("n_samples") = static_cast<double>(count_),
      Rcpp::Named("weight_ess") = ess);
  }
};
#endif
