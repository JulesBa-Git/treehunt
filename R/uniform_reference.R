#' Effective sample size based on importance weights
#'
#' Measures weight concentration using \eqn{(\sum_i w_i)^2 / \sum_i w_i^2}.
#' Multiplying all weights by a common positive constant leaves this value
#' unchanged. Computation is rescaled to avoid overflow and underflow.
#'
#' @param weights Non-negative finite weights, or log weights when `log = TRUE`.
#'   Log weights may include `-Inf` (zero weight).
#' @param log Whether `weights` are on the natural-log scale.
#' @return A numeric scalar between 1 and the number of supplied weights.
#'   Returns `NA_real_` for an empty vector; all-zero weights are an error.
#' @details This diagnostic measures the concentration of importance weights.
#'   It does not account for serial dependence in an MCMC trajectory and is not
#'   a convergence diagnostic or an estimate of the number of independent draws.
#' @examples
#' weight_ess(c(1, 1, 1))
#' weight_ess(c(1, 0.25, 0.25, 0.25, 0.25))
#' weight_ess(c(-1000, -1001), log = TRUE)
#' @export
weight_ess <- function(weights, log = FALSE) {
  if (!is.numeric(weights) || !is.null(dim(weights)) || anyNA(weights)) {
    stop("weights must be a numeric vector without missing values.", call. = FALSE)
  }
  if (!is.logical(log) || length(log) != 1L || is.na(log)) {
    stop("log must be TRUE or FALSE.", call. = FALSE)
  }
  if (!length(weights)) return(NA_real_)
  if (log) {
    if (any(weights == Inf) || all(weights == -Inf)) {
      stop("Log weights must exclude +Inf and include a finite value.", call. = FALSE)
    }
    scaled <- exp(weights - max(weights))
  } else {
    if (any(!is.finite(weights)) || any(weights < 0) || max(weights) == 0) {
      stop("weights must be finite, non-negative and not all zero.", call. = FALSE)
    }
    scaled <- weights / max(weights)
  }
  min(length(weights), max(1, sum(scaled)^2 / sum(scaled^2)))
}

#' Recover the uniform score reference from a tempered MCMC trajectory
#'
#' Estimates the distribution of scores when combinations are sampled uniformly
#' from the state space explored by a chain targeting
#' \eqn{f_T(C) \propto \exp(S(C)/T)}. Scores are weighted before any aggregation.
#'
#' @param scores Numeric vector of scores from successive retained MCMC states,
#'   including repeats after rejection. Remove burn-in before calling this function.
#' @param temperature Finite, positive temperature used by the chain.
#' @param max_score Upper score cap used by the chain. Use `Inf` only for a chain
#'   that used uncapped scores. For `run_mcmc()` output, scores in `trace` are
#'   already capped, so this argument may be omitted.
#' @return A list with `distribution` (columns `score`, `probability`, `cdf`, and
#'   inclusive `upper_tail`), normalized per-iteration `weights`, `n_samples`, and
#'   `weight_ess`. Empty input returns an empty distribution and `NA_real_` ESS.
#' @details With \eqn{s_i = \min(S_i, \mathrm{max\_score})}, the normalized weights
#'   are proportional to \eqn{\exp(-s_i/T)}. For a threshold \eqn{x}, the uniform
#'   upper-tail estimate is \eqn{\sum_i w_i 1(s_i \ge x)/\sum_i w_i}.
#'   This is a self-normalized importance-sampling estimate; its consistency
#'   requires an ergodic chain with the stated stationary distribution.
#'
#'   Uniformity refers to combinations, not to distinct score values. Retaining
#'   only accepted moves or unique combinations changes the sampling law and
#'   invalidates these weights. Legacy binned counts alone are insufficient for
#'   exact reweighting because scores within a bin can have different weights.
#'   The returned distribution concerns capped scores when a cap is used.
#'
#'   Filtering the retained trajectory by a support criterion estimates the
#'   uniform reference conditional on that criterion. Compute its ESS using only
#'   those retained iterations. See [weight_ess()] for the diagnostic's scope.
#' @examples
#' # Two combinations with scores 0 and log(4) have target probabilities 1/5, 4/5.
#' ref <- uniform_score_reference(c(0, rep(log(4), 4)), temperature = 1)
#' ref$distribution # Equal estimated mass under the uniform combination law
#' ref$weight_ess
#' @seealso [run_mcmc()], [run_mcmc_df_tree()], [run_pwp_mcmc()]
#' @export
uniform_score_reference <- function(scores, temperature = 1, max_score = Inf) {
  if (!is.numeric(scores) || !is.null(dim(scores)) || anyNA(scores)) {
    stop("scores must be a numeric vector without missing values.", call. = FALSE)
  }
  if (!is.numeric(temperature) || length(temperature) != 1L ||
      !is.finite(temperature) || temperature <= 0) {
    stop("temperature must be finite and positive.", call. = FALSE)
  }
  if (!is.numeric(max_score) || length(max_score) != 1L ||
      is.na(max_score) || max_score <= 0) {
    stop("max_score must be positive (possibly Inf).", call. = FALSE)
  }
  scores <- pmin(scores, max_score)
  if (any(!is.finite(scores))) {
    stop("Scores must be finite after applying max_score.", call. = FALSE)
  }
  n <- length(scores)
  distribution <- data.frame(score = numeric(), probability = numeric(),
                             cdf = numeric(), upper_tail = numeric())
  if (!n) return(list(distribution = distribution, weights = numeric(),
                      n_samples = 0L, weight_ess = NA_real_))
  delta <- (scores - min(scores)) / temperature
  overflow <- is.infinite(scores - min(scores))
  if (any(overflow)) delta[overflow] <- scores[overflow] / temperature - min(scores) / temperature
  log_weights <- -delta
  scaled <- exp(log_weights)
  weights <- scaled / sum(scaled)
  values <- sort(unique(scores))
  masses <- as.numeric(rowsum(weights, group = match(scores, values), reorder = TRUE))
  cdf <- pmin(1, cumsum(masses))
  tails <- pmin(1, rev(cumsum(rev(masses))))
  cdf[length(cdf)] <- 1
  tails[1L] <- 1
  list(distribution = data.frame(score = values, probability = masses,
                                  cdf = cdf, upper_tail = tails),
       weights = weights, n_samples = n, weight_ess = weight_ess(log_weights, log = TRUE))
}
