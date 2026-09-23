# treehunt

**treehunt** is an R package that implements evolutionary optimization algorithms to search for optimal node combinations in hierarchical tree structures. Originally developed for pharmacovigilance applications to identify drug (ATC code) combinations associated with adverse events, the package generalizes to any tree hierarchy and optimization metric.

## 🎯 Overview

Finding optimal combinations of elements from a hierarchical structure is a challenging combinatorial problem. **treehunt** addresses this by providing two algorithms:

- **MCMC (Modified Metropolis-Hastings)**: A stochastic search that explores the solution space through random walks in order to estimate the score distribution among a fixed size of node cocktails, accepting or rejecting proposals based on score improvements
- **Genetic Algorithm**: An evolutionary approach that evolves a population of solutions through selection, crossover, and mutation in order to find node combinations that optimize a score function

Both algorithms are implemented in C++ via Rcpp for high performance and support multiple scoring functions.

## 📦 Installation

### From GitHub (development version)

```r
# install.packages("devtools")
devtools::install_github("JulesBa-Git/treehunt")
```


## 🚀 Quick Start

### Basic Example

```r
library(treehunt)

# Create example patient data
set.seed(42)
n_patients <- 500

patient_df <- data.frame(
  patient_id = 1:n_patients,
  adverse_event = rbinom(n_patients, 1, 0.2)
)

# Each patient has a list of tree node indices
patient_df$drug_codes <- lapply(1:n_patients, function(i) {
  sample(5:25, size = sample(1:4, 1), replace = FALSE)
})

# Define tree structure: depth vector where index = node, value = depth
# Example: 3-level tree with 1 root, 5 level-2 nodes, each with 4 level-3 children
 tree_depth <- c(1,  # Node 0: root
				 2,  3,  3,  3,  3,  # Node 1 (child of 0) + its 4 children
				 2,  3,  3,  3,  3,  # Node 6 (child of 0) + its 4 children
				 2,  3,  3,  3,  3,  # Node 11 (child of 0) + its 4 children 
				 2,  3,  3,  3,  3,  # Node 16 (child of 0) + its 4 children
				 2,  3,  3,  3,  3)  # Node 21 (child of 0) + its 4 children

# Run MCMC search
mcmc_results <- run_mcmc(

  patient_data = patient_df,
  node_column = "drug_codes", # or by index : 3
  target_column = "adverse_event", # or by index : 2
  tree_depth = tree_depth,
  epochs = 10000,
  cocktail_size = 2,
  score_type = "hypergeometric",
  seed = 4601L,
  verbose = TRUE
)

# View top results
print(mcmc_results$top_scores)
print(mcmc_results$top_solutions)
```

### Uniform reference and weight ESS

The MCMC target favours high scores: `f_T(C)` is proportional to
`exp(S(C) / temperature)`. Here `S` is capped at `max_score`, as in the
Metropolis-Hastings acceptance step. The reference below gives equal weight to
each observed combination of the chosen size; different score values need not
have equal probabilities.

```r
reference_run <- run_mcmc(
  patient_data = patient_df, node_column = "drug_codes",
  target_column = "adverse_event", tree_depth = tree_depth,
  epochs = 10000, burn_in = 1000, cocktail_size = 2,
  temperature = 2, max_score = 200, beta = 4, seed = 4601L,
  store_trace = TRUE
)

reference_run$uniform_reference$distribution
reference_run$uniform_reference$weight_ess
reference_run$uniform_reference_filtered$weight_ess

# Exact score values instead of histogram bins, using the optional trace:
reference <- uniform_score_reference(reference_run$trace$score, temperature = 2)
reference$distribution
```

At every retained iteration, including repeats after rejection, the inverse
target weight is proportional to `exp(-S / temperature)`. The engine computes
these weights before binning, then returns normalized bin probabilities, their
CDF and inclusive upper tails. The filtered reference conditions on coverage
strictly greater than `beta` (distinct patients for PWP). An empty filtered
reference returns `n_samples = 0` and `NA` probabilities and ESS.

`weight_ess = sum(w)^2 / sum(w^2)` measures weight concentration. It does not
account for MCMC autocorrelation or establish convergence. For large score
ranges relative to temperature, a longer run or a higher temperature may be
needed to cover low-score regions. Check mixing separately. With a finite cap,
the reference describes capped scores; scores beyond the cap are grouped together.

By default `store_trace = FALSE`: weighted histograms and ESS are accumulated
without storing the full trajectory. `burn_in` is a number of initial iterations
discarded from distributions, top solutions and trace. Acceptance statistics
still describe all iterations. Legacy histograms alone cannot be reweighted
exactly because they do not retain individual scores within each bin.

The MCMC state space contains observed fixed-size subsets of distinct nodes,
including ancestor-descendant pairs. Node vectors in observations, MCMC outputs
and traces are **zero-based** tree indices. The support condition and tree must
match when comparing a sampled reference with
`mcmc_size_2_true_score_distribution()`; that function enumerates supported pairs
with unit weights. The GA enforces its own ancestor-descendant validity rule.

### Using the Genetic Algorithm

```r
# Run genetic algorithm search
ga_results <- run_genetic_algorithm(
  patient_data = patient_df,
  node_column = "drug_codes",
  target_column = "adverse_event",
  tree_depth = tree_depth,
  population_size = 100,
  epochs = 500,
  mutation_rate = 0.1,
  elite_count = 2,
  score_type = "hypergeometric",
  seed = 4602L,
  verbose = TRUE
)

# Inspect the ten highest scores and their zero-based node indices
top <- head(order(ga_results$final_scores, decreasing = TRUE), 10)
data.frame(
  score = ga_results$final_scores[top],
  nodes = vapply(ga_results$final_population[top], paste,
                 collapse = ",", FUN.VALUE = character(1))
)
```

The generic GA and MCMC wrappers seed their C++ random-number generators
directly. Supply `seed` to reproduce a stochastic run; the default `NULL`
preserves non-deterministic initialization. Calling `set.seed()` alone does not
set these C++ streams.

## 📊 Scoring Functions

### Hypergeometric Score

Tests whether the node (drug) combination is significantly over-represented among patients with the adverse event using the hypergeometric distribution:

$$ \text{score} = -\log(P(X ≥ k)) $$


where k is the number of patients with both the combination and the adverse event.

### Relative Risk

Calculates the ratio of adverse event probability in patients with the combination versus those without:

$$ RR =\frac{\mathbb{P}(\text{event}| \text{combination})}{\mathbb{P}(\text{event} | \text{no combination})} $$

### Want a Different Scoring Function? 

If you need a scoring function that isn't currently implemented, we'd appreciate to hear from you! You can:
 -  **Open an issue** describing your use case and the scoring function you need 
 -  **Submit a pull request** with your implementation

## 🌳 Tree Structure

The tree structure is defined by a depth vector where:
- Index represents the node ID (DFS algorithm)
- Value represents the depth level (1 = root son's, 2 = first children, etc.)
- Nodes must be ordered such that children immediately follow their parent
- A node following a node of depth $k$ can have a depth in $\{ 1, \dots,k+1\}$

**Example**: A tree with structure:
```
       0 (root)
      /|\
     1 4 5
    /|   |
   2 3   6
```

Would be defined as:
```r
tree_depth <- c(1, 2, 3, 3, 2, 2, 3)
```

### Using ATC Classification

For pharmacovigilance applications with ATC codes:

```r
# ATC codes have five levels (code lengths 1, 3, 4, 5, 7).
# Example path, already in depth-first order:
atc_data <- data.frame(code = c("A", "A01", "A01A", "A01AA", "A01AA01"))
atc_data$depth <- match(nchar(atc_data$code), c(1L, 3L, 4L, 5L, 7L))
```

## PWP Cox Rao searches

`treehunt` can rank hierarchical code combinations for recurrent-event data
with an exact Efron score test from a stratified PWP gap-time Cox model. Fit the
reduced model once, then reuse its context throughout the search:

```r
context <- fit_pwp_rao_context(
  data = pwp_data,
  status_column = "status_any", # or "status_emergency"
  covariates = c("age_start", "index_emergency", "year_start"),
  min_covered_patients = 20L,
  min_covered_events = 5L
)

ga <- run_pwp_genetic_algorithm(
  data = pwp_data,
  context = context,
  tree = icd10_tree,
  node_column = "icd10_indices",
  seed = 20260816L
)

details <- score_pwp_combinations(
  combinations = ga$final_population,
  data = pwp_data,
  context = context,
  tree = icd10_tree,
  index_base = 0L
)
```

The interval list-column and tree `UpperBound` must both use zero-based node
indexes. In particular, when a CSV stores an inclusive one-based subtree bound,
use `UpperBound <- UpperBound - 1L`. The screening fitness is
`max(0, signed_z)`; `details` retains the signed Rao statistic, model-based
p-value, coverage counts, and efficient information. Use
`refit_pwp_combinations()` to obtain patient-clustered robust Cox estimates for
the selected candidates.

## 📈 Output Structure

### MCMC Results

```r
mcmc_results <- run_mcmc(...)

mcmc_results$top_solutions          		# List of node vectors
mcmc_results$top_scores             		# Corresponding scores
mcmc_results$top_solutions_filtered 		# Solutions meeting beta threshold
mcmc_results$score_distribution     		# Histogram of visited scores
mcmc_results$score_distribution_filtered 	# Histogram of visited scores 
											# meeting beta threshold
mcmc_results$statistics             		# Acceptance rates, move counts
```

### GA Results

```r
ga_results <- run_genetic_algorithm(...)

ga_results$final_population # Complete final population
ga_results$statistics       # Generations, cache hits
```

## 🔬 Algorithm Details

For more details about the method you can watch our preprint:

[Bangard, J., Holsbø, E., Svendsen, K., Perduca, V., & Birmelé, E. (2025). Detecting adverse high-order drug interactions from individual case safety reports using computational statistics on disproportionality measures. _arXiv preprint arXiv:2504.00646_](https://arxiv.org/abs/2504.00646)


## 🤝 Contributing

Contributions are welcome! Please feel free to submit issues and pull requests.

1. Fork the repository
2. Create your feature branch (`git checkout -b feature/feature_name`)
3. Commit your changes (`git commit -m 'Add feature_name: details'`)
4. Push to the branch (`git push origin feature/feature_name`)
5. Open a Pull Request

## License

GPL-3. See the repository LICENSE file.
