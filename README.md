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
  verbose = TRUE
)

# View top results
print(mcmc_results$top_scores)
print(mcmc_results$top_solutions)
```

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
  verbose = TRUE
)

# Get summary of top solutions
top_solutions_summary(ga_results)
```

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
# ATC codes have 5 levels (lengths 1, 3, 4, 5, 7)
# Map to depths 1-5
atc_data$depth <- case_when(
  nchar(atc_data$code) == 1 ~ 1,
  nchar(atc_data$code) == 3 ~ 2,
  nchar(atc_data$code) == 4 ~ 3,
  nchar(atc_data$code) == 5 ~ 4,
  nchar(atc_data$code) == 7 ~ 5
)
```

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
