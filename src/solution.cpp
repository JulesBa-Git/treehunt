#include "solution.h"
#include "tree_structure.h"
#include "patient_data.h"

size_t Solution::compute_hash_internal() const{
  size_t seed = selected_nodes_.size();
  
  for(int node : selected_nodes_){
    // Boost hash_combine formula
    seed ^= std::hash<int>{}(node) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
  }
  
  return seed;
}

Solution::Solution(std::vector<int> nodes) : selected_nodes_{std::move(nodes)},
 score_{0.0}, score_computed_{false}, hash_{0}, hash_computed_{false} {
   
   if(selected_nodes_.empty())
     Rcpp::stop("Cannot create a solution from an empty vector");
   
   std::sort(selected_nodes_.begin(), selected_nodes_.end());
 }

Solution::Solution(std::initializer_list<int> nodes)
  : Solution(std::vector<int>(nodes)) {}

Solution Solution::create_random_valid(const tree_structure& tree,
                                       std::mt19937& rng,
                                       size_t target_size,
                                       int max_attempts,
                                       bool exact_size){
  
  const auto& depths = tree.get_depth();
  
  if (depths.empty()) {
    Rcpp::stop("Cannot create a solution, tree is empty");
  }

  std::uniform_int_distribution<int> dist(0, depths.size() - 1);
  std::vector<int> nodes; nodes.reserve(target_size);
  
  for (int attempt = 0; attempt < max_attempts; ++attempt) {
    
    size_t num_nodes;
    if (exact_size) {
      num_nodes = target_size;
    } else {
      std::uniform_int_distribution<size_t> size_dist(1, target_size);
      num_nodes = size_dist(rng);
    } 
    
    for (size_t i = 0; i < num_nodes; ++i) {
      nodes.push_back(dist(rng));
    }
    
    nodes.shrink_to_fit();
    Solution candidate(nodes);  
    
    if (candidate.is_valid(tree)) {
      return candidate;
    }
    
    nodes.clear();
  }
  // Failed to find valid solution
  Rcpp::warning("Could not generate valid solution in max_attempts. Returning single node.");
  return Solution{std::vector<int>{dist(rng)}};
}

bool Solution::is_valid(const tree_structure& tree) const{
  if(selected_nodes_.empty())
    return false;
  
  // Check for duplicates
  for (size_t i = 1; i < selected_nodes_.size(); ++i) {
    if (selected_nodes_[i] == selected_nodes_[i-1]) {
      return false; 
    }
  }
  
  const auto& upper_bound = tree.get_upper_bound();
  
  for(size_t i = 0; i < selected_nodes_.size()-1; ++i){
    int node_a = selected_nodes_[i];
    for(size_t j = i+1; j < selected_nodes_.size() ; ++j){
      int node_b = selected_nodes_[j];
      
      //node b below node a, l.e. and g.e. to check equality
      if(node_b >= node_a && node_b <= upper_bound[node_a])
        return false;
      
      //node a below node b, l.e. and g.e. to check equality
      if(node_a >= node_b && node_a <= upper_bound[node_b])
        return false;
    }
  }
  
  return true;
}

std::vector<std::pair<int,int>> Solution::determine_vertex(
    const tree_structure& tree) const{
  std::vector<std::pair<int,int>> vertex;
  const auto& depth = tree.get_depth();
  const auto& upper_bound = tree.get_upper_bound();
  
  for(int node : selected_nodes_){
    //if we are below a root node
    if(depth[node] > 1){
      int idx = node-1;
      //we find the father
      while(depth[idx] != depth[node] - 1){
        --idx;
      }
      vertex.push_back(std::make_pair(node,idx));
    }
    
    //if this node have son we add all of his sons
    if(node != upper_bound[node]){
      for(int i = node+1; i <= upper_bound[node]; ++i){
        if(depth[i] == depth[node]+1){
          vertex.push_back(std::make_pair(node,i));
        }
      }
    }
  }
  return vertex;
}

Solution Solution::mutate_swap_type2(const tree_structure& tree,
                                     std::mt19937& rng) const{
  if(empty()){
    return *this;
  }
  auto vertex = determine_vertex(tree);
  
  std::uniform_int_distribution<size_t> vertex_range(0, vertex.size()-1);
  size_t chosen_index = vertex_range(rng);
  
  std::vector<int> new_nodes;
  new_nodes.reserve(selected_nodes_.size());
  
  int node_to_remove = vertex[chosen_index].first;
  int node_to_add = vertex[chosen_index].second;
  
  for(int node : selected_nodes_){
    if(node != node_to_remove)
      new_nodes.push_back(node);
  }
  
  new_nodes.push_back(node_to_add);
  
  return Solution{std::move(new_nodes)};
}

Solution Solution::mutate_add_remove_type1(const tree_structure& tree, 
                                           double alpha,
                                           std::mt19937& rng) const{
  std::uniform_real_distribution<double> runif(0.0,1.0);
  double prob = runif(rng);
  
  auto new_nodes = selected_nodes_;
  //we add if size of solution is 1 or with probability alpha / solution_size
  if(selected_nodes_.size() <= 1 || prob <= alpha / selected_nodes_.size()){
    std::uniform_int_distribution<int> nodes_sampler(0, tree.get_depth().size()-1);
    int new_node = -1;
    int max_attempts = 100;
    for (int attempt = 0; attempt < max_attempts; ++attempt) {
      int candidate = nodes_sampler(rng);
      if(std::find(selected_nodes_.begin(), selected_nodes_.end(),
                   candidate) == selected_nodes_.end()){
       new_node = candidate;
       break; 
      }
    }
    if(new_node == -1){
      Rcpp::warning("In mutation 1, couldn't add a candidate to the solution in "
                      + std::to_string(max_attempts) + "attempts");
      return *this;
    }
    new_nodes.push_back(new_node);
  }else{
    //we remove with probability 1- (alpha/p)
    std::uniform_int_distribution<size_t> index_sampler(0, selected_nodes_.size()-1);
    size_t index_to_remove = index_sampler(rng);
    new_nodes.erase(new_nodes.begin() + index_to_remove);
  }
  return Solution(std::move(new_nodes));
}

Solution Solution::mutate_replace_type1(const tree_structure& tree,
                                        std::mt19937& rng) const{
  return create_random_valid(tree, rng, selected_nodes_.size());
}

Solution Solution::mutate_genetic_algorithm(const tree_structure& tree,
                                        double alpha, std::mt19937& rng) const{
  double type_probability = R::runif(0.0,1.0);
  
  if(type_probability < .5)
    return mutate_add_remove_type1(tree, alpha, rng);
  else
    return mutate_swap_type2(tree, rng);
  
}

std::pair<Solution, Solution> Solution::crossover_single_point(
    const Solution& parent1,
    const Solution& parent2,
    const tree_structure& tree,
    std::mt19937& rng){
  const auto& depths = tree.get_depth();
  const auto& upper_bound = tree.get_upper_bound();
  
  std::vector<int> non_leaf_nodes;
  non_leaf_nodes.reserve(depths.size());
  
  for(size_t i = 0 ; i < depths.size(); ++i){
    if(i != upper_bound[i])
      non_leaf_nodes.push_back(i);
  }
  non_leaf_nodes.shrink_to_fit();
  
  //if no internal node in the tree
  if(non_leaf_nodes.empty())
    return {parent1, parent2};
  
  std::uniform_int_distribution<int> 
    internal_node_sampler(0, non_leaf_nodes.size()-1);
  auto selected_node = non_leaf_nodes[internal_node_sampler(rng)];
  
  std::vector<int> selected_node_sol_1, selected_node_sol_2;
  for(auto node : parent1.get_nodes()){
    if(node >= selected_node && node <= upper_bound[selected_node]){
      selected_node_sol_2.push_back(node);
    }else{
      selected_node_sol_1.push_back(node);
    }
  }
  
  for(auto node : parent2.get_nodes()){
    if(node >= selected_node && node <= upper_bound[selected_node]){
      selected_node_sol_1.push_back(node);
    }else{
      selected_node_sol_2.push_back(node);
    }
  }
  
  if(selected_node_sol_1.empty()){
    selected_node_sol_1.push_back(selected_node);
  }
  if(selected_node_sol_2.empty()){
    selected_node_sol_2.push_back(selected_node);
  }
  
  return std::make_pair(Solution(std::move(selected_node_sol_1)),
                        Solution(std::move(selected_node_sol_2)));
}

void Solution::print() const {
  Rcpp::Rcout << "Solution(nodes=[";
  for (size_t i = 0; i < selected_nodes_.size(); ++i) {
    Rcpp::Rcout << selected_nodes_[i];
    if (i < selected_nodes_.size() - 1) {
      Rcpp::Rcout << ", ";
    }
  }
  Rcpp::Rcout << "], score=";
  if (score_computed_) {
    Rcpp::Rcout << score_;
  } else {
    Rcpp::Rcout << "NOT_COMPUTED";
  }
  Rcpp::Rcout << ", hash=";
  if (hash_computed_) {
    Rcpp::Rcout << hash_;
  } else {
    Rcpp::Rcout << "NOT_COMPUTED";
  }
  Rcpp::Rcout << ")\n";
}
