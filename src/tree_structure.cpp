#include "tree_structure.h"

tree_structure::tree_structure(const Rcpp::DataFrame &Tree, SEXP depth,
                               SEXP upper_bound, SEXP name) {

  depth_ =
      Rcpp::as<std::vector<int>>(get_column<Rcpp::IntegerVector>(Tree, depth));
  check_depth();

  if (depth_.empty()) {
    Rcpp::stop(
        "No value in depth vector, the given column index might be wrong");
  } else {
    max_depth_ = *std::max_element(depth_.begin(), depth_.end());
  }

  if (upper_bound != R_NilValue) {
    upper_bound_ = Rcpp::as<std::vector<int>>(
        get_column<Rcpp::IntegerVector>(Tree, upper_bound));
  } else {
    initialize_upper_bound();
  }

  if (name != R_NilValue) {
    name_ = Rcpp::as<std::vector<std::string>>(
        get_column<Rcpp::StringVector>(Tree, name));
    has_name_ = true;
  } else {
    has_name_ = false;
  }
  compute_father();
  has_father_ = true;
}

tree_structure::tree_structure(const Rcpp::IntegerVector& depth): depth_{},
upper_bound_{}, name_{}, father_{}, has_father_{false},has_name_{false},
max_depth_{} {
  depth_ = Rcpp::as<std::vector<int>>(depth);
  
  if (depth_.empty()) {
    Rcpp::stop(
      "No value in depth vector, the given vector is wrong");
  } else {
    max_depth_ = *std::max_element(depth_.begin(), depth_.end());
  }
  
  check_depth();
  initialize_upper_bound();
}

void tree_structure::initialize_upper_bound() {
  auto current_father = std::vector<int>(max_depth_ + 1);
  current_father[0] = -1;
  int prev_dep = 0;

  std::vector<int> tmp_upper_bound(depth_.size());

  for (size_t i = 0; i < tmp_upper_bound.size() - 1; ++i) {

    if (depth_[i] < prev_dep) {
      for (int j = depth_[i]; j < prev_dep; ++j) {
        tmp_upper_bound[current_father[j]] = current_father[prev_dep];
      }
    }

    if (depth_[i] >= depth_[i + 1]) {
      tmp_upper_bound[i] = i;
    }
    prev_dep = depth_[i];
    current_father[depth_[i]] = i;
  }
  
  // Soit, je suis un noeud de depth inférieur à prev_dep, dans ce cas je mets à 
  // jour les upper_bound de (depth à prev_depth) = tmp_upper_bound.size()-2 
  // puis je mets à jours les noeuds au dessus de moi et moi à tmp_upper_bound.size()-1
  // Soit je suis sur un noeud de depth >= prev_depth et je mets à jours tous les noeuds
  // au dessus de moi jusqu'à moi à tmp_upper_bound.size()-1
  auto current_depth = depth_[depth_.size()-1];
  if(current_depth < prev_dep){
    for(size_t i = current_depth; i <= prev_dep; ++i){
      tmp_upper_bound[current_father[i]] = tmp_upper_bound.size()-2;
    }
  }
  for(size_t i = 1; i < current_depth; ++i){
    tmp_upper_bound[current_father[i]] = tmp_upper_bound.size()-1;
  }
  
  tmp_upper_bound[depth_.size() - 1] = depth_.size() - 1;

  upper_bound_ = std::move(tmp_upper_bound);
}

void tree_structure::add_names(const Rcpp::CharacterVector& names){
  if(names.size() != depth_.size())
    Rcpp::stop("There must be as many name nodes as there are nodes in the tree.");
  
  name_ = Rcpp::as<std::vector<std::string>>(names);
  has_name_ = true;
}

void tree_structure::compute_father(){
  if(!has_father_){
    std::vector<int>current_father(max_depth_ + 1);
    current_father[0] = -1;
    
    std::vector<int> father; father.reserve(depth_.size());
    for(size_t i = 0 ; i < depth_.size(); ++i){
      father.push_back(current_father[depth_[i] - 1]);
      current_father[depth_[i]] = i;
    }
    has_father_ = true;
  }
}

bool tree_structure::check_depth() const {
  int prev_dep = 0;
  for (const auto dep : depth_) {
    // if we go on a son, his depth should be strictly one more than our
    if (dep > prev_dep && dep != prev_dep + 1) {
      Rcpp::stop("The son of a node must be in the next depth of the tree.");
    }
    prev_dep = dep;
  }
  return true;
}
