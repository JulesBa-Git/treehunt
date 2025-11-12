#include <Rcpp.h>
#include "tree_structure.h"

//' @keywords internal
// [[Rcpp::export(.test_create_tree_constructor1)]]
SEXP test_create_tree_constructor1(const Rcpp::DataFrame& df, SEXP depth, 
                      SEXP upper_bound = R_NilValue,
                      SEXP name = R_NilValue){
  tree_structure *tree = new tree_structure(df, depth, upper_bound, name);
  
  return Rcpp::XPtr<tree_structure>(tree, true);
}

//' @keywords internal
 // [[Rcpp::export(.test_create_tree_constructor2)]]
 SEXP test_create_tree_constructor2(const Rcpp::IntegerVector& vect){
   tree_structure *tree = new tree_structure(vect);
   
   return Rcpp::XPtr<tree_structure>(tree, true);
 }

//' @keywords internal
// [[Rcpp::export(.test_tree_info)]]
Rcpp::List test_tree_info(SEXP ptr){
  Rcpp::XPtr<tree_structure> tree(ptr);
  
  Rcpp::List tree_list = Rcpp::List::create(
    Rcpp::Named("max_depth") = tree->max_depth(),
    Rcpp::Named("has_name") = tree->has_name(),
    Rcpp::Named("n_nodes") = tree->get_depth().size(),
    Rcpp::Named("depth") = tree->get_depth(),
    Rcpp::Named("upper_bound") = tree->get_upper_bound()
  );
  
  if(tree->has_name())
    tree_list["name"] = tree->get_name();
  
  return tree_list;
}