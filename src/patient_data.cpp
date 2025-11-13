#include "patient_data.h"

template<typename TargetType>
std::vector<int> PatientData<TargetType>::parse_node_string(
    const std::string& node_str) const {
  
  std::vector<int> nodes;
  std::stringstream ss(node_str);
  std::string item;
  
  while (std::getline(ss, item, ',')) {
    item.erase(0, item.find_first_not_of(" \t\r\n"));
    item.erase(item.find_last_not_of(" \t\r\n") + 1);
    
    if (!item.empty()) {
      try {
        nodes.push_back(std::stoi(item));
      } catch (const std::exception& e) {
        Rcpp::warning("Could not parse node value: " + item);
      }
    }
  }
  
  return nodes;
}

template <typename TargetType>
std::vector<TargetType> PatientData<TargetType>::extract_column(
    const Rcpp::DataFrame& df, SEXP col_target) const{
  SEXP data_col_sexp;
  
  if(TYPEOF(col_target) == STRSXP){
    data_col_sexp = df[Rcpp::as<Rcpp::String>(col_target)];
  }else if(TYPEOF(col_target) == INTSXP || TYPEOF(col_target) == REALSXP){
    data_col_sexp = df[Rcpp::as<int>(col_target)-1];
  }else{
    Rcpp::stop("Column must be specified as a string or integer");
  }
  
  if constexpr (std::is_same_v<TargetType, int>) {
    Rcpp::IntegerVector col_data(data_col_sexp);
    return Rcpp::as<std::vector<int>>(col_data);
  } else if constexpr (std::is_same_v<TargetType, double>) {
    Rcpp::NumericVector col_data(data_col_sexp);
    return Rcpp::as<std::vector<double>>(col_data);
  } else {
    Rcpp::stop("Unsupported target type. Use int (bool) or double.");
  } 
}

template <typename TargetType>
PatientData<TargetType>::PatientData(const Rcpp::DataFrame& df, SEXP node_column, 
                                     SEXP target_column, const tree_structure& tree):
  n_patients_{0}, tree_{tree} {
  
  if(df.nrows() == 0)
    Rcpp::stop("DataFrame is empty");
  
  n_patients_ = df.nrows();
  
  targets_ = extract_column(df, target_column);
  
  SEXP node_col_data;
  if (TYPEOF(node_column) == STRSXP) {
    node_col_data = df[Rcpp::as<Rcpp::String>(node_column)];
  } else {
    node_col_data = df[Rcpp::as<int>(node_column)-1];
  } 
  
  patient_nodes_.reserve(n_patients_);
  
  //we handle 2 format, list of integer vector and vector of string
  if(TYPEOF(node_col_data) == VECSXP){
    Rcpp::List nodes_list(node_col_data);
    
    for(int i = 0; i < nodes_list.size(); ++i){
      if(TYPEOF(nodes_list[i]) == INTSXP || TYPEOF(nodes_list[i]) == REALSXP){
        patient_nodes_.push_back(Rcpp::as<std::vector<int>>(nodes_list[i]));
      }else{
        Rcpp::stop("List elements must be integer or numeric vectors at row " + 
          std::to_string(i + 1));
      }
  
    }
  } else if (TYPEOF(node_col_data) == STRSXP){
    Rcpp::CharacterVector nodes_vec(node_col_data);
    
    for (int i = 0; i < nodes_vec.size(); ++i) {
      std::string node_str = Rcpp::as<std::string>(nodes_vec[i]);
      patient_nodes_.push_back(parse_node_string(node_str));
    }
  }else{
    Rcpp::stop("Node column must be one of:\n"
                 "  1. List of integer vectors: list(c(1,2), c(3), c(4,5))\n"
                 "  2. Character vector with comma-separated values: c('1,2', '3', '4,5')\n");
  }
  
  if (targets_.size() != n_patients_) {
    Rcpp::stop("Target column size doesn't match number of patients");
  }
  
}

template class PatientData<int>;
template class PatientData<double>;