#ifndef PATIENT_DATA_H
#define PATIENT_DATA_H

#include <Rcpp.h>
#include <vector>
#include "tree_structure.h"

template <typename TargetType = int>
class PatientData{
private :
  std::vector<std::vector<int>> patient_nodes_;
  std::vector<TargetType> targets_;
  // Note: Patients id are only integer for now (should cover most of application),
  // could be a template ?
  std::vector<int> patients_id_; 
  size_t n_patients_;
  const tree_structure& tree_;
  
  std::vector<int> parse_node_string(const std::string& node_str) const;
  
  template <typename T>
  std::vector<T> extract_column(const Rcpp::DataFrame& df,
                                SEXP col_target) const;
  
public:
  PatientData() = delete;
  
  PatientData(const Rcpp::DataFrame& df, SEXP node_column, SEXP target_column,
              const tree_structure& tree, SEXP id_column = R_NilValue);
  
  const std::vector<int>& get_patient_nodes(size_t i) const{
    if (i >= n_patients_) {
      Rcpp::stop("Patient index out of bounds");
    }
    return patient_nodes_[i];
  }
  
  TargetType get_target(size_t i) const{
    if( i >= n_patients_)
      Rcpp::stop("Patient index out of bounds");
    return targets_[i];
  }
  
  int get_id(size_t i) const{
    if( i >= n_patients_)
      Rcpp::stop("Patient index out of bounds");
    return patients_id_[i];
  }
  
  void set_ids(const Rcpp::DataFrame& df, SEXP id_column){
    patients_id_ = extract_column<int>(df, id_column);
  }
  
  size_t size() const { return n_patients_; }
  
  const tree_structure& get_tree() const {
    return tree_;
  }
  
  double mean_patient_nodes() const{
    double total_size = 0.0;
    for(const auto& node : patient_nodes_){
      total_size += node.size();
    }
    return total_size / static_cast<double>(patient_nodes_.size());
  }
  
  bool patient_has_combination(size_t i,
                               const std::vector<int>& solution_nodes) const{
    if(i >= n_patients_)
      return false;
    
    const auto& patient_nodes = patient_nodes_[i];
    const auto& upper_bound = tree_.get_upper_bound();
    
    //if the solution is larger than the patient intake, it can not contain it
    if(patient_nodes.size() < solution_nodes.size())
      return false;
    
    for(int solution_node : solution_nodes){
      bool is_covered = false;
      
      for(int patient_node : patient_nodes){
        if(patient_node >= solution_node &&
           patient_node <= upper_bound[solution_node]){
          is_covered = true;
          break;
        }
      }
      if(!is_covered)
        return false;
    }
    return true;
  }
};

#endif