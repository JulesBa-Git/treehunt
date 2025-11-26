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
  size_t n_patients_;
  const tree_structure& tree_;
  
  std::vector<int> parse_node_string(const std::string& node_str) const;
  
  std::vector<TargetType> extract_column(const Rcpp::DataFrame& df,
                                         SEXP col_target) const;
  
public:
  PatientData() = delete;
  
  PatientData(const Rcpp::DataFrame& df, SEXP node_column, SEXP target_column,
              const tree_structure& tree);
  
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
    
    for(int patient_node : patient_nodes){
      bool is_covered = false;
      
      for(int solution_node : solution_nodes){
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
  
  std::pair<std::vector<std::vector<int>>, std::vector<int>> 
    dissimilarity_input_computation() const {
      std::vector<std::vector<int>> M;
      std::vector<int> index;
      const auto& father = tree_.get_father();
      const auto& depth = tree_.get_depth();
      
      int node_count = 0;
      for(const auto& nodes : patient_nodes_){
        node_count+= nodes.size();
      }
      index.reserve(node_count);
      M.resize(node_count);
      for(size_t i = 0; i < M.size(); ++i)
        M[i].reserve(tree_.max_depth() + 1);
      
      size_t row = 0;
      size_t current_patient_node_index = 0;
      for(const std::vector<int>& nodes_vector : patient_nodes_){
        for(const auto& node : nodes_vector){
          //add the node
          M[row].push_back(node);
          int current_father = father[node];
          //add the father of this node
          while(depth[current_father] != 1){
            M[row].push_back(current_father);
            current_father = father[current_father];
          }// push the father of depth 1
          M[row].push_back(current_father);
        }
        index.push_back(current_patient_node_index);
        current_patient_node_index += nodes_vector.size();
      }
      
      return std::make_pair(std::move(M), std::move(index));
    }
};

#endif