#include "helpers.h"

double normalized_distance(size_t i, size_t j,
                           const std::vector<std::vector<int>>& M,
                           const std::vector<int>& index, 
                           int max_depth){
  size_t begin_M_index_i = index[i], begin_M_index_j = index[j];
  size_t end_M_index_i = index[i+1];
  size_t end_M_index_j = j == index.size() -1 ? M.size() : index[j+1];
  double insertion_cost = static_cast<double>(max_depth)/2.0;
  
  std::unordered_set<int> index_i;
  for(int k = begin_M_index_i; k < end_M_index_i; ++k )
    index_i.insert(k);
  std::vector<int> index_i_to_delete;
  
  std::unordered_set<int> index_j;
  for(int k = begin_M_index_j; k < end_M_index_j; ++k)
    index_j.insert(k);
  double cost = 0;
  double initial_length = index_i.size() + index_j.size();
  int depth = 0;
  
  while(!index_i.empty() && !index_j.empty()){
    for(const auto& idx_i : index_i){
      int matched_index_j = -1;
      
      for(const auto& idx_j : index_j){
        if(M[idx_i][depth] == M[idx_j][depth]){
          matched_index_j = idx_j;
          break;
        }
      }
      
      if(matched_index_j != -1){
        if(depth == 0){
          cost += 0;
        }else{
          size_t count_i = std::count(M[idx_i].begin(),M[idx_i].end(), M[idx_i][0]);
          size_t count_j = std::count(M[matched_index_j].begin(),
                                      M[matched_index_j].end(),
                                      M[matched_index_j][0]);
          cost += static_cast<double>(depth) - std::min(count_i, count_j) + 1.0;
        }
        index_j.erase(matched_index_j);
        index_i_to_delete.push_back(idx_i);
      }
    }
    for(int idx_i : index_i_to_delete)
      index_i.erase(idx_i);
    ++depth;
  }
  
  double additional_cost = (index_j.size() + index_i.size()) * insertion_cost;
  
  return (cost + additional_cost) / (static_cast<double>(initial_length) * insertion_cost);
}