#include "helpers.h"

double normalized_distance(size_t i, size_t j,
                           const std::vector<std::vector<int>>& M,
                           const std::vector<int>& index, 
                           int max_depth){
  size_t begin_M_index_i = index[i], begin_M_index_j = index[j];
  size_t end_M_index_i = index[i+1];
  size_t end_M_index_j = j == index.size() -1 ? M.size() : index[j+1];
  
  std::unordered_set<int> index_i;
  for(int k = begin_M_index_i; k < end_M_index_i; ++k )
    index_i.insert(k);
  
  std::unordered_set<int> index_j;
  for(int k = begin_M_index_j; k < end_M_index_j; ++k)
    index_j.insert(k);
  double cost = 0;
  
  return 0.0;
}