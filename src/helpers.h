#ifndef HELPERS_TREEHUNT_H
#define HELPERS_TREEHUNT_H

#include <numeric>
#include <vector>
#include <unordered_set>

double normalized_distance(size_t i, size_t j,
                           const std::vector<std::vector<int>>& M,
                           const std::vector<int>& index, 
                           int max_depth);

#endif