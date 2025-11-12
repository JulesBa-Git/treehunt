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
