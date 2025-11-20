#ifndef SCORE_FUNCTIONS_H
#define SCORE_FUNCTIONS_H

#include <Rcpp.h>
#include <cmath>
#include <limits>
#include "solution.h"
#include "patient_data.h"

enum class ScoreType{
  HYPERGEOMETRIC,
  RELATIVE_RISK
};

template<typename TargetType>
class ScoreFunctions{
  
public:
  struct ScoreData{
    double score;
    size_t covered_patients;
    size_t covered_nonzero_target;
  };
  //function to use with TargetType being bool (or 0-1 integer)
  /*static ScoreData compute_hypergeometric_with_data(const PatientData<TargetType>& data,
                                     const Solution& solution){
    ScoreData result;
    result.covered_patients = 0;
    result.covered_nonzero_target = 0;
    
    size_t nonzero_target = 0;
    const auto& nodes = solution.get_nodes();
    
    for(size_t i = 0; i < data.size(); ++i){
      bool patient_have_solution = data.patient_has_combination(i, nodes);
      
      if(patient_have_solution){
        ++result.covered_patients; //K
        if(data.get_target(i)){
          ++result.covered_nonzero_target; //Q
        }
      }
      if(data.get_target(i)){
        ++nonzero_target;
      }
    }
    
    size_t zero_target = data.size() - nonzero_target;
    
    if (result.covered_patients == 0) {
      result.score = 0.0; 
      return result;
    }
    
    auto log_phyper = R::phyper(result.covered_nonzero_target-1,
                                nonzero_target, 
                                zero_target,
                                result.covered_patients, 
                                false, 
                                true);
    result.score = -log_phyper;
    
    return result;
  }*/
  
  static ScoreData compute_hypergeometric_with_data(
      const PatientData<TargetType>& data,
      const Solution& solution,
      bool parallel = true) {
    
    ScoreData result;
    result.covered_patients = 0;
    result.covered_nonzero_target = 0;
    
    size_t nonzero_target = 0;
    const auto& nodes = solution.get_nodes();
    
    size_t n = data.size();
    
    bool use_parallel = parallel && n > 1000;  // Only parallelize for large datasets
    
#ifdef _OPENMP
    if (use_parallel) {
      Rcpp::Rcout << "CPU parallelization enabled\n";
      // Parallel version with OpenMP
      size_t local_covered = 0;
      size_t local_covered_nonzero = 0;
      size_t local_nonzero = 0;
      
#pragma omp parallel for reduction(+:local_covered,local_covered_nonzero,local_nonzero) schedule(static)
      for (size_t i = 0; i < n; ++i) {
        bool has_combination = data.patient_has_combination(i, nodes);
        bool has_nonzero_target = (data.get_target(i) != static_cast<TargetType>(0));
        
        if (has_combination) {
          local_covered++;
          if (has_nonzero_target) {
            local_covered_nonzero++;
          }
        }
        
        if (has_nonzero_target) {
          local_nonzero++;
        }
      }
      
      result.covered_patients = local_covered;
      result.covered_nonzero_target = local_covered_nonzero;
      nonzero_target = local_nonzero;
      
    } else {
#endif
      // Serial version (small datasets or OpenMP not available)
      for (size_t i = 0; i < n; ++i) {
        bool has_combination = data.patient_has_combination(i, nodes);
        
        if (has_combination) {
          ++result.covered_patients;
          if (data.get_target(i)) {
            ++result.covered_nonzero_target;
          }
        }
        
        if (data.get_target(i)) {
          ++nonzero_target;
        }
      }
#ifdef _OPENMP
    }
#endif
    
    if (result.covered_patients == 0) {
      result.score = 0.0;
      return result;
    }
    
    size_t zero_target = n - nonzero_target;
    
    auto log_phyper = R::phyper(
      result.covered_nonzero_target - 1,
      nonzero_target,
      zero_target,
      result.covered_patients,
      false,
      true
    );
    
    result.score = -log_phyper;
    
    return result;
  }
  
  static double compute_hypergeometric(const PatientData<TargetType>& data,
                                       const Solution& solution,
                                       bool parallel = true) {
    return compute_hypergeometric_with_data(data, solution, parallel).score;
  }
  
  static ScoreData compute_relative_risk_with_data(const PatientData<TargetType>& data,
                                                    const Solution& solution){
    ScoreData result;
    result.covered_patients = 0;
    result.covered_nonzero_target = 0;
    
    size_t noncovered_nonzero_target = 0;
    size_t noncovered = 0;
    const auto& nodes = solution.get_nodes();
    
    for(size_t i = 0; i < data.size(); ++i){
      bool patient_have_solution = data.patient_has_combination(i, nodes);
      
      if(patient_have_solution){
        ++result.covered_patients;
        if(data.get_target(i)){
          ++result.covered_nonzero_target;
        }
      }else{
        ++noncovered;
        if(data.get_target(i)){
          ++noncovered_nonzero_target;
        }
      }
    }
    
    if (result.covered_patients == 0 || noncovered == 0) {
      result.score = 0.0; 
      return result;
    } 
    
    
    double prop_non_zero_knowing_covered = 
      static_cast<double>(result.covered_nonzero_target) / 
      static_cast<double>(result.covered_patients);
      
      double prop_non_zero_knowing_notcovered = 
        static_cast<double>(noncovered_nonzero_target) / 
        static_cast<double>(noncovered);
    
    if (prop_non_zero_knowing_notcovered == 0.0) {
      if (prop_non_zero_knowing_covered > 0.0) {
        // infinite relative risk
        result.score = std::numeric_limits<double>::infinity();
      } else {
        // both are 0 - undefined ? return 1.0 ?
        result.score = 1.0;
      }
      return result;
    }
    
    result.score = prop_non_zero_knowing_covered / prop_non_zero_knowing_notcovered;
    
    return result;
  }
  
  static double compute_relative_risk(const PatientData<TargetType>& data,
                                      const Solution& solution) {
    return compute_relative_risk_with_data(data, solution).score;
  }
};
#endif