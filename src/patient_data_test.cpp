#include <Rcpp.h>
#include "patient_data.h"
#include "tree_structure.h"

//' @keywords internal
// [[Rcpp::export(.test_create_patient_data_int)]]
 SEXP test_create_patient_data_int(Rcpp::DataFrame df,
                                  SEXP node_column,
                                  SEXP target_column,
                                  SEXP tree_ptr) {
   Rcpp::XPtr<tree_structure> tree(tree_ptr);
   
   PatientData<int>* data = new PatientData<int>(df, node_column, target_column, *tree);
   return Rcpp::XPtr<PatientData<int>>(data, true);
 }

//' @keywords internal
// [[Rcpp::export(.test_create_patient_data_double)]]
 SEXP test_create_patient_data_double(Rcpp::DataFrame df,
                                      SEXP node_column,
                                      SEXP target_column,
                                      SEXP tree_ptr) {
   Rcpp::XPtr<tree_structure> tree(tree_ptr);
   
   PatientData<double>* data = new PatientData<double>(df, node_column, target_column, *tree);
   return Rcpp::XPtr<PatientData<double>>(data, true);
 }

//' @keywords internal
// [[Rcpp::export(.test_patient_data_size)]]
 int test_patient_data_size(SEXP data_ptr) {
   Rcpp::XPtr<PatientData<int>> ptr(data_ptr);
   return ptr->size();
 }

//' @keywords internal
// [[Rcpp::export(.test_patient_data_get_nodes)]]
 Rcpp::IntegerVector test_patient_data_get_nodes(SEXP data_ptr, int patient_idx) {
   Rcpp::XPtr<PatientData<int>> ptr(data_ptr);
   return Rcpp::wrap(ptr->get_patient_nodes(patient_idx));
 }

//' @keywords internal
// [[Rcpp::export(.test_patient_data_get_target)]]
 int test_patient_data_get_target(SEXP data_ptr, int patient_idx) {
   Rcpp::XPtr<PatientData<int>> ptr(data_ptr);
   return ptr->get_target(patient_idx);
 }

//' @keywords internal
// [[Rcpp::export(.test_patient_data_get_target_double)]]
 double test_patient_data_get_target_double(SEXP data_ptr, int patient_idx) {
   Rcpp::XPtr<PatientData<double>> ptr(data_ptr);
   return ptr->get_target(patient_idx);
 }

//' @keywords internal
// [[Rcpp::export(.test_patient_data_has_combination)]]
 bool test_patient_data_has_combination(SEXP data_ptr,
                                        int patient_idx,
                                        Rcpp::IntegerVector combination) {
   Rcpp::XPtr<PatientData<int>> ptr(data_ptr);
   std::vector<int> combo = Rcpp::as<std::vector<int>>(combination);
   return ptr->patient_has_combination(patient_idx, combo);
 }
