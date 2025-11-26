#ifndef TREE_STRUCTURE
#define TREE_STRUCTURE

#include <Rcpp.h>
#include <optional>
#include <string>
#include <vector>

// [[Rcpp::plugins(cpp17)]]
class tree_structure{
private:
  std::vector<int> depth_;
  std::vector<int> upper_bound_;
  std::vector<std::string> name_;
  std::vector<int> father_;
  bool has_father_;
  bool has_name_;
  int max_depth_; 
public:
  tree_structure() = delete;
  tree_structure(const Rcpp::DataFrame& Tree, SEXP depth, 
                 SEXP upper_bound = R_NilValue,
                 SEXP name = R_NilValue);
  tree_structure(const Rcpp::IntegerVector& depth);
  
  inline const std::vector<int>& get_depth() const{
    return depth_;
  }
  
  inline const std::vector<int>& get_upper_bound() const{
    return upper_bound_;
  }
  
  inline std::vector<std::string> get_name() const{
    if(!has_name_)
      Rcpp::stop("Tree has no name column.");
    
    return name_;
  }
  
  inline const std::vector<int>& get_father() const{
    if(has_father_)
      return father_;
    else
      Rcpp::stop("Father of the tree has not been computed.");
  }
  
  inline bool has_name() const{
    return has_name_;
  }
  
  inline int max_depth() const{
    return max_depth_;
  }
  
  void initialize_upper_bound();
  void add_names(const Rcpp::CharacterVector& names);
  void compute_father();
  
private:
  
  template <typename T>
  T get_column(const Rcpp::DataFrame& df, SEXP col){
    if(TYPEOF(col) == STRSXP){
      return df[Rcpp::as<Rcpp::String>(col)];
    } else if (TYPEOF(col) == INTSXP || TYPEOF(col) == REALSXP){
      return df[Rcpp::as<int>(col)-1]; // minus one since we expect user to give R index
    } 
    Rcpp::stop("Columns must be index or name of column");
  }
  
  bool check_depth() const;
  
};

#endif