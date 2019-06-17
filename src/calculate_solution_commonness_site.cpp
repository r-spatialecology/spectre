#include <Rcpp.h>
#include <RcppParallel.h>

using namespace Rcpp;

// [[Rcpp::export]]
IntegerMatrix calculate_solution_commonness_rcpp(IntegerMatrix solution_matrix) {
  
  int ncols = solution_matrix.ncol();
  IntegerMatrix commonness_s_mat(ncols, ncols, NA) ;
  
  
  for (int i = 0; i < ncols; i++) {
    for (int j = 0; j < ncols; j++) {
      commonness_s_mat(j,i) = sum((solution_matrix(_, j) + solution_matrix(_, i)) > 1);
  }

  return commonness_s_mat ;
  
}

// [[Rcpp::export]]
IntegerMatrix calculate_solution_commonness_site_rcpp(IntegerMatrix solution_matrix, 
                                                 IntegerMatrix solution_commonness, 
                                                 int site) {

  int nrows = solution_commonness.nrow();
  site = site - 1; // because C++ starts indexing at zero
  

  for (int j = 0; j < site; j++) {
    solution_commonness(site, j) = sum((solution_matrix(_, j) + solution_matrix(_, site)) > 1);
  }
  
  for (int j = site; j < nrows; j++) {
    solution_commonness(j, site) = sum((solution_matrix(_, j) + solution_matrix(_, site)) > 1);
  }
  
  return solution_commonness ;

}

// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::export]]
IntegerMatrix calculate_solution_commonness_site_rcpp_par(IntegerMatrix solution_matrix, 
                                                      IntegerMatrix solution_commonness, 
                                                      int site) {
  
  int nrows = solution_commonness.nrow();
  site = site - 1; // because C++ starts indexing at zero
  
  
  for (int j = 0; j < site; j++) {
    solution_commonness(site, j) = sum((solution_matrix(_, j) + solution_matrix(_, site)) > 1);
  }
  
  for (int j = site; j < nrows; j++) {
    solution_commonness(j, site) = sum((solution_matrix(_, j) + solution_matrix(_, site)) > 1);
  }
  
  return solution_commonness ;
  
}
/*** R
dim(new_solution)

new_solution <- matrix(sample(0:1, 480000, replace = TRUE), 6000, 8000)

solution_commonness <- matrix(sample(0:8, 640000, replace = TRUE), 8000, 8000)

random_col <- 68

calculate_solution_commonness_site <- function(solution_matrix, solution_commonness, site){
  
  for(j in 1:nrow(solution_commonness)){
    if (j < site){
      solution_commonness[site,j] <- sum((solution_matrix[,j] + solution_matrix[,site]) > 1)
    }
    if (j > site){
      solution_commonness[j,site] <- sum((solution_matrix[,j] + solution_matrix[,site]) > 1)
    }
    
  }
  return(solution_commonness)
}

calculate_solution_commonness_site_new <- function(solution_matrix, solution_commonness, site){
  
  for (j in 1:site) {
    solution_commonness[site, j] = sum((solution_matrix[, j] + solution_matrix[, site]) > 1)
  }
  
  for (j in site:nrow(solution_commonness)) {
    solution_commonness[j, site] = sum((solution_matrix[, j] + solution_matrix[, site]) > 1)
  }
  
  return(solution_commonness)
}

calculate_solution_commonness_site(new_solution, solution_commonness, random_col) == calculate_solution_commonness_site_rcpp(new_solution, solution_commonness, random_col) 
calculate_solution_commonness_site(new_solution, solution_commonness, random_col) == calculate_solution_commonness_site_new(new_solution, solution_commonness, random_col)

bench::mark(
  calculate_solution_commonness_site(new_solution, solution_commonness, random_col),
  calculate_solution_commonness_site_rcpp(new_solution, solution_commonness, random_col),
  calculate_solution_commonness_site_new(new_solution, solution_commonness, random_col),
  iterations = 100
)

*/