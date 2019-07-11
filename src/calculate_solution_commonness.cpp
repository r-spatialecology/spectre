#include "calculate_solution_commonness.h"

// approx 2x faster than old variant
// [[Rcpp::plugins(cpp11)]]
IntegerMatrix calculate_solution_commonness_rcpp(const IntegerMatrix solution_matrix) {

    // create results matrix, filled with NAs
    int nsites = solution_matrix.ncol();
    IntegerMatrix commonness_mat(nsites, nsites);
    std::fill(commonness_mat.begin(), commonness_mat.end(), NumericVector::get_na() ) ;

    // iterate over the lower triangular matrix without diagonal
    // and calculate # of common species between sites
    for (int col = 0; col < nsites - 1; col++) {
        for (int row = col + 1; row < nsites; row++) {
            commonness_mat(row, col) = sum((solution_matrix(_, row) +
                                            solution_matrix(_, col)) > 1);
        }
    }

    return commonness_mat ;

}

IntegerMatrix calculate_solution_commonness_rcpp_old(IntegerMatrix solution_matrix) {

    int ncols = solution_matrix.ncol();
    IntegerMatrix commonness_s_mat(ncols, ncols);

    for (int i = 0; i < ncols; i++) {
        for (int j = 0; j < ncols; j++) {
            commonness_s_mat(j,i) = sum((solution_matrix(_, j) +
                                         solution_matrix(_, i)) > 1);
        }
    }

    return commonness_s_mat ;

}


IntegerMatrix calculate_solution_commonness_site_rcpp(const IntegerMatrix solution_matrix,
                                                      IntegerMatrix solution_commonness,
                                                      int site) {

    const int nrows = solution_commonness.nrow();
    const int site_ = site - 1; // because C++ starts indexing at zero


    for (int j = 0; j < site_; j++) {
        solution_commonness(site_, j) = sum((solution_matrix(_, j) +
                                             solution_matrix(_, site_)) > 1);
    }

    for (int j = site_ + 1; j < nrows; j++) {
        solution_commonness(j, site_) = sum((solution_matrix(_, j) +
                                             solution_matrix(_, site_)) > 1);
    }

    return solution_commonness;

}

void update_solution_commonness_site_rcpp(const IntegerMatrix solution_matrix,
                                          IntegerMatrix &solution_commonness,
                                          const int site) {

    const int nrows = solution_commonness.nrow();
    const int site_ = site - 1; // because C++ starts indexing at zero


    for (int j = 0; j < site_; j++) {
        solution_commonness(site_, j) = sum((solution_matrix(_, j) +
                                             solution_matrix(_, site_)) > 1);
    }

    for (int j = site_ + 1; j < nrows; j++) {
        solution_commonness(j, site_) = sum((solution_matrix(_, j) +
                                             solution_matrix(_, site_)) > 1);
    }

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
