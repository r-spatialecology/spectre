#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerMatrix calculate_solution_commonness_site_rcpp(IntegerMatrix solution_matrix, 
                                                 IntegerMatrix solution_commonness, 
                                                 int site) {

  int nrows = solution_commonness.nrow();
  site = site - 1; // because C++ starts indexing at zero
  
  int j = 0;
  
  for (j = 0; j < site; j++) {
    solution_commonness(site, j) = sum((solution_matrix(_, j) + solution_matrix(_, site)) > 1);
  }
  
  j++; // just to avoid comparing with itself
  
  for (; j < nrows; j++) {
    solution_commonness(j, site) = sum((solution_matrix(_, j) + solution_matrix(_, site)) > 1);
  }
  
  return solution_commonness ;

}


/*** R
dim(new_solution)

new_solution <- matrix(sample(0:1, 48, replace = TRUE), 6, 8)

solution_commonness <- matrix(sample(0:8, 64, replace = TRUE), 8, 8)

random_col <- 6

calculate_solution_commonness_site(new_solution, solution_commonness, random_col)

calculate_solution_commonness_site_rcpp(new_solution, solution_commonness, random_col)


bench::mark(
  calculate_solution_commonness_site(new_solution, solution_commonness, random_col),
  
  calculate_solution_commonness_site_rcpp(new_solution, solution_commonness, random_col)
  
  
)
*/