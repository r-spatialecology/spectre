#include <Rcpp.h>
#include "rcpp_sample.h"
using namespace Rcpp;

//' rcpp_sample
//'
//' @description Rcpp sample function
//'
//' @param x Vector of elements to sample from.
//' @param size Size of the sample.
//' @param replace Sample with replacement.
//'
//' @details
//' \code{Rcpp} implementation of the \code{sample} function.
//'
//' @seealso
//' \code{\link{sample}}
//'
//' @return vector
//'
//' @name rcpp_sample
//' @export
// [[Rcpp::export]]
IntegerVector rcpp_sample(IntegerVector x, int size, bool replace = false) {

  return sample(x, size, replace);

}

// [[Rcpp::export]]
IntegerVector which_not(IntegerVector x, int y) {
  IntegerVector v = seq(0, x.size()-1);
  return v[x != y];
}

//// [[Rcpp::export]]
//IntegerVector get_swap_rows_rcpp(const IntegerVector v) {
//  RNGScope scope;
//  IntegerVector retval(2);
//  int n_rows = v.length();
//  int species1 = R::runif(1, 140);
//  for (int i = 0; i < 100; i++) {
//      std::cout << R::runif(1, 140) << std::endl;
//  }
//  IntegerVector xxx = sample(n_rows, 139);
//  auto x = as<std::vector<int> >(xxx);
//  auto yy = runif(139, 1, 140);
//  auto y = as<std::vector<int> >(yy);
//  int species2 = rcpp_sample(which_not(v, species1), 10)[0];
//  retval[0] = species1; retval[1] = species2;
//  assert(v[species1] != v[species2]);
//  return(retval);
//}

// [[Rcpp::export]]
IntegerVector get_swap_rows_rcpp(const IntegerVector v) {
  RNGScope scope; ///TODO: understand why this is needed
  unsigned species1 = sample(IntegerVector(seq_len(v.length())), 1)[0];
  unsigned species2 = sample(which_not(v, v[species1]), 1)[0];
  return(IntegerVector::create(species1 + 1, species2 + 1)); // R starts counting at 1
}

// [[Rcpp::export]]
IntegerVector get_swap_rows_rcpp_bruteforce(const IntegerVector v) {
  RNGScope scope; ///TODO: understand why this is needed
  IntegerVector rows(seq_len(v.length()));
  unsigned species1 = sample(rows, 1)[0];
  unsigned species2 = sample(rows, 1)[0];

  while ((v[species1] + v[species2]) != 1) {
    species2 = sample(rows, 1)[0];
  }
  return(IntegerVector::create(species1 + 1, species2 + 1)); // R starts counting at 1
}

/*** R
rcpp_sample(x = 1:100, size = 5)
rcpp_sample(x = 1:5, size = 10, replace = TRUE)

bench::mark(rcpp_sample(x = 1:100, size = 5),
            sample(x = 1:100, size = 5),
            check = FALSE, relative = TRUE, iterations = 10000)

bench::mark(rcpp_sample(x = 1:10, size = 20, replace = TRUE),
            sample(x = 1:10, size = 20, replace = TRUE),
            check = FALSE, relative = TRUE, iterations = 10000)

mat <- matrix(sample(0:1, 480000, replace = TRUE), 6000, 8000)

x <- spectre:::get_swap_rows(mat[, site])
mat[x[1], site] != mat[x[2], site]
x <- spectre:::get_swap_rows_rcpp(mat[, site])
mat[x[1], site] != mat[x[2], site]
x <- spectre:::get_swap_rows_rcpp_bruteforce(mat[, site])
mat[x[1], site] != mat[x[2], site]

bench::mark(spectre:::get_swap_rows(mat[, site]),
            spectre:::get_swap_rows_rcpp(mat[, site]),
            spectre:::get_swap_rows_rcpp_bruteforce(mat[, site]),
            check = FALSE)
# 1 spectre:::get_swap_rows(mat[, site])                 321µs 652.8µs      909.   134.2KB        0   100     0
# 2 spectre:::get_swap_rows_rcpp(mat[, site])            162µs 270.1µs     2942.   143.4KB        0   100     0
# 3 spectre:::get_swap_rows_rcpp_bruteforce(mat[, site])  63µs  86.8µs     9964.    72.9KB        0   100     0
*/
