#include <Rcpp.h>
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
Rcpp::NumericVector rcpp_sample(Rcpp::NumericVector x, int size, bool replace = false) {
  
  return Rcpp::sample(x, size, replace);
  
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
*/