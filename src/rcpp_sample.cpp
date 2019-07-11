#include "rcpp_sample.h"

Rcpp::IntegerVector rcpp_sample(Rcpp::IntegerVector x, int size, bool replace) {

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
