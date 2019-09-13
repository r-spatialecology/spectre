#ifndef RCPP_SAMPLE_H
#define RCPP_SAMPLE_H
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
Rcpp::IntegerVector rcpp_sample(Rcpp::IntegerVector x, int size, bool replace = false);

// [[Rcpp::export]]
IntegerVector get_swap_rows_rcpp(const IntegerVector &v, bool cpp_index);

// [[Rcpp::export]]
IntegerVector get_swap_rows_rcpp_bruteforce(const IntegerVector &v, bool cpp_index);

// [[Rcpp::export]]
IntegerVector which_not_vec(const IntegerVector x, const IntegerVector y);

// [[Rcpp::export]]
IntegerVector which_not(const IntegerVector x, const int y);
#endif // RCPP_SAMPLE_H
