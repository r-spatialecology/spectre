#ifndef CALCULATE_SOLUTION_COMMONNESS_H
#define CALCULATE_SOLUTION_COMMONNESS_H
#include <Rcpp.h>
//#include <RcppParallel.h>

using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
IntegerMatrix calculate_solution_commonness_rcpp(IntegerMatrix solution_matrix);

// [[Rcpp::export]]
IntegerMatrix calculate_solution_commonness_site_rcpp(IntegerMatrix solution_matrix,
                                                      IntegerMatrix solution_commonness,
                                                      int site);
#endif // CALCULATE_SOLUTION_COMMONNESS_H
