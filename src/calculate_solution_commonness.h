#ifndef CALCULATE_SOLUTION_COMMONNESS_H
#define CALCULATE_SOLUTION_COMMONNESS_H
#include <Rcpp.h>
//#include <RcppParallel.h>

using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
IntegerMatrix calculate_solution_commonness_rcpp(IntegerMatrix solution_matrix);

// [[Rcpp::export]]
IntegerMatrix calculate_solution_commonness_rcpp_p(IntegerMatrix solution_matrix);

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
IntegerMatrix calculate_solution_commonness_rcpp_old(IntegerMatrix solution_matrix);
#include <omp.h>
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
IntegerMatrix calculate_solution_commonness_site_rcpp(const IntegerMatrix solution_matrix,
                                                      IntegerMatrix solution_commonness,
                                                      const int site);
// [[Rcpp::export]]
IntegerMatrix calculate_solution_commonness_site_rcpp_p(const IntegerMatrix solution_matrix,
                                                      IntegerMatrix solution_commonness,
                                                      const int site);
#endif // CALCULATE_SOLUTION_COMMONNESS_H
