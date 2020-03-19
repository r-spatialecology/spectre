#ifndef CALCULATE_SOLUTION_COMMONNESS_H
#define CALCULATE_SOLUTION_COMMONNESS_H
#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix calculate_solution_commonness_rcpp(const NumericMatrix solution_matrix);

IntegerMatrix calculate_solution_commonness_rcpp(std::vector<double> solution_matrix);

// [[Rcpp::export]]
NumericMatrix calculate_solution_commonness_site_rcpp(const NumericMatrix solution_matrix,
                                                      const NumericMatrix solution_commonness,
                                                      const int site);
// [[Rcpp::export]]
void update_solution_commonness_site_rcpp(const IntegerMatrix solution_matrix,
                                          IntegerMatrix &solution_commonness,
                                          const unsigned site);
#endif // CALCULATE_SOLUTION_COMMONNESS_H
