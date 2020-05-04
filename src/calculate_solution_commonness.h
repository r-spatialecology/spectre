#ifndef CALCULATE_SOLUTION_COMMONNESS_H
#define CALCULATE_SOLUTION_COMMONNESS_H
#include <Rcpp.h>

using namespace Rcpp;

//NumericMatrix calculate_solution_commonness_rcpp(const NumericMatrix solution_matrix);

// [[Rcpp::export]]
IntegerMatrix calculate_solution_commonness_rcpp(IntegerMatrix solution_matrix);

std::vector<int> calculate_solution_commonness(const std::vector<int> solution_matrix,
                                               const unsigned n_sites,
                                               const unsigned n_species);

// [[Rcpp::export]]
IntegerMatrix calculate_solution_commonness_site_rcpp(const IntegerMatrix solution_matrix,
                                                      const IntegerMatrix solution_commonness,
                                                      const int site);
std::vector<int> calculate_solution_commonness_site(const std::vector<int> solution_matrix,
                                                    const std::vector<int> solution_commonness,
                                                    const unsigned n_sites,
                                                    const unsigned n_species,
                                                    const unsigned site);

// [[Rcpp::export]]
IntegerMatrix calculate_solution_commonness_species_site_rcpp(const IntegerMatrix solution_matrix,
                                                      const IntegerMatrix solution_commonness,
                                                      const int site, const int species);

// [[Rcpp::export]]
void update_solution_commonness_site_rcpp(const IntegerMatrix solution_matrix,
                                          IntegerMatrix &solution_commonness,
                                          const unsigned site);
void update_solution_commonness_site(const std::vector<int> solution_matrix,
                                     std::vector<int> &solution_commonness,
                                     const unsigned n_sites,
                                     const unsigned n_species,
                                     const unsigned site);
#endif // CALCULATE_SOLUTION_COMMONNESS_H
