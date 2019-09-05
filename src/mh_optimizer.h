#ifndef MH_OPTIMIZER_H
#define MH_OPTIMIZER_H
#include "Rcpp.h"

using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
List mh_optimizer(IntegerVector alpha_list,
                  const int total_gamma,
                  IntegerMatrix target,
                  const double acceptance_rate_threshold = 0.2,
                  const unsigned max_iterations = 1000,
                  const unsigned burn_in = 200,
                  unsigned long seed = 0);

// [[Rcpp::export]]
double calc_energy(const IntegerVector target,
                          const IntegerVector solution_commonness);

// [[Rcpp::export]]
void species_swap_rcpp(IntegerMatrix &mat, const IntegerVector species, unsigned site);
#endif // MH_OPTIMIZER_H
