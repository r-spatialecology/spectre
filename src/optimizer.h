#ifndef MH_OPTIMIZER_H
#define MH_OPTIMIZER_H
#include "Rcpp.h"

using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
List optimizer_min_conf(IntegerVector alpha_list, const unsigned total_gamma,
                         IntegerMatrix target, IntegerMatrix fixed_species,
                         IntegerMatrix partial_solution,
                         const unsigned max_iterations,
                         const double energy_threshold = 0.0,
                         unsigned long seed = 0, bool verbose = true);
// [[Rcpp::export]]
IntegerMatrix calculate_solution_commonness_rcpp(IntegerMatrix solution_matrix);
#endif // MH_OPTIMIZER_H
