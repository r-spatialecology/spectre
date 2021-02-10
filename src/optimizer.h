#ifndef MH_OPTIMIZER_H
#define MH_OPTIMIZER_H
#include "Rcpp.h"

using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
List optimizer_min_conf(const IntegerVector alpha_list,
                        const unsigned total_gamma, const IntegerMatrix target,
                        const unsigned max_iterations,
                        const IntegerMatrix partial_solution,
                        const IntegerMatrix fixed_species,
                        const unsigned long seed = 0, const bool verbose = true,
                        const bool interruptible = true);
// [[Rcpp::export]]
IntegerMatrix
calculate_solution_commonness_rcpp(const IntegerMatrix solution_matrix);
#endif // MH_OPTIMIZER_H
