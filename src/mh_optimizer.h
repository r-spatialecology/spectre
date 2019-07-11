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
                  const unsigned max_iterations = 100000,
                  const unsigned burn_in = 20000,
                  unsigned long seed = 0);

#endif // MH_OPTIMIZER_H
