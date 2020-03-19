#ifndef MH_OPTIMIZER_H
#define MH_OPTIMIZER_H
#include "Rcpp.h"

using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
List optimizer(IntegerVector alpha_list,
               const unsigned total_gamma,
               NumericMatrix target,
               const unsigned max_iterations = 20000,
               unsigned long seed = 0, bool verbose = true,
               const double increment = 0.01);

// [[Rcpp::export]]
double calc_energy(const NumericMatrix solution_commonness,
                   const NumericMatrix solution_commonness_target);

// Helper functions
NumericMatrix gen_init_solution(const IntegerVector &alpha_list,
                                const unsigned gamma_diversity);

void update_metrics(double &accepted,
                    const double energy, const unsigned iter,
                    std::vector<double> &energy_vector,
                    std::vector<unsigned> &iteration_count,
                    std::vector<double> &acceptance_rate);

#endif // MH_OPTIMIZER_H
