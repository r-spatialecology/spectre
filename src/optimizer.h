#ifndef MH_OPTIMIZER_H
#define MH_OPTIMIZER_H
#include "Rcpp.h"

using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
List optimizer_min_conf(IntegerVector alpha_list,
                        const unsigned total_gamma,
                        IntegerMatrix target,
                        const unsigned max_iterations = 2000,
                        const double energy_theshold = 0.1,
                        unsigned long seed = 0, bool verbose = true);

// [[Rcpp::export]]
List optimizer_backtracking(IntegerVector alpha_list,
               const unsigned total_gamma,
               IntegerMatrix target,
               const unsigned max_iterations = 2000,
               unsigned long seed = 0, bool verbose = true);

// [[Rcpp::export]]
double calc_energy(const IntegerMatrix solution_commonness,
                   const IntegerMatrix solution_commonness_target);

// Helper functions
// [[Rcpp::export]]
IntegerMatrix gen_init_solution(const IntegerVector alpha_list,
                                const unsigned gamma_diversity);

std::vector<unsigned> calc_min_conflict_species(const unsigned site,
                                                const IntegerMatrix current_solution,
                                                const IntegerMatrix target);

void update_metrics(const double energy, const unsigned iter,
                    std::vector<double> &energy_vector,
                    std::vector<unsigned> &iteration_count);

#endif // MH_OPTIMIZER_H
