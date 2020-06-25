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
                        const double energy_threshold = 0.0,
                        unsigned long seed = 0, bool verbose = true);

// [[Rcpp::export]]
List optimizer_min_conf0(IntegerVector alpha_list, const unsigned total_gamma,
                         IntegerMatrix target, IntegerMatrix fixed_species,
                         IntegerMatrix partial_solution,
                         const unsigned max_iterations,
                         const double energy_threshold = 0.0,
                         unsigned long seed = 0, bool verbose = true, std::string norm = "sum",
                         const unsigned p = 1);

// [[Rcpp::export]]
List optimizer_min_conf1(IntegerVector alpha_list, const unsigned total_gamma,
                         IntegerMatrix target, IntegerMatrix fixed_species,
                         IntegerMatrix partial_solution,
                         const unsigned max_iterations,
                         const double energy_threshold = 0.0,
                         unsigned long seed = 0, bool verbose = true, std::string norm = "sum");


// [[Rcpp::export]]
List optimizer_min_conf2(IntegerVector alpha_list, const unsigned total_gamma,
                         IntegerMatrix target, IntegerMatrix fixed_species,
                         IntegerMatrix partial_solution,
                         const unsigned max_iterations,
                         const double energy_threshold = 0.0,
                         unsigned long seed = 0, bool verbose = true, std::string norm = "sum");

// [[Rcpp::export]]
List optimizer_backtracking(IntegerVector alpha_list,
                            const unsigned total_gamma,
                            IntegerMatrix target,
                            const unsigned max_iterations = 2000,
                            bool verbose = true, std::string norm = "sum");

// [[Rcpp::export]]
double calc_energy(const IntegerMatrix solution_commonness,
                   const IntegerMatrix solution_commonness_target);

// [[Rcpp::export]]
double calc_random_energy(unsigned n,
                          IntegerVector alpha_list,
                          const unsigned total_gamma,
                          IntegerMatrix target,
                          unsigned long seed = 0, std::string norm = "sum");

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
