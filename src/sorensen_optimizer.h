#ifndef SORENSEN_OPTIMIZER_H
#define SORENSEN_OPTIMIZER_H
#include "Rcpp.h"

using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]


// [[Rcpp::export]]
List optimizer_min_conf0_sorensen(IntegerVector alpha_list, const unsigned total_gamma,
                                  NumericMatrix target, IntegerMatrix fixed_species,
                                  IntegerMatrix partial_solution,
                                  const unsigned max_iterations,
                                  const unsigned patience = 2000,
                                  const double energy_threshold = 0.0,
                                  unsigned long seed = 0, bool verbose = true, std::string norm = "sum",
                                  const unsigned p = 1);

std::vector<unsigned> calc_min_conflict_species_sorensen(const unsigned site,
                                                         const IntegerMatrix current_solution,
                                                         const NumericMatrix target);


// [[Rcpp::export]]
double calc_energy_sorensen(const NumericMatrix solution_commonness,
                            const NumericMatrix solution_commonness_target);

/* commmented out on 2020-08-08, WIP, may not be a good idea... 




// [[Rcpp::export]]
double calc_random_energy_sorensen(unsigned n,
                          IntegerVector alpha_list,
                          const unsigned total_gamma,
                          IntegerMatrix target,
                          unsigned long seed = 0, std::string norm = "sum");


*/

/* 

// Helper functions
// [[Rcpp::export]]
IntegerMatrix gen_init_solution(const IntegerVector alpha_list,
                                const unsigned gamma_diversity);


void update_metrics(const double energy, const unsigned iter,
                    std::vector<double> &energy_vector,
                    std::vector<unsigned> &iteration_count);
 */

#endif // SORENSEN_OPTIMIZER_H
