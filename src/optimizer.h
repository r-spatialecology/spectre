#ifndef MH_OPTIMIZER_H
#define MH_OPTIMIZER_H
#include "Rcpp.h"

using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
List mh_optimizer(IntegerVector alpha_list,
                  const unsigned total_gamma,
                  IntegerMatrix target,
                  const unsigned max_iterations = 1000,
                  unsigned long seed = 0, bool verbose = true,
                  double base_probability_jump = 1.0);

// [[Rcpp::export]]
List mh_optimizer_neutral(IntegerVector alpha_list,
                          const unsigned total_gamma,
                          IntegerMatrix solution_commonness_target,
                          const unsigned max_iterations = 1000,
                          unsigned long seed = 0, bool verbose = true);

// [[Rcpp::export]]
double calc_energy(const IntegerVector solution_commonness,
                   const IntegerVector solution_commonness_target);

// [[Rcpp::export]]
void species_swap_rcpp(IntegerMatrix &mat, const IntegerVector species, unsigned site);

// Helper functions
IntegerMatrix gen_init_solution(const unsigned n_row, const unsigned n_col,
                                const IntegerVector &alpha_list);

IntegerMatrix gen_init_solution(const IntegerVector &alpha_diversities,
                                const unsigned gamma_diversity);

unsigned long factorial_loop(unsigned n);
unsigned long factorial_ln(unsigned n);
unsigned long factorial_rec(unsigned n);
unsigned long factorial(unsigned n);

void update_metrics(double &accepted,
                    const double energy, const unsigned iter,
                    std::vector<double> &energy_vector,
                    std::vector<unsigned> &iteration_count,
                    std::vector<double> &acceptance_rate);

double jump_probability(const double new_energy, const double energy, double base_prob);

#endif // MH_OPTIMIZER_H
