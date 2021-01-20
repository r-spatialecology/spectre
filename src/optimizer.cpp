#include "optimizer.h"
#include <random>
#include <chrono>
#include <limits>
#include <cmath>
#include "minconf.h"

List optimizer_min_conf(IntegerVector alpha_list, const unsigned total_gamma,
                        IntegerMatrix target, IntegerMatrix fixed_species,
                        IntegerMatrix partial_solution,
                        const unsigned max_iterations,
                        unsigned long seed, bool verbose, bool interruptible)
{
    MinConf mc(as<std::vector<unsigned> >(alpha_list),
               total_gamma,
               as<std::vector<int> >(target),
               as<std::vector<int> >(fixed_species),
               as<std::vector<int> >(partial_solution));

    long iter = max_iterations - mc.optimize(max_iterations,
                                             seed, verbose, interruptible);


    if (iter == max_iterations - mc.RET_ABORT) {
        Rcout << "The processing was aborted by the user. \n";
        return List();
    }

    const unsigned n_sites = alpha_list.size();
    IntegerMatrix solution(total_gamma, n_sites);

    for (unsigned site = 0; site < n_sites; site++) {
        for (unsigned species = 0; species < total_gamma; species++) {
            solution(species, site) = mc.solution[site][species];
        }
    }

    // generate dataframe for i and error
    DataFrame measures_df = DataFrame::create(_["i"] = mc.iteration_count,
            _["error"] = mc.error_vector);

    List results = List::create(Rcpp::Named("optimized_grid") = solution,
                                Rcpp::Named("error") = measures_df);
    if (verbose) {
        double best_error = *std::min_element(mc.error_vector.begin(), mc.error_vector.end());
        double worst_error = *std::max_element(mc.error_vector.begin(), mc.error_vector.end());
        Rcout << "\n > Optimization finished with lowest absolute error = " << best_error << " %"
              << " (highest absolute error was: " << worst_error << " %, improved by: "
              << worst_error - best_error << " %)";
    }

    if (!mc.solution_has_best_error) {
        Rcout << "\n Warning: this solution does not neccessarily correnpond to the lowest absolute error. \n";
    }

    return(results);
}

IntegerMatrix calculate_solution_commonness_rcpp(const IntegerMatrix solution_matrix) {
    // create results matrix, filled with NAs (== -1)
    // and translate solution_matrix to 2-D std::vector
    const unsigned nsites = solution_matrix.ncol();
    const unsigned gamma_div = solution_matrix.nrow();
    std::vector<std::vector<int> > solution_mat(nsites, std::vector<int>(gamma_div, -1));

    for (unsigned site = 0; site < nsites; site++) {
        for (unsigned species = 0; species < gamma_div; species++) {
            solution_mat[site][species] = solution_matrix(species, site);
        }
    }

    IntegerMatrix result(nsites, nsites);
    std::fill(result.begin(), result.end(), NumericVector::get_na() );

    const auto commonness_mat = MinConf::calculate_commonness(solution_mat, nsites);
    for (unsigned site = 0; site < nsites; site++) {
        for (unsigned other_site = 0; other_site < nsites; other_site++) {
            if(commonness_mat[site][other_site] != -1) { // NA
                result(site, other_site) = commonness_mat[site][other_site];
            }
        }
    }
    return result;
}

unsigned calc_error_random_solution(const unsigned n, IntegerVector alpha_list, const unsigned total_gamma,
                                     IntegerMatrix target, unsigned long seed)
{
    MinConf mc(as<std::vector<unsigned> >(alpha_list), total_gamma, as<std::vector<int> >(target));
    mc.setSeed(seed);
    return mc.calc_error_random_solution(n);
}
