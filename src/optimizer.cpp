#include "optimizer.h"
#include <random>
#include <chrono>
#include <limits>
#include <cmath>
#include "minconf.h"
//omp_set_nested(1);

List optimizer_min_conf(IntegerVector alpha_list, const unsigned total_gamma,
                        IntegerMatrix target, IntegerMatrix fixed_species,
                        IntegerMatrix partial_solution,
                        const unsigned max_iterations,
                        const double energy_threshold,
                        unsigned long seed, bool verbose, bool interruptible)
{
    MinConf mc(as<std::vector<unsigned> >(alpha_list),
               total_gamma,
               as<std::vector<int> >(target),
               as<std::vector<int> >(fixed_species),
               as<std::vector<int> >(partial_solution));

    long iter = max_iterations - mc.optimize(max_iterations,
                                             energy_threshold,
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

    const auto energy_normalizer = mc.calc_energy_random_solution(1000);
    for (unsigned i = 0; i < mc.energy_vector.size(); i++) {
        mc.energy_vector[i] /= energy_normalizer;
    }

    // generate dataframe for i and energy
    DataFrame measures_df = DataFrame::create(_["i"] = mc.iteration_count,
            _["energy"] = mc.energy_vector);

    List results = List::create(Rcpp::Named("optimized_grid") = solution,
                                Rcpp::Named("energy") = measures_df);
    if (verbose) {
        double best_energy = *std::min_element(mc.energy_vector.begin(), mc.energy_vector.end());
        double worst_energy = *std::max_element(mc.energy_vector.begin(), mc.energy_vector.end());
        Rcout << "\n > Optimization finished with lowest energy = " << best_energy << " %"
              << " (highest energy was: " << worst_energy << " %, improved by: "
              << worst_energy - best_energy << " %)";
    }

    if (!mc.solution_has_best_energy) {
        Rcout << "\n Warning: this solution does not neccessarily correnpond to the lowest energy. \n";
    }

    return(results);
}

IntegerMatrix calculate_solution_commonness_rcpp(const IntegerMatrix solution_matrix) {
    // create results matrix, filled with NAs (== -1)
    // and translate solution_matrix to 2-D std::vector
    const unsigned nsites = solution_matrix.ncol();
    const unsigned gamma_div = solution_matrix.nrow();
    std::vector<std::vector<int> > solution_mat(nsites, std::vector<int>(gamma_div, -1));
    std::vector<std::vector<int> > commonness_mat(nsites, std::vector<int>(nsites));
    for (unsigned site = 0; site < nsites; site++) {
        for (unsigned species = 0; species < gamma_div; species++) {
            solution_mat[site][species] = solution_matrix(species, site);
        }
    }
    IntegerMatrix result(nsites, nsites);
    std::fill(result.begin(), result.end(), NumericVector::get_na() );

    for (unsigned site = 0; site < nsites; site++) {
        MinConf::update_solution_commonness_site(solution_mat, commonness_mat,
                                                 nsites, gamma_div, site);
        for (unsigned other_site = 0; other_site < nsites; other_site++) {
            if(commonness_mat[site][other_site] != -1) { // NA
                result(site, other_site) = commonness_mat[site][other_site];
            }
        }
    }
    return result;
}
