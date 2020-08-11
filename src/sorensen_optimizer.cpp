#include "sorensen_optimizer.h"
#include <random>
#include <chrono>
#include <limits>
#include <cmath>
#include "calculate_solution_commonness.h"
#include "backtracking.h"
#include "sorensen_minconf.h"
#include "rcpp_sample.h"
//omp_set_nested(1);


// MSP start 
List optimizer_min_conf0_sorensen(IntegerVector alpha_list, const unsigned total_gamma,
                                  NumericMatrix target, IntegerMatrix fixed_species,
                                  IntegerMatrix partial_solution,
                                  const unsigned max_iterations,
                                  const unsigned patience,
                                  const double energy_threshold,
                                  unsigned long seed, bool verbose, std::string norm,
                                  const unsigned p)
{
    Sorensen_MinConf mc(as<std::vector<unsigned> >(alpha_list),
                        total_gamma,
                        as<std::vector<double> >(target),
                        as<std::vector<int> >(fixed_species),
                        as<std::vector<int> >(partial_solution),
                        norm);
    mc.p = p;
    long iter = max_iterations - mc.optimize0_sorensen(max_iterations, energy_threshold, seed, patience);
    const unsigned n_sites = alpha_list.size();
    IntegerMatrix solution(total_gamma, n_sites);
    
    for (unsigned site = 0; site < n_sites; site++) {
        for (unsigned species = 0; species < total_gamma; species++) {
            solution(species, site) = mc.solution[site][species];
        }
    }
    
    const auto energy_normalizer = mc.calc_energy_random_solution_sorensen(1000);
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
        Rcout << "\n Warning: this solution does not neccessarily correspond to the lowest energy. \n";
    }
    
    return(results);
}


std::vector<unsigned> calc_min_conflict_species_sorensen(const unsigned site,
                                                         const IntegerMatrix current_solution,
                                                         const NumericMatrix target)
{
    const double epsilon = 0.00001;
    const auto gamma_div = current_solution.nrow();
    const NumericMatrix sorensen = calculate_solution_sorensen_rcpp(current_solution);
    auto solution = current_solution;
    double energy = std::numeric_limits<double>::max(); // makes sure that the first energy_ is smaller
    std::vector<unsigned> min_conflict_species;
    
    // try for each species at this site where the enery would be minimal
    for (unsigned species = 0; species < gamma_div; species++) {
        if (solution(species, site) == 1) {
            continue;
        }
        solution(species, site) = 1;
        
        const NumericMatrix sorensen_new =
            calculate_solution_sorensen_site_rcpp(solution, sorensen, site + 1); // _rcpp functions start indexing at 1 //
            // MSP before... calculate_solution_sorensen_rcpp(solution); 
        // MSP my approach here differs from Sebastian work. Please double-check! 
        
        double energy_ = calc_energy_sorensen(sorensen_new, target);
        
        if (energy_ < energy) {
            min_conflict_species.clear(); // found better fitting species delete other
            min_conflict_species.push_back(species);
            energy = energy_;
        } else if (fabs(energy_ - energy) <= epsilon) {
            min_conflict_species.push_back(species); // as good as others, add this species
        }
        solution(species, site) = 0;
    }
    
    return min_conflict_species;
}

/*
double calc_random_energy_sorensen(unsigned n, IntegerVector alpha_list, const unsigned total_gamma, NumericMatrix target, unsigned long seed, std::string norm)
{
    Sorensen_MinConf mc(as<std::vector<unsigned> >(alpha_list),
               total_gamma,
               as<std::vector<double> >(target),
               std::vector<int> (),
               std::vector<int> (),
               norm);
    mc.setSeed(seed);
    return mc.calc_energy_random_solution_sorensen(n);
}
 */ 

double calc_energy_sorensen(const NumericMatrix solution_commonness, 
                   const NumericMatrix solution_commonness_target)
{
    // calculate the difference between target and current solution
    return(sum(abs(na_omit(solution_commonness - solution_commonness_target))) /
           (double)sum(na_omit(solution_commonness_target)));
}

/* // This function is not used by min_conf_0, ... 
IntegerMatrix gen_init_solution(const IntegerVector alpha_list,
                                const unsigned gamma_diversity)
{
    const auto n_sites = alpha_list.size(); // number of cols
    
    IntegerMatrix current_solution(gamma_diversity, n_sites);
    
    for (int site = 0; site < n_sites; site++) {
        const unsigned alpha = static_cast<unsigned>(alpha_list[site]);
        IntegerVector alpha_div(gamma_diversity);
        for (unsigned i = 0; i < alpha; i++) {
            alpha_div[i] = 1;
        }
        std::random_shuffle(alpha_div.begin(), alpha_div.end()); ///BUG: seed is ignored, here
        current_solution(_, site) = alpha_div;
    }
    
    return current_solution;
}
 */ 
    