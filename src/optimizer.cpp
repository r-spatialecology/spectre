#include "optimizer.h"
#include <random>
#include <chrono>
#include <limits>
#include <cmath>
#include "calculate_solution_commonness.h"
#include "backtracking.h"
#include "minconf.h"
#include "rcpp_sample.h"
//omp_set_nested(1);

List optimizer_min_conf(IntegerVector alpha_list, const unsigned total_gamma,
                        IntegerMatrix target, const unsigned max_iterations,
                        const double energy_threshold,
                        unsigned long seed, bool verbose)
{
    RNGScope scope; // needed for debugging in Qt Creator
    
    IntegerMatrix current_solution = gen_init_solution(alpha_list, total_gamma);
    
    auto solution_ = as<std::vector<int> >(current_solution);
    // get dimensions of matrix
    const unsigned n_col = alpha_list.length();
    
    // check for sanity to avoid an infinite loop later on
    // assert() is actually bad practice in Rcpp because it also kills R...
    const unsigned check = sum(alpha_list);
    assert(check > 0);
    assert(check < total_gamma * n_col);
    
    std::vector<double> energy_vector;
    energy_vector.reserve(max_iterations + 1);
    std::vector<unsigned> iteration_count;
    iteration_count.reserve(max_iterations + 1);
    std::vector<double> acceptance_rate;
    acceptance_rate.reserve(max_iterations + 1);
    
    // Random number generator
    if (seed == 0) {
        seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    }
    std::mt19937 rng(seed);
    std::uniform_int_distribution<unsigned> site_dist(0, current_solution.ncol() - 1);
    
    for (unsigned iter = 0; iter < max_iterations; iter++) {
        const auto site = site_dist(rng);
        
        solution_ = as<std::vector<int> >(current_solution);
        
        // Get index of all species of this site
        std::vector<unsigned> species_idx;
        for (unsigned species = 0; species < total_gamma; species++) {
            if (current_solution(species, site)) {
                species_idx.push_back(species);
            }
        }
        
        if (species_idx.size()) { // if the site is not empty
            // Choose a random species and remove it
            std::uniform_int_distribution<unsigned> rnd_species(0, species_idx.size() - 1);
            auto random_species = rnd_species(rng);
            current_solution(species_idx[random_species], site) = 0;
            
            // calculate the best-fitting species for this site
            auto min_conflict_species = calc_min_conflict_species(site,
                                                                  current_solution,
                                                                  target);
            // choose a random best-fitting species and add it
            std::uniform_int_distribution<unsigned> rnd_new_species(0, min_conflict_species.size() - 1);
            const auto s = rnd_new_species(rng);
            current_solution(min_conflict_species[s], site) = 1;
        }
        
        solution_ = as<std::vector<int> >(current_solution);
        
        const auto current_commonness = calculate_solution_commonness_rcpp(current_solution);
        const auto current_energy = calc_energy(current_commonness, target);
        update_metrics(current_energy, iter, energy_vector, iteration_count);
        
        if (current_energy <= energy_threshold) {
            break;
        }
    }
    
    // generate dataframe for i and energy
    DataFrame measures_df = DataFrame::create(_["i"] = iteration_count,
                                              _["energy"] = energy_vector);
    
    List results = List::create(Rcpp::Named("optimized_grid") = current_solution,
                                Rcpp::Named("energy") = measures_df);
    if (verbose) {
        double best_energy = *std::min_element(energy_vector.begin(), energy_vector.end());
        double worst_energy = *std::max_element(energy_vector.begin(), energy_vector.end());
        Rcout << "\n > Optimization finished with lowest energy = " << best_energy << " %"
              << " (highest energy was: " << worst_energy << " %, improved by: "
              << worst_energy - best_energy << " %)";
    }
    return(results);
}



List optimizer_min_conf0(IntegerVector alpha_list, const unsigned total_gamma,
                         IntegerMatrix target, IntegerMatrix fixed_species,
                         IntegerMatrix partial_solution,
                         const unsigned max_iterations,
                         const unsigned patience,
                         const double energy_threshold,
                         unsigned long seed, bool verbose, std::string norm,
                         const unsigned p)
{
    MinConf mc(as<std::vector<unsigned> >(alpha_list),
               total_gamma,
               as<std::vector<int> >(target),
               as<std::vector<int> >(fixed_species),
               as<std::vector<int> >(partial_solution),
               norm);
    mc.p = p;
    long iter = max_iterations - mc.optimize0(max_iterations, energy_threshold, seed, patience);
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

List optimizer_min_conf1(IntegerVector alpha_list, const unsigned total_gamma,
                         IntegerMatrix target, IntegerMatrix fixed_species,
                         IntegerMatrix partial_solution,
                         const unsigned max_iterations,
                         const double energy_threshold,
                         unsigned long seed, bool verbose, std::string norm)
{
    MinConf mc(as<std::vector<unsigned> >(alpha_list),
               total_gamma,
               as<std::vector<int> >(target),
               as<std::vector<int> >(fixed_species),
               as<std::vector<int> >(partial_solution),
               norm);
    long iter = max_iterations - mc.optimize1(max_iterations, energy_threshold, seed);
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
        Rcout << "\n > Optimization finished with lowest energy = " << best_energy
              << " (highest energy was: " << worst_energy << " , improved by: "
              << worst_energy - best_energy << " )";
    }
    
    return(results);
}

List optimizer_min_conf2(IntegerVector alpha_list, const unsigned total_gamma,
                         IntegerMatrix target, IntegerMatrix fixed_species,
                         IntegerMatrix partial_solution,
                         const unsigned max_iterations,
                         const double energy_threshold,
                         unsigned long seed, bool verbose, std::string norm)
{
    MinConf mc(as<std::vector<unsigned> >(alpha_list),
               total_gamma,
               as<std::vector<int> >(target),
               as<std::vector<int> >(fixed_species),
               as<std::vector<int> >(partial_solution),
               norm);
    long iter = max_iterations - mc.optimize2(max_iterations, energy_threshold, seed);
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

List optimizer_backtracking(IntegerVector alpha_list, const unsigned total_gamma,
                            IntegerMatrix target, const unsigned max_iterations, bool verbose, std::string norm)
{
    RNGScope scope; // needed for debugging in Qt Creator
    
    Backtracking bt(as<std::vector<unsigned> >(alpha_list),
                    total_gamma,
                    as<std::vector<int> >(target),
                    std::vector<int>(),
                    std::vector<int>(),
                    norm);
    long iter = max_iterations - bt.optimize(max_iterations);
    const unsigned n_sites = alpha_list.size();
    IntegerMatrix solution(total_gamma, n_sites);
    
    for (unsigned site = 0; site < n_sites; site++) {
        for (unsigned species = 0; species < total_gamma; species++) {
            solution(species, site) = bt.solution[site][species];
        }
    }
    
    IntegerVector solved_sites;
    solved_sites.assign(bt.solved_sites.begin(), bt.solved_sites.end());
    
    const auto solution_commonness = calculate_solution_commonness_rcpp(solution);
    const double energy = calc_energy(solution_commonness, target);
    List results = List::create(Rcpp::Named("optimized_grid") = solution,
                                Rcpp::Named("energy") = energy,
                                Rcpp::Named("solved_sites") = solved_sites);
    if (verbose) {
        Rcout << "\n > Optimization stopped after " << iter
              << " steps";
    }
    return(results);
}

double calc_energy(const IntegerMatrix solution_commonness,
                   const IntegerMatrix solution_commonness_target)
{
    // calculate the difference between target and current solution
    return(sum(abs(na_omit(solution_commonness - solution_commonness_target))) /
           (double)sum(na_omit(solution_commonness_target)));
}

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

void update_metrics(const double energy, const unsigned iter,
                    std::vector<double> &energy_vector,
                    std::vector<unsigned> &iteration_count)
{
    energy_vector.push_back(energy);
    iteration_count.push_back(iter + 1);
}



std::vector<unsigned> calc_min_conflict_species(const unsigned site,
                                                const IntegerMatrix current_solution,
                                                const IntegerMatrix target)
{
    const double epsilon = 0.00001;
    const auto gamma_div = current_solution.nrow();
    const IntegerMatrix commoness = calculate_solution_commonness_rcpp(current_solution);
    auto solution = current_solution;
    double energy = std::numeric_limits<double>::max(); // makes sure that the first energy_ is smaller
    std::vector<unsigned> min_conflict_species;
    
    // try for each species at this site where the enery would be minimal
    for (unsigned species = 0; species < gamma_div; species++) {
        if (solution(species, site) == 1) {
            continue;
        }
        solution(species, site) = 1;
        const IntegerMatrix commoness_new =
            calculate_solution_commonness_site_rcpp(solution, commoness, site + 1); // _rcpp functions start indexing at 1
        
        double energy_ = calc_energy(commoness_new, target);
        
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

double calc_random_energy(unsigned n, IntegerVector alpha_list, const unsigned total_gamma, IntegerMatrix target, unsigned long seed, std::string norm)
{
    MinConf mc(as<std::vector<unsigned> >(alpha_list),
               total_gamma,
               as<std::vector<int> >(target),
               std::vector<int> (),
               std::vector<int> (),
               norm);
    mc.setSeed(seed);
    return mc.calc_energy_random_solution(n);
}
