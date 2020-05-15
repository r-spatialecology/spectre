#include "optimizer.h"
#include <random>
#include <chrono>
#include <limits>
#include <math.h>
#include "calculate_solution_commonness.h"
#include "constraint_satisfaction_problem.h"
#include "rcpp_sample.h"
//omp_set_nested(1);

List optimizer_min_conf(IntegerVector alpha_list, const unsigned total_gamma,
                        IntegerMatrix target, const unsigned max_iterations,
                        const double energy_theshold,
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

        if (current_energy <= energy_theshold) {
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

List optimizer_backtracking(IntegerVector alpha_list, const unsigned total_gamma,
               IntegerMatrix target, const unsigned max_iterations,
               unsigned long seed, bool verbose)
{
    RNGScope scope; // needed for debugging in Qt Creator

    Constraint_satisfaction_problem csr(as<std::vector<unsigned> >(alpha_list),
                                        total_gamma,
                                        as<std::vector<int> >(target));
    long iter = max_iterations - csr.optimize(max_iterations);
    const unsigned n_sites = alpha_list.size();
    IntegerMatrix solution(total_gamma, n_sites);

    for (unsigned site = 0; site < n_sites; site++) {
        for (unsigned species = 0; species < total_gamma; species++) {
            solution(species, site) = csr.solution[site][species];
        }
    }

    const auto solution_commonness = calculate_solution_commonness_rcpp(solution);
    const double energy = calc_energy(solution_commonness, target);
    List results = List::create(Rcpp::Named("optimized_grid") = solution,
                                Rcpp::Named("energy") = energy);
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

        //        auto commoness_ = as<std::vector<int> >(commoness_new);
        //        auto target_ = as<std::vector<int> >(target);
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



std::vector<unsigned> calc_constraints(const IntegerMatrix target)
{
    const unsigned n_sites = target.nrow();
    std::vector<unsigned> constraints(n_sites);

    for (unsigned site = 0; site < n_sites; site++) {
        for (int j = 0; j < site - 1; j++) {
            constraints[site] += target(site, j);
        }

        for (int j = site + 1; j < n_sites; j++) {
            constraints[site] += target(j, site);
        }
    }
}
