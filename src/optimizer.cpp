#include "optimizer.h"
#include <random>
#include <chrono>
#include <cmath>
#include "calculate_solution_commonness.h"
#include "rcpp_sample.h"
//omp_set_nested(1);

List mh_optimizer(IntegerVector alpha_list,
                  const unsigned total_gamma,
                  IntegerMatrix solution_commonness_target,
                  NumericVector species_prop,
                  const unsigned max_iterations,
                  const double energy_theshold,
                  double base_probability_jump,
                  unsigned long seed,
                  bool verbose)
{
    RNGScope scope; // needed for debugging
    // get dimensions of matrix
    const unsigned n_species = total_gamma;
    const unsigned n_sites = alpha_list.length();

    // check for sanity to avoid an infinite loop later on
    // assert() is actually bad practice in Rcpp because it also kills R...
    const unsigned check = sum(alpha_list);
    assert(check > 0);
    assert(check < total_gamma * n_sites);

    std::vector<double> energy_vector;
    energy_vector.reserve(max_iterations + 1);
    std::vector<unsigned> iteration_count;
    iteration_count.reserve(max_iterations + 1);
    std::vector<double> acceptance_rate;
    acceptance_rate.reserve(max_iterations + 1);

    // Random number generator
    if (seed < 1) {
        seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    }
    std::mt19937 rng(seed);
    std::uniform_int_distribution<unsigned> site_dist(0, n_sites - 1);
    std::uniform_real_distribution<double> real_dist_01(0, 1);

    // generate initial solution, loop changes the alpha number of species in each site to present (i.e. 1)
    ///TODO: do this in parallel (chains)
    IntegerMatrix current_solution = gen_init_solution(n_species, n_sites, alpha_list,
                                                       species_prop, rng, real_dist_01);
    IntegerMatrix best_solution = current_solution;

    // calculate the site x site commonness for the current solution
    auto solution_commonness = calculate_solution_commonness_rcpp(current_solution);

    // calculate the difference between target and current solution
    double energy = calc_energy(solution_commonness, solution_commonness_target);

    energy_vector.push_back(energy);
    acceptance_rate.push_back(0.0);
    iteration_count.push_back(0);
    double accepted = 0;

    for (unsigned iter = 0; iter < max_iterations; iter++) {
        // swap species occurrence at a random site
        bool find_site = false;
        IntegerVector species_to_swap;
        unsigned random_site = 0;
        while (!find_site) {
            random_site = site_dist(rng);
            species_to_swap = get_swap_rows_rcpp(
                        current_solution(_, random_site), true);

            double swap_prob = 1.0; // - (1-base_probability_jump);
            double swap_prob_0 = species_prop[species_to_swap[0]];
            double swap_prob_1 = species_prop[species_to_swap[1]];
            if (current_solution(species_to_swap[0], random_site) == 1) {
                if (swap_prob_0 > swap_prob_1) {
                    swap_prob_0 = 1 - swap_prob_0; // i.e. prob. to become a zero
                    swap_prob = 0.5 * (swap_prob_0 + swap_prob_1);
                }
            } else if (swap_prob_0 < swap_prob_1){
                swap_prob_1 = 1 - swap_prob_1;
                swap_prob = 0.5 * (swap_prob_0 + swap_prob_1);
            }

            const double random = real_dist_01(rng);
            if (random < swap_prob) {
                find_site = true;
            }
        }


        std::swap(current_solution(species_to_swap[0], random_site),
                current_solution(species_to_swap[1], random_site));

        // calculate commoness and energy for new solution
        const auto new_solution_commonness = calculate_solution_commonness_site_rcpp(
                    current_solution, solution_commonness, random_site + 1);

        const double new_energy = calc_energy(new_solution_commonness, solution_commonness_target);

        // if proposed solution is better or same, accept, else calculate jump probability
        if (new_energy <= energy) {
            best_solution = current_solution;
            energy = new_energy;
            accepted++;
        } else if (real_dist_01(rng) < jump_probability(new_energy, energy, base_probability_jump)) {
            // accept anyway?
            // best_solution = current_solution; // not the best solution, however
            energy = new_energy;
            accepted++;
        } else {
            // not accepted, undo changes
            std::swap(current_solution(species_to_swap[0], random_site),
                    current_solution(species_to_swap[1], random_site));
        }
        //        update_metrics(accepted, energy, new_energy, iter, energy_vector,
        //                       iteration_count, acceptance_rate);
        // print progress
        if (verbose) {
            Rcout << " Progress: max_runs: " << iter + 1 << "/" << max_iterations <<
                     " || energy =  " << energy << " %\r";
            Rcout.flush();
        }

        update_metrics(accepted, energy, iter, energy_vector,
                       iteration_count, acceptance_rate);
        if (energy <= energy_theshold) {
            break;
        }
    }

    // generate dataframe for i and energy
    DataFrame measures_df = DataFrame::create(_["i"] = iteration_count,
            _["energy"] = energy_vector);

    List results = List::create(Rcpp::Named("optimized_grid") = best_solution,
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

double calc_energy(const IntegerVector solution_commonness, const IntegerVector solution_commonness_target)
{
    // calculate the difference between target and current solution
    unsigned n_diff = 0;
    unsigned n_all = 0;
    for (int i = 0; i < solution_commonness.size(); i++) {
        if (solution_commonness_target[i] == NA_INTEGER)
            continue;
        n_diff += std::abs(solution_commonness[i] - solution_commonness_target[i]);
        n_all += solution_commonness_target[i];
    }
    return(static_cast<double>(n_diff) / static_cast<double>(n_all));
    //            return(sum(abs(na_omit(solution_commonness_target - solution_commonness))) /
    //                   (double)sum(na_omit(solution_commonness_target)));
}

void species_swap_rcpp(IntegerMatrix &mat, const IntegerVector species, unsigned site) {
    std::swap(mat(species[0] - 1, site - 1), mat(species[1] - 1, site - 1));
}

IntegerMatrix gen_init_solution(const unsigned n_row, const unsigned n_col, const IntegerVector &alpha_list)
{
    IntegerMatrix current_solution = IntegerMatrix(n_row, n_col);

    const IntegerVector species_seq = seq_len(n_row);
    for (unsigned site = 0; site < n_col; site++) {
        const IntegerVector change_species = sample(species_seq, alpha_list[site], false);
        for (auto species : change_species) {
            current_solution(species, site) = 1;
        }
    }
    return current_solution;
}

IntegerMatrix gen_init_solution(const unsigned n_species, const unsigned n_sites,
                                const IntegerVector &alpha_list, NumericVector w,
                                std::mt19937 &mt, std::uniform_real_distribution<double> &real_dist_01)
{
    IntegerMatrix current_solution = IntegerMatrix(n_species, n_sites);
    std::uniform_int_distribution<unsigned> species_dist(0, n_species - 1);
    auto alpha = alpha_list;
    for (unsigned site = 0; site < n_sites; site++) {
        while (alpha[site])
        {
            const unsigned species = species_dist(mt);
            const double random = real_dist_01(mt);
            const double prob = w[species];
            if (random < w[species]) {
                current_solution(species, site) = 1;
                alpha[site]--;
            }
        }
    }
    return current_solution;
}

IntegerMatrix gen_init_solution(const IntegerVector &alpha_diversities, const unsigned gamma_diversity)
{
    IntegerMatrix current_solution(gamma_diversity,
                                   alpha_diversities.length());

    const IntegerVector species_seq = seq_len(gamma_diversity);
    for (unsigned site = 0; site < alpha_diversities.length(); site++) {
        const IntegerVector change_species = sample(species_seq,
                                                    alpha_diversities[site],
                                                    false);
        for (auto species : change_species) {
            current_solution(species, site) = 1;
        }
    }
    return current_solution;
}

void update_metrics(double &accepted, const double energy, const unsigned iter,
                    std::vector<double> &energy_vector,
                    std::vector<unsigned> &iteration_count,
                    std::vector<double> &acceptance_rate)
{
    energy_vector.push_back(energy);
    iteration_count.push_back(iter + 1);
    acceptance_rate.push_back(accepted / (iter + 1));
}



double jump_probability(const double new_energy, const double energy, double base_prob)
{
    double result = 1.0;
    if (new_energy > energy) {
        result = base_prob * (energy / new_energy);
    }

    return result; // - (1-base_prob));
}

List mh_optimizer_neutral(IntegerVector alpha_list,
                          const unsigned total_gamma,
                          IntegerMatrix solution_commonness_target,
                          NumericVector species_prop,
                          const unsigned max_iterations,
                          unsigned long seed,
                          bool verbose)
{
    RNGScope scope; // needed for debugging
    // get dimensions of matrix
    const unsigned n_row = total_gamma;
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
    if (seed < 1) {
        seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    }
    std::mt19937 rng(seed);
    std::uniform_int_distribution<unsigned> site_dist(0, n_col - 1);
    std::uniform_real_distribution<double> real_dist_01(0, 1);

    // generate initial solution, loop changes the alpha number of species in each site to present (i.e. 1)
    ///TODO: do this in parallel (chains)
    IntegerMatrix current_solution = gen_init_solution(n_row, n_col, alpha_list,
                                                       species_prop, rng, real_dist_01);
    IntegerMatrix best_solution = current_solution;

    // calculate the site x site commonness for the current solution
    auto solution_commonness = calculate_solution_commonness_rcpp(current_solution);

    // calculate the difference between target and current solution
    double energy = calc_energy(solution_commonness, solution_commonness_target);

    energy_vector.push_back(energy);
    acceptance_rate.push_back(0.0);
    iteration_count.push_back(0);
    double accepted = 0;

    for (unsigned iter = 0; iter < max_iterations; iter++) {

        // new solution
        current_solution = gen_init_solution(n_row, n_col, alpha_list,
                                             species_prop, rng, real_dist_01);

        // calculate commoness and energy for new solution
        const auto new_solution_commonness = calculate_solution_commonness_rcpp(current_solution);

        const double new_energy = calc_energy(new_solution_commonness, solution_commonness_target);

        // if proposed solution is better or same, accept, else calculate jump probability
        if (new_energy < energy) {
            best_solution = current_solution;
            energy = new_energy;
        }
        update_metrics(accepted, energy, iter, energy_vector,
                       iteration_count, acceptance_rate);
        // print progress
        if (verbose) {
            Rcout << " Progress: max_runs: " << iter + 1 << "/" << max_iterations <<
                     " || energy = " << energy << " %\r";
            Rcout.flush();
        }
    }

    // generate dataframe for i and energy
    DataFrame measures_df = DataFrame::create(_["i"] = iteration_count,
            _["energy"] = energy_vector);

    List results = List::create(Rcpp::Named("optimized_grid") = best_solution,
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

List optimizer(const IntegerVector alpha,
               const IntegerMatrix beta,
               const unsigned gamma)
{
    RNGScope scope; // needed for debugging
    // get dimensions of matrix
    const unsigned n_row = gamma;
    const unsigned n_col = alpha.length();

    // 0. create a solution based on alpha and gamma diversity
    IntegerMatrix current_solution = gen_init_solution(alpha,
                                                       gamma);

    // 1. Find the most efficient order, sorted by alpha diversity
    std::vector<unsigned> alpha_ordered = as<std::vector<unsigned> >(alpha);
    std::sort(alpha_ordered.begin(), alpha_ordered.end()); // ascending
    std::reverse(alpha_ordered.begin(), alpha_ordered.end()); // now descending order

    // 2. Number of possible combinations
    unsigned combinations = factorial(gamma) / (
                factorial(alpha[alpha_ordered[0]]) *
            factorial(gamma - alpha[alpha_ordered[0]])
            );

}

unsigned long factorial_loop(unsigned n)
{
    unsigned fact = 1;
    if (n == 1 || n == 0) {
        return fact;
    } else {
        for(unsigned i = 1; i <= n; i++) {
            fact *= i;
        }
    }
    return fact;
}

unsigned long factorial_ln(unsigned n)
{

}

unsigned long factorial_rec(unsigned n)
{
    return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

unsigned long factorial(unsigned n)
{
    const std::vector<unsigned long> factorials =
    {1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880,  3628800,  39916800,
     479001600,  6227020800,  87178291200, 1307674368000,  20922789888000,
     355687428096000,  6402373705728000, 121645100408832000,  2432902008176640000};

    if (n < factorials.size()) {
        return factorials[n];
    }

    return factorial(n);
}
