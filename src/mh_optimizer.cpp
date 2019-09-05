#include "mh_optimizer.h"
#include <random>
#include <chrono>
#include "calculate_solution_commonness.h"
#include "rcpp_sample.h"
//omp_set_nested(1);

List mh_optimizer(IntegerVector alpha_list,
                  const unsigned total_gamma,
                  IntegerMatrix target,
                  const double acceptance_rate_threshold,
                  const unsigned max_iterations,
                  const unsigned burn_in, unsigned long seed)
{
    RNGScope scope; // needed for debugging
    // get dimensions of matrix
    const unsigned n_row = total_gamma;
    const unsigned n_col = alpha_list.length();

    // check for sanity to avoid an infinite loop later on
    const unsigned check = sum(alpha_list);
    assert(check > 0);
    assert(check < total_gamma * n_col);

    std::vector<double> energy_vector;
    energy_vector.reserve(max_iterations);
    std::vector<unsigned> iteration_count;
    iteration_count.reserve(max_iterations);
    std::vector<double> acceptance_rate;
    acceptance_rate.reserve(max_iterations);

    // generate initial solution
    ///TODO: do this in parallel (chains)
    IntegerMatrix current_solution = IntegerMatrix(n_row, n_col);

    if (seed < 1) {
        seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    }

    std::mt19937 rng(seed);
    std::uniform_int_distribution<unsigned> site_dist(0, current_solution.ncol() - 1);
    //std::uniform_int_distribution<unsigned> species_dist(0, current_solution.nrow() - 1);
    std::uniform_real_distribution<double> real_dist_01(0, 1);

    // loop changes the alpha number of species in each site to present (i.e. 1)
    const IntegerVector species_seq = seq_len(n_row);
    for (unsigned site = 0; site < n_col; site++) {
        const IntegerVector change_species = sample(species_seq, alpha_list[site], false);
        for (auto species : change_species) {
            current_solution(species, site) = 1;
        }
    }

    // calculate the site x site commonness for the current solution
    auto solution_commonness = calculate_solution_commonness_rcpp(current_solution);

    // calculate the difference between target and current solution
    double energy = calc_energy(target, solution_commonness);
    energy_vector.push_back(energy);
    iteration_count.push_back(burn_in);
    acceptance_rate.push_back(0.0);

    double accepted = 0;

    for (unsigned iter = 0; iter < max_iterations + burn_in; iter++) {
        // swap species occurrence in a random site
        const unsigned site = site_dist(rng);
        auto species = get_swap_rows_rcpp(current_solution(_, site));

        // swap species
        std::swap(current_solution(species[0], site), current_solution(species[1], site));
///BUG: swapping seems not to work!
        // calculate commoness for both solutions
        auto new_solution_commonness = calculate_solution_commonness_site_rcpp(current_solution, solution_commonness, site + 1);
       // update_solution_commonness_site_rcpp(current_solution, new_solution_commonness, site + 1);

        const double new_energy = calc_energy(target, new_solution_commonness);

        // if proposed solution is better or same, accept, else calculate jump probability
        if (new_energy <= energy) {
            if (iter > burn_in) {
                accepted++;
                energy = new_energy;
                energy_vector.push_back(energy);
                iteration_count.push_back(iter);
                acceptance_rate.push_back(accepted / iter);
                //                if (acceptance_rate.back() < acceptance_rate_threshold) {
                //                    break;
                //                }
            }
            continue;
        }

        // accept anyway?
        const double jump_probability = energy; ///TODO: adjust
        if (jump_probability > real_dist_01(rng)) {
            if (iter > burn_in) {
                accepted++;
                energy = new_energy;
                energy_vector.push_back(energy);
                iteration_count.push_back(iter);
                acceptance_rate.push_back(accepted / iter);
                //                if (acceptance_rate.back() < acceptance_rate_threshold) {
                //                    break;
                //                }
            }
            continue;
        }

        // not accepted, undo changes
        std::swap(current_solution(species[0], site), current_solution(species[1], site));
    }

    // generate dataframe for i and energy
    DataFrame measures_df = DataFrame::create(_["i"] = iteration_count,
            _["energy"] = energy, _["acceptance_rate"] = acceptance_rate);

    List results = List::create(current_solution, measures_df);
    //std::cout << "Optimization finished with an energy = " << energy << std::endl;
    return(results);
}

double calc_energy(const IntegerVector target, const IntegerVector solution_commonness)
{
    // calculate the difference between target and current solution
    return(sum(abs(na_omit(target - solution_commonness))) /
            (double)sum(na_omit(target)));
}

void species_swap_rcpp(IntegerMatrix &mat, const IntegerVector species, unsigned site) {
    std::swap(mat(species[0] - 1, site - 1), mat(species[1] - 1, site - 1));
}
