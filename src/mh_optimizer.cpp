#include "mh_optimizer.h"
#include <random>
#include <chrono>
#include "calculate_solution_commonness.h"
#include "rcpp_sample.h"

List mh_optimizer(IntegerVector alpha_list,
                  const int total_gamma,
                  IntegerMatrix target,
                  const double acceptance_rate_threshold,
                  const unsigned max_iterations,
                  const unsigned burn_in, unsigned long seed)
{
    // get dimensions of matrix
    const int n_row = total_gamma;
    const int n_col = alpha_list.length();

    // check for sanity to avoid an infinite loop later on
    const int check = sum(alpha_list);
    assert(check > 0);
    assert(check < total_gamma * n_col);

    std::vector<unsigned> energy_vector;
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
    std::uniform_int_distribution<unsigned> species_dist(0, current_solution.nrow() - 1);
    std::uniform_real_distribution<double> real_dist_01(0, 1);

    // loop changes the alpha number of species in each site to present (i.e. 1)
    const IntegerVector species_seq = seq_len(n_row);
    for (int site = 0; site < n_col; site++) {
        const IntegerVector change_species = rcpp_sample(species_seq, alpha_list[site]);
        for (auto species : change_species) {
            current_solution(species, site) = 1;
        }
    }

    // calculate the site x site commonness for the current solution
    auto solution_commonness = calculate_solution_commonness_rcpp(current_solution);

    // calculate the difference between target and current solution
    unsigned energy = abs(sum(na_omit(solution_commonness - target)));
    energy_vector.push_back(energy);
    iteration_count.push_back(0);

    unsigned accepted = 0;

    for (unsigned iter = 0; iter < max_iterations + burn_in; ++iter) {
        // swap species occurrence in a random site
        const unsigned site = site_dist(rng);
        const unsigned species1 = species_dist(rng);
        const unsigned species2 = species_dist(rng);

        // check if species1 and species2 are neither both present nor both absent
        // if one species is present (=1) and the other not (=0),
        // the sum of both is 1 which C++ interprets as true:
        if (!(current_solution(species1, site) + current_solution(species2, site))) {
            iter--;
            continue; // bad luck, try again...
        }

        // change values
        std::swap(current_solution(site, species1), current_solution(site, species2));

        // calculate commoness for both solutions
        ///TODO: use site function
        const auto new_solution_commonness =
                //calculate_solution_commonness_site_rcpp(current_solution,)
               calculate_solution_commonness_rcpp(current_solution);
        const unsigned new_energy = abs(sum(na_omit(new_solution_commonness - target)));

        // if proposed solution is better, accept, else calculate jump probability
        if (new_energy <= energy) {
            if (iter > burn_in) {
                accepted++;
                energy = new_energy;
                energy_vector.push_back(energy);
                iteration_count.push_back(iter);
                acceptance_rate.push_back(accepted / iter);
                if (acceptance_rate.back() < acceptance_rate_threshold) {
                    break;
                }
            }
            continue;
        }

        // accept anyway?
        const double jump_probability = energy / new_energy;
        if (jump_probability > real_dist_01(rng)) {
            if (iter > burn_in) {
                accepted++;
                energy = new_energy;
                energy_vector.push_back(energy);
                iteration_count.push_back(iter);
                acceptance_rate.push_back(accepted / iter);
                if (acceptance_rate.back() < acceptance_rate_threshold) {
                    break;
                }
            }
            continue;
        }

        // not accepted, undo changes
        std::swap(current_solution(site, species1), current_solution(site, species2));
    }

    // generate dataframe for i and energy
    DataFrame measures_df = DataFrame::create(_["i"] = iteration_count,
            _["energy"] = energy, _["acceptance_rate"] = acceptance_rate);

    List results = List::create(current_solution, measures_df);

    return(results);
}
