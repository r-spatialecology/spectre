#include "optimizer.h"
#include <random>
#include <chrono>
#include "calculate_solution_commonness.h"
#include "rcpp_sample.h"
//omp_set_nested(1);

List optimizer(IntegerVector alpha_list,
               const unsigned total_gamma,
               NumericMatrix solution_commonness_target,
               const unsigned max_iterations,
               unsigned long seed,
               bool verbose,
               const double increment)
{
    RNGScope scope; // needed for debugging
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

    // generate initial solution, loop changes the alpha number of species in each site to present (i.e. 1)
    ///TODO: do this in parallel (chains)
    NumericMatrix current_solution = gen_init_solution(alpha_list, total_gamma);
    NumericMatrix best_solution = current_solution;
    auto best_solution_ = as<std::vector<double> >(current_solution);
    std::vector<std::vector<unsigned> > considered_species(alpha_list.size(),
                                                           std::vector<unsigned>(total_gamma));
    for (unsigned site = 0; site < alpha_list.size(); site++) {
        for (unsigned species = 0; species < total_gamma; species++) {
            considered_species[site][species] = species;
        }
    }

    // Random number generator
    if (seed < 1) {
        seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    }
    std::mt19937 rng(seed);
    std::uniform_int_distribution<unsigned> random_site(0, alpha_list.size() - 1);
    std::uniform_real_distribution<double> random_double(0, 1);

    // calculate the site x site commonness for the current solution
    auto solution_commonness = calculate_solution_commonness_rcpp(current_solution);
    auto target_solution_ = as<std::vector<double> >(solution_commonness_target);
    auto solution_commonness_ = as<std::vector<double> >(solution_commonness);

    // calculate the difference between target and current solution
    double energy = calc_energy(solution_commonness, solution_commonness_target);

    energy_vector.push_back(energy);
    acceptance_rate.push_back(0.0);
    iteration_count.push_back(0);
    double accepted = 0;

    for (unsigned iter = 0; iter < max_iterations; iter++) {
        const unsigned site = random_site(rng);
        if(considered_species[site].empty())
            continue;

        std::uniform_int_distribution<unsigned>
                random_species(0, considered_species[site].size() - 1);
        const unsigned idx = random_species(rng);
        const unsigned cur_species = considered_species[site][idx];

        best_solution_ = as<std::vector<double> >(current_solution);

        double unassigned_incr = 0.0;

        current_solution[site * total_gamma + cur_species] -= 2 * increment;

        for (unsigned species : considered_species[site]) {
            current_solution[site * total_gamma + cur_species] += increment;

            if(current_solution[site * total_gamma + species] > 0.95 ||
                    current_solution[site * total_gamma + species] < 0.05) {
                const double val_ = current_solution[site * total_gamma + species];
                const double new_val_ =
                        round(current_solution[site * total_gamma + species]);
                unassigned_incr +=
                        current_solution[site * total_gamma + species] -
                        new_val_;
                current_solution[site * total_gamma + species] = new_val_;

                // remove species from considered_species
                considered_species[site].erase(std::remove(
                                                   considered_species[site].begin(),
                                                   considered_species[site].end(), species),
                                               considered_species[site].end());
            }
            auto current_solution_ = as<std::vector<double> >(current_solution);
            auto val = current_solution[site * total_gamma + species];


            if (unassigned_incr > 0.0) {
                for (unsigned species : considered_species[site]) {
                    current_solution[site * total_gamma + cur_species] += unassigned_incr;
                }
            }

            best_solution_ = as<std::vector<double> >(current_solution);
            // calculate commoness and energy for new solution
            const auto new_solution_commonness = calculate_solution_commonness_site_rcpp(
                        current_solution, solution_commonness, site + 1);
            const auto new_solution_commonness_ = as<std::vector<double> >(new_solution_commonness);

            const double new_energy = calc_energy(new_solution_commonness, solution_commonness_target);

            // if proposed solution is better or same, accept, else calculate jump probability
            if (new_energy <= energy) {
                best_solution = current_solution;
                best_solution_ = as<std::vector<double> >(current_solution);
                solution_commonness = new_solution_commonness;
                solution_commonness_ = as<std::vector<double> >(solution_commonness);
                energy = new_energy;
                //        } else if (jump_probability(new_energy, energy, base_probability_jump) > random_double(rng)) {
                //            // accept anyway?
                //            energy = new_energy;
                update_metrics(accepted, energy, iter, energy_vector,
                               iteration_count, acceptance_rate);
            } else {
                // not accepted, undo changes
                current_solution = best_solution;
            }
        }
        //        update_metrics(accepted, energy, new_energy, iter, energy_vector,
        //                       iteration_count, acceptance_rate);
        // print progress
        if (verbose) {
            Rcout << " Progress: max_runs: " << iter + 1 << "/" << max_iterations <<
                     " || energy = " << energy << " %\r";
            Rcout.flush();
        }
        //        update_metrics(accepted, energy, iter, energy_vector,
        //                       iteration_count, acceptance_rate);
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


//List optimizer(IntegerVector alpha_list,
//               const unsigned total_gamma,
//               NumericMatrix solution_commonness_target,
//               const unsigned max_iterations,
//               unsigned long seed,
//               bool verbose,
//               double base_probability_jump)
//{
//    RNGScope scope; // needed for debugging
//    // get dimensions of matrix
//    const unsigned n_col = alpha_list.length();

//    // check for sanity to avoid an infinite loop later on
//    // assert() is actually bad practice in Rcpp because it also kills R...
//    const unsigned check = sum(alpha_list);
//    assert(check > 0);
//    assert(check < total_gamma * n_col);

//    std::vector<double> energy_vector;
//    energy_vector.reserve(max_iterations + 1);
//    std::vector<unsigned> iteration_count;
//    iteration_count.reserve(max_iterations + 1);
//    std::vector<double> acceptance_rate;
//    acceptance_rate.reserve(max_iterations + 1);

//    // generate initial solution, loop changes the alpha number of species in each site to present (i.e. 1)
//    ///TODO: do this in parallel (chains)
//    NumericMatrix current_solution = gen_init_solution(alpha_list, total_gamma);
//    NumericMatrix best_solution = current_solution;
//    auto best_solution_ = as<std::vector<double> >(current_solution);
//    std::vector<std::vector<unsigned> > considered_species(alpha_list.size(),
//                                                           std::vector<unsigned>(total_gamma));
//    for (unsigned site = 0; site < alpha_list.size(); site++) {
//        for (unsigned species = 0; species < total_gamma; species++) {
//            considered_species[site][species] = species;
//        }
//    }

//    // Random number generator
//    if (seed < 1) {
//        seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
//    }
//    std::mt19937 rng(seed);
//    std::uniform_int_distribution<unsigned> random_site(0, alpha_list.size() - 1);
//    std::uniform_real_distribution<double> random_double(0, 1);

//    // calculate the site x site commonness for the current solution
//    auto solution_commonness = calculate_solution_commonness_rcpp(current_solution);
//    auto target_solution_ = as<std::vector<double> >(solution_commonness_target);
//    auto solution_commonness_ = as<std::vector<double> >(solution_commonness);

//    // calculate the difference between target and current solution
//    double energy = calc_energy(solution_commonness, solution_commonness_target);

//    energy_vector.push_back(energy);
//    acceptance_rate.push_back(0.0);
//    iteration_count.push_back(0);
//    double accepted = 0;

//    for (unsigned iter = 0; iter < max_iterations; iter++) {
//        for  (unsigned site = 0; site < alpha_list.size(); site++) { ///TODO: Make random

//            //  const unsigned site = random_site(rng);
//            if(considered_species[site].empty())
//                continue;

//            std::uniform_int_distribution<unsigned>
//                    random_species(0, considered_species[site].size() - 1);
//            const unsigned idx = random_species(rng);
//            const unsigned cur_species = considered_species[site][idx];

//            best_solution_ = as<std::vector<double> >(current_solution);
//            const double new_val = round(current_solution[site * total_gamma + cur_species]);
//            const double incr =
//                    (current_solution[site * total_gamma + cur_species] -
//                    new_val) / (considered_species[site].size() - 1);
//            double unassigned_incr = 0.0;

//            current_solution[site * total_gamma + cur_species] = new_val;

//            for (unsigned species : considered_species[site]) {
//                if(current_solution[site * total_gamma + species] > 0.95 ||
//                        current_solution[site * total_gamma + species] < 0.05) {
//                    const double new_val_ =
//                            round(current_solution[site * total_gamma + species]);
//                    const double val_ = current_solution[site * total_gamma + species];
//                    unassigned_incr +=
//                            current_solution[site * total_gamma + species] -
//                            new_val_;
//                    current_solution[site * total_gamma + species] = new_val_;

//                    // remove species from considered_species
////                    considered_species[site].erase(std::remove(
////                                                       considered_species[site].begin(),
////                                                       considered_species[site].end(), species),
////                                                   considered_species[site].end());
//                } else {
//                    current_solution[site * total_gamma + species] += incr;
//                }
//                auto current_solution_ = as<std::vector<double> >(current_solution);
//                auto val = current_solution[site * total_gamma + species];
//            }

//            for (unsigned species : considered_species[site]) {
//               current_solution[site * total_gamma + cur_species] += unassigned_incr;
//            }

//            best_solution_ = as<std::vector<double> >(current_solution);
//            // calculate commoness and energy for new solution
//            const auto new_solution_commonness = calculate_solution_commonness_site_rcpp(
//                        current_solution, solution_commonness, site + 1);
//            const auto new_solution_commonness_ = as<std::vector<double> >(new_solution_commonness);

//            const double new_energy = calc_energy(new_solution_commonness, solution_commonness_target);

//            // if proposed solution is better or same, accept, else calculate jump probability
//            if (new_energy <= energy) {
//                best_solution = current_solution;
//                best_solution_ = as<std::vector<double> >(current_solution);
//                solution_commonness = new_solution_commonness;
//                solution_commonness_ = as<std::vector<double> >(solution_commonness);
//                energy = new_energy;
//                //        } else if (jump_probability(new_energy, energy, base_probability_jump) > random_double(rng)) {
//                //            // accept anyway?
//                //            energy = new_energy;
//                update_metrics(accepted, energy, iter, energy_vector,
//                               iteration_count, acceptance_rate);
//            } else {
//                // not accepted, undo changes
//                current_solution = best_solution;
//            }
//            //        update_metrics(accepted, energy, new_energy, iter, energy_vector,
//            //                       iteration_count, acceptance_rate);
//            // print progress
//            if (verbose) {
//                Rcout << " Progress: max_runs: " << iter + 1 << "/" << max_iterations <<
//                         " || energy = " << energy << " %\r";
//                Rcout.flush();
//            }
//            //        update_metrics(accepted, energy, iter, energy_vector,
//            //                       iteration_count, acceptance_rate);
//        }
//    }

//    // generate dataframe for i and energy
//    DataFrame measures_df = DataFrame::create(_["i"] = iteration_count,
//            _["energy"] = energy_vector);

//    List results = List::create(Rcpp::Named("optimized_grid") = best_solution,
//                                Rcpp::Named("energy") = measures_df);
//    if (verbose) {
//        double best_energy = *std::min_element(energy_vector.begin(), energy_vector.end());
//        double worst_energy = *std::max_element(energy_vector.begin(), energy_vector.end());
//        Rcout << "\n > Optimization finished with lowest energy = " << best_energy << " %"
//              << " (highest energy was: " << worst_energy << " %, improved by: "
//              << worst_energy - best_energy << " %)";
//    }
//    return(results);
//}

double calc_energy(const NumericMatrix solution_commonness,
                   const NumericMatrix solution_commonness_target)
{
    // calculate the difference between target and current solution
    // return(sum(abs(na_omit(solution_commonness_target - solution_commonness))));
    return(sum((na_omit(solution_commonness_target - solution_commonness)) *
               (na_omit(solution_commonness_target - solution_commonness))));
}

NumericMatrix gen_init_solution(const IntegerVector &alpha_list,
                                const unsigned gamma_diversity)
{
    const auto n_sites = alpha_list.size(); // number of cols

    NumericMatrix current_solution(gamma_diversity, n_sites);

    for (int site = 0; site < n_sites; site++) {
        NumericVector proportion(gamma_diversity,
                                 static_cast<double>(alpha_list[site]) /
                                 gamma_diversity);

        current_solution(_, site) = proportion;
    }

    return current_solution;
}

void update_metrics(double &accepted, const double energy, const unsigned iter,
                    std::vector<double> &energy_vector,
                    std::vector<unsigned> &iteration_count,
                    std::vector<double> &acceptance_rate)
{
    accepted++;
    energy_vector.push_back(energy);
    iteration_count.push_back(iter + 1);
    acceptance_rate.push_back(accepted / iter);
}


