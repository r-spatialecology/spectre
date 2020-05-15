//#include "mh_optimizer.h"
//#include <random>
//#include <chrono>
//#include "calculate_solution_commonness.h"
//#include "rcpp_sample.h"
////omp_set_nested(1);

//List mh_optimizer(IntegerVector alpha_list,
//                  const unsigned total_gamma,
//                  IntegerMatrix solution_commonness_target,
//                  const unsigned max_iterations,
//                  unsigned long seed,
//                  bool verbose,
//                  double base_probability_jump)
//{
//    RNGScope scope; // needed for debugging
//    // get dimensions of matrix
//    const unsigned n_row = total_gamma;
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
//    IntegerMatrix current_solution = gen_init_solution(n_row, n_col, alpha_list);
//    IntegerMatrix best_solution = current_solution;

//    // Random number generator
//    if (seed < 1) {
//        seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
//    }
//    std::mt19937 rng(seed);
//    std::uniform_int_distribution<unsigned> site_dist(0, current_solution.ncol() - 1);
//    std::uniform_real_distribution<double> real_dist_01(0, 1);

//    // calculate the site x site commonness for the current solution
//    auto solution_commonness = calculate_solution_commonness_rcpp(current_solution);

//    // calculate the difference between target and current solution
//    double energy = calc_energy(solution_commonness, solution_commonness_target);

//    energy_vector.push_back(energy);
//    acceptance_rate.push_back(0.0);
//    iteration_count.push_back(0);
//    double accepted = 0;

//    for (unsigned iter = 0; iter < max_iterations; iter++) {

//        // swap species occurrence at a random site
//        const unsigned random_site = site_dist(rng);
//        IntegerVector species_to_swap = get_swap_rows_rcpp(
//                    current_solution(_, random_site), true);

//        std::swap(current_solution(species_to_swap[0], random_site),
//                current_solution(species_to_swap[1], random_site));

//        // calculate commoness and energy for new solution
//        const auto new_solution_commonness = calculate_solution_commonness_site_rcpp(
//                    current_solution, solution_commonness, random_site + 1);

//        const double new_energy = calc_energy(new_solution_commonness, solution_commonness_target);

//        // if proposed solution is better or same, accept, else calculate jump probability
//        if (new_energy <= energy) {
//            best_solution = current_solution;
//            energy = new_energy;
//        } else if (jump_probability(new_energy, energy, base_probability_jump) > real_dist_01(rng)) {
//            // accept anyway?
//            energy = new_energy;
//        } else {
//            // not accepted, undo changes
//            std::swap(current_solution(species_to_swap[0], random_site),
//                    current_solution(species_to_swap[1], random_site));
//        }
////        update_metrics(accepted, energy, new_energy, iter, energy_vector,
////                       iteration_count, acceptance_rate);
//        // print progress
//        if (verbose) {
//            Rcout << " Progress: max_runs: " << iter + 1 << "/" << max_iterations <<
//                     " || energy = " << energy << " %\r";
//            Rcout.flush();
//        }
//        update_metrics(accepted, energy, iter, energy_vector,
//                       iteration_count, acceptance_rate);
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

//double calc_energy(const IntegerVector solution_commonness, const IntegerVector solution_commonness_target)
//{
//    // calculate the difference between target and current solution
//    return(sum(abs(na_omit(solution_commonness_target - solution_commonness))) /
//           (double)sum(na_omit(solution_commonness_target)));
//}

//void species_swap_rcpp(IntegerMatrix &mat, const IntegerVector species, unsigned site) {
//    std::swap(mat(species[0] - 1, site - 1), mat(species[1] - 1, site - 1));
//}

//IntegerMatrix gen_init_solution(const unsigned n_row, const unsigned n_col, const IntegerVector &alpha_list)
//{
//    IntegerMatrix current_solution = IntegerMatrix(n_row, n_col);

//    const IntegerVector species_seq = seq_len(n_row);
//    for (unsigned site = 0; site < n_col; site++) {
//        const IntegerVector change_species = sample(species_seq, alpha_list[site], false);
//        for (auto species : change_species) {
//            current_solution(species, site) = 1;
//        }
//    }
//    return current_solution;
//}

//void update_metrics(double &accepted, const double energy, const unsigned iter,
//                    std::vector<double> &energy_vector,
//                    std::vector<unsigned> &iteration_count,
//                    std::vector<double> &acceptance_rate)
//{
//    accepted++;
//    energy_vector.push_back(energy);
//    iteration_count.push_back(iter + 1);
//    acceptance_rate.push_back(accepted / iter);
//}



//double jump_probability(const double new_energy, const double energy, double base_prob)
//{
//    double result = 1.0;
//    if (new_energy > energy) {
//        result = energy / new_energy;
//    }

//    return result - (1.0 - base_prob);
//}

//List mh_optimizer_neutral(IntegerVector alpha_list,
//                          const unsigned total_gamma,
//                          IntegerMatrix solution_commonness_target,
//                          const unsigned max_iterations,
//                          unsigned long seed,
//                          bool verbose)
//{
//    RNGScope scope; // needed for debugging
//    // get dimensions of matrix
//    const unsigned n_row = total_gamma;
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
//    IntegerMatrix current_solution = gen_init_solution(n_row, n_col, alpha_list);
//    IntegerMatrix best_solution = current_solution;

//    // Random number generator
//    if (seed < 1) {
//        seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
//    }
//    std::mt19937 rng(seed);
//    std::uniform_int_distribution<unsigned> site_dist(0, current_solution.ncol() - 1);
//    std::uniform_real_distribution<double> real_dist_01(0, 1);

//    // calculate the site x site commonness for the current solution
//    auto solution_commonness = calculate_solution_commonness_rcpp(current_solution);

//    // calculate the difference between target and current solution
//    double energy = calc_energy(solution_commonness, solution_commonness_target);

//    energy_vector.push_back(energy);
//    acceptance_rate.push_back(0.0);
//    iteration_count.push_back(0);
//    double accepted = 0;

//    for (unsigned iter = 0; iter < max_iterations; iter++) {

//        // new solution
//        current_solution = gen_init_solution(n_row, n_col, alpha_list);
//        // calculate commoness and energy for new solution
//        const auto new_solution_commonness = calculate_solution_commonness_rcpp(current_solution);

//        const double new_energy = calc_energy(new_solution_commonness, solution_commonness_target);

//        // if proposed solution is better or same, accept, else calculate jump probability
//        if (new_energy < energy) {
//            best_solution = current_solution;
//            energy = new_energy;
//        }
//        update_metrics(accepted, energy, iter, energy_vector,
//                       iteration_count, acceptance_rate);
//        // print progress
//        if (verbose) {
//            Rcout << " Progress: max_runs: " << iter + 1 << "/" << max_iterations <<
//                     " || energy = " << energy << " %\r";
//            Rcout.flush();
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