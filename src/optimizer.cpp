#include "optimizer.h"
#include "minconf.h"

List optimizer_min_conf(const IntegerVector alpha_list,
                        const unsigned total_gamma, const IntegerMatrix target,
                        const unsigned max_iterations,
                        const IntegerMatrix partial_solution,
                        const IntegerMatrix fixed_species,
                        const unsigned long seed, const bool verbose,
                        const bool interruptible) {
  MinConf mc(as<std::vector<unsigned>>(alpha_list), total_gamma,
             as<std::vector<int>>(target),
             as<std::vector<int>>(partial_solution),
             as<std::vector<int>>(fixed_species), seed, NA_INTEGER);

  const long iter =
      max_iterations - mc.optimize(max_iterations, verbose, interruptible);

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
    float best_error =
        *std::min_element(mc.error_vector.begin(), mc.error_vector.end());
    float worst_error =
        *std::max_element(mc.error_vector.begin(), mc.error_vector.end());
    Rcout << "\n > Optimization finished with lowest absolute error = "
          << best_error << " (highest absolute error was: " << worst_error
          << " improved by: " << worst_error - best_error << ")";
  }

  return (results);
}

IntegerMatrix
calculate_solution_commonness_rcpp(const IntegerMatrix solution_matrix) {
  // create results matrix, filled with NAs
  // and translate solution_matrix to 2-D std::vector
  const unsigned n_sites = solution_matrix.ncol();
  const unsigned gamma_div = solution_matrix.nrow();
  std::vector<int> fixed_mat(n_sites * gamma_div, 1); // all fixed

  MinConf mc(std::vector<unsigned>(n_sites), gamma_div,
             std::vector<int>(n_sites * n_sites),
             as<std::vector<int>>(solution_matrix), fixed_mat, 0, NA_INTEGER);

  // zero iterations, just to calculate the commonness
  mc.optimize(0, false, false);

  IntegerMatrix result(n_sites, n_sites);
  std::fill(result.begin(), result.end(), NA_INTEGER);

  for (unsigned site = 0; site < n_sites; site++) {
    for (unsigned other_site = 0; other_site < n_sites; other_site++) {
      if (mc.commonness[site][other_site] != NA_INTEGER) {
        result(site, other_site) = mc.commonness[site][other_site];
      }
    }
  }

  return result;
}

NumericMatrix
calculate_solution_bc_rcpp(const IntegerMatrix solution_matrix, IntegerVector alpha_vector) {
    auto commonness = calculate_solution_commonness_rcpp(solution_matrix);
    const unsigned n_sites = solution_matrix.ncol();
    const unsigned gamma_div = solution_matrix.nrow();
    std::vector<int> fixed_mat(n_sites * gamma_div, 1); // all fixed

    MinConf mc(as<std::vector<unsigned>>(alpha_vector), gamma_div,
               std::vector<int>(n_sites * n_sites),
               as<std::vector<int>>(solution_matrix), fixed_mat, 0, NA_INTEGER);

    // zero iterations, just to calculate the commonness
    mc.optimize(0, false, false);

    NumericMatrix result(n_sites, n_sites);
    std::fill(result.begin(), result.end(), NA_REAL);

    for (unsigned site = 0; site < n_sites; site++) {
        for (unsigned other_site = 0; other_site < n_sites; other_site++) {
            if (mc.bray_curtis[site][other_site] != NA_INTEGER) {
                result(site, other_site) = mc.bray_curtis[site][other_site];
            }
        }
    }

    return result;
}
