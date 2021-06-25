#include "minconf.h"
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>

/**
 * @brief MinConf
 * @param alpha_list list of alpha diversity for each site
 * @param gamma_div total number of species
 * @param target matrix of the target beta diversity
 * @param partial_solution optional solution matrix as a start for the
 * optimization. Species are added to fit the alpha diveristy of each site
 * @param fixed_species matrix to fix specific species presences/absences (0:
 * not fixed, !=0 fixed)
 * @param seed random seed
 * @param na_val NA value (default: Rcpp's NA -2147483648)
 */
MinConf::MinConf(const std::vector<unsigned> &alpha_list,
                 const unsigned gamma_div, const std::vector<int> &target,
                 const std::vector<int> &partial_solution,
                 const std::vector<int> &fixed_species,
                 const unsigned long seed, const int na_val)
    : NA(na_val), alpha_list(alpha_list), gamma_div(gamma_div),
      n_sites(alpha_list.size()) {
  // Random number generator
  rng = std::mt19937(seed);

  solution.resize(n_sites);
  this->target.resize(n_sites);
  commonness.resize(n_sites);

  for (unsigned site = 0; site < n_sites; site++) {
    solution[site].resize(gamma_div);
    commonness[site].resize(n_sites, NA);

    // convert target matrix to a more convenient format
    this->target[site].resize(n_sites, NA);
    // iterate over the upper triagonal matrix, only (w/o diagonal)
    for (unsigned other_site = site + 1; other_site < n_sites; other_site++) {
      if (target[other_site * n_sites + site] == NA) {
        this->target[site][other_site] = NA;
      } else {
        this->target[site][other_site] = target[other_site * n_sites + site];
      }
    }
  }

  if (partial_solution.size() > 1) {
    if (partial_solution.size() != n_sites * gamma_div) {
      Rcpp::stop("The size of the partial_solution vector does not match "
                 "n_sites * gamma_div. "
                 "partial_solution ignored.");
    } else {
      for (unsigned site = 0; site < n_sites; site++) {
        for (unsigned species = 0; species < gamma_div; species++) {
          if (partial_solution[site * gamma_div + species]) {
            solution[site][species] = 1;
          }
        }
      }
    }
  }

  if (fixed_species.size() > 1) {
    if (fixed_species.size() != n_sites * gamma_div) {
      Rcpp::stop("The size of the fixed_species vector does not match "
                 "n_sites * gamma_div. "
                 "fixed_species ignored.");
    } else {
      this->fixed_species.resize(n_sites);
      for (unsigned site = 0; site < n_sites; site++) {
        this->fixed_species[site].resize(gamma_div);
        for (unsigned species = 0; species < gamma_div; species++) {
          if (fixed_species[site * gamma_div + species]) {
            this->fixed_species[site][species] = 1;
          }
        }
      }
    }
  }

  gen_init_solution();
}

/**
 * @brief optimize runs the actual optimization algorithm
 * @param max_steps
 * @param verbose enables a progress bar to the R terminal. This will slow
 * down the optimization.
 * @param interruptible makes the optimization interruptable from R (this is
 * computationally a rather cheap option)
 * @return
 */
int MinConf::optimize(const long max_steps, bool verbose, bool interruptible) {
  Progress p(max_steps, verbose);
  std::uniform_int_distribution<unsigned> site_dist(0, n_sites - 1);

  auto iter = max_steps;

  // calculate and save the start values
  update_solution_commonness();
  unsigned error = calc_error();
  iteration_count.push_back(0);
  error_vector.push_back(error);

  // optimize
  while (iter-- > 0) {
    // update progress bar
    p.increment();
    if (interruptible) {
      if (Progress::check_abort()) {
        return RET_ABORT;
      }
    }

    // choose random site
    const auto site = site_dist(rng);

    // remove a random species at this site
    // (and skip iteration if it was not possible to remove a species)
    if (remove_random_species(site) == false) {
      continue;
    }

    // add min conf species
    add_species_min_conf(site);

    update_solution_commonness();
    error = calc_error();

    iteration_count.push_back(max_steps - iter);
    error_vector.push_back(error);

    if (error == 0) {
      return iter;
    }
  }
  return iter;
}

/**
 * @brief MinConf::calc_missing_species calculates the number of species the
 * algorithm needs to add to equal the number of species with alpha diversity.
 * The actual number of missing species species might differ for two reasons: If
 * the number of species present is higher than the number of species estimated
 * by alpha diversity (due to the partial_solution parameter), then the number
 * of missing species is 0. If the number of missing species is higher than the
 * number of potential species (because potential species are blocked by the
 * fixed_species parameter), then the number of missing_species is equal to the
 * number of potential species.
 * @return vector of the number of species to be assigned for each site.
 */
std::vector<unsigned> MinConf::calc_missing_species() {
  auto missing_species = alpha_list;
  for (unsigned site = 0; site < n_sites; site++) {
    const auto present_species = present_species_index(site, false).size();
    const auto potential_species = absent_species_index(site).size();
    missing_species[site] = (missing_species[site] < present_species)
                                ? 0
                                : missing_species[site] - present_species;

    if (missing_species[site] > potential_species) {
      missing_species[site] = potential_species;
    }
  }

  return missing_species;
}

/**
 * @brief MinConf::present_species_index
 * @param site
 * @param omit_fixed_species
 * @return a vector with the indices of each species present.
 */
std::vector<unsigned> MinConf::present_species_index(unsigned site,
                                                     bool omit_fixed_species) {
  std::vector<unsigned> species_idx;
  for (unsigned species = 0; species < gamma_div; species++) {
    if (omit_fixed_species && fixed_species.size()) {
      if (fixed_species[site][species]) {
        continue;
      }
    }
    if (solution[site][species]) {
      species_idx.push_back(species);
    }
  }
  return species_idx;
}

/**
 * @brief MinConf::absent_species_index omits fixed species
 * @param site
 * @return a vector with the indices of each species absent.
 */
std::vector<unsigned> MinConf::absent_species_index(unsigned site) {
  std::vector<unsigned> absent_species_idx;
  for (unsigned species = 0; species < gamma_div; species++) {
    if (!solution[site][species]) {
      if (fixed_species.size()) {
        if (fixed_species[site][species]) {
          continue;
        }
      }
      absent_species_idx.push_back(species);
    }
  }
  return absent_species_idx;
}

/**
 * @brief MinConf::remove_random_species
 * @param site
 * @return true on success
 */
bool MinConf::remove_random_species(const unsigned site) {
  std::vector<unsigned> species_idx =
      present_species_index(site, true); // omit fixed_species == true
  if (species_idx.size() == 0) {
    return false; // all species fixed, nothing to remove here
  }
  std::shuffle(species_idx.begin(), species_idx.end(), rng);
  const unsigned species = species_idx.back();
  solution[site][species] = 0;

  return true;
}

/**
 * @brief MinConf::add_species_min_conf adds a species at `site` that results in
 * the lowest error. If there is more than one possibility, a random choice is
 * drawn.
 * @param site
 */
void MinConf::add_species_min_conf(unsigned site) {
  // calculate the best-fitting species for this site
  auto min_conflict_species = calc_min_conflict_species(site);

  if (min_conflict_species.size() < 1) {
    Rcpp::Rcerr << "no species found to add at add_species_min_conf, site: "
                << site << std::endl; // something odd happened
  } else {
    std::shuffle(min_conflict_species.begin(), min_conflict_species.end(), rng);
    const auto species = min_conflict_species[0];
    solution[site][species] = 1;
  }
}

/**
 * @brief MinConf::calc_min_conflict_species adds and removes all *absent*
 * species one-by-one to figure out which ones would cause the lowest error.
 * @param site
 * @return all *absent* species that would lead to the lowest error if they were
 * present.
 */
std::vector<unsigned> MinConf::calc_min_conflict_species(const unsigned site) {
  // Get index of all non-present species of this site
  std::vector<unsigned> absent_species_idx = absent_species_index(site);

  // makes sure that the first actual error_ is smaller
  unsigned error = std::numeric_limits<unsigned>::max();
  std::vector<unsigned> min_conflict_species;

  // try for each species at this site where the enery would be minimal
  for (unsigned species_idx = 0; species_idx < absent_species_idx.size();
       species_idx++) {
    const unsigned species = absent_species_idx[species_idx];
    solution[site][species] = 1; // assign species (will be un-done later)
    update_solution_commonness();
    unsigned error_ = calc_error();

    if (error_ < error) {
      min_conflict_species.clear(); // found better fitting species delete other
      min_conflict_species.push_back(species);
      error = error_;
    } else if (error_ == error) {
      // as good as other species, add this species
      min_conflict_species.push_back(species);
    }
    solution[site][species] = 0; // undo assign species
  }
  return min_conflict_species;
}

/**
 * @brief MinConf::gen_init_solution generates a random solution
 */
void MinConf::gen_init_solution() {
  const auto missing_species = calc_missing_species();

  for (unsigned site = 0; site < n_sites; site++) {
    auto absent_species = absent_species_index(site);

    std::shuffle(absent_species.begin(), absent_species.end(), rng);
    for (unsigned s = 0; s < missing_species[site]; s++) {
      const unsigned species = absent_species[s];
      solution[site][species] = 1;
    }
  }
}

/**
 * @brief MinConf::update_solution_commonness
 */
void MinConf::update_solution_commonness() {
  for (unsigned site = 0; site < n_sites - 1; site++) {
    for (unsigned other_site = site + 1; other_site < n_sites; other_site++) {
      commonness[site][other_site] =
          std::inner_product(solution[site].begin(), solution[site].end(),
                             solution[other_site].begin(), 0);
    }
  }
}

/**
 * @brief MinConf::calc_error the error is the Hamming distance between the
 * commonness and the target
 */
unsigned MinConf::calc_error() {
  unsigned sum_diff = 0;
  for (unsigned site = 0; site < n_sites; site++) {
    for (unsigned other_site = site + 1; other_site < n_sites; other_site++) {
      if (target[site][other_site] == NA) {
        continue;
      }
      auto tmp =
          std::abs(commonness[site][other_site] - target[site][other_site]);
      sum_diff += tmp;
    }
  }
  return sum_diff;
}
