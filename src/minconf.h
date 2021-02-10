#ifndef MINCONF_H
#define MINCONF_H
#include <random>
#include <string>
#include <vector>

class MinConf {
public:
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
  MinConf(const std::vector<unsigned> &alpha_list, const unsigned gamma_div,
          const std::vector<int> &target,
          const std::vector<int> &partial_solution,
          const std::vector<int> &fixed_species, const unsigned long seed,
          const int na_val = -2147483648);

  /**
   * @brief optimize runs the actual optimization algorithm
   * @param max_steps
   * @param verbose enables a progress bar to the R terminal. This will slow
   * down the optimization.
   * @param interruptible makes the optimization interruptable from R (this is
   * computationally a rather cheap option)
   * @return
   */
  int optimize(const long max_steps, bool verbose, bool interruptible);

  std::vector<std::vector<int>> solution;
  std::vector<std::vector<int>> commonness;
  std::vector<int> iteration_count;
  std::vector<unsigned> error_vector;
  bool solution_has_best_error = true;
  const int RET_ABORT = -999;
  const int NA;

protected:
  std::mt19937 rng;
  const double epsilon = 0.00001;

  std::vector<std::vector<int>> target;
  const std::vector<unsigned> alpha_list;
  std::vector<std::vector<int>> fixed_species;
  const unsigned gamma_div;
  const unsigned n_sites;

  void gen_init_solution(std::vector<unsigned> missing_species);
  std::vector<unsigned> calc_missing_species();
  std::vector<unsigned> present_species_index(unsigned site,
                                              bool omit_fixed_species = true);
  std::vector<unsigned>
  present_species_index(unsigned site,
                        const std::vector<std::vector<int>> partial_solution);
  std::vector<unsigned> absent_species_index(unsigned site);
  bool remove_random_species(const unsigned site);
  void add_species_min_conf(unsigned site,
                            const std::vector<std::vector<int>> &target);
  std::vector<unsigned>
  calc_min_conflict_species(const unsigned site,
                            const std::vector<unsigned> free_species,
                            const std::vector<std::vector<int>> &target);

  unsigned calc_error(const std::vector<std::vector<int>> &commonness,
                      const std::vector<std::vector<int>> &target);
  void update_solution_commonness();
};

#endif // MINCONF_H
