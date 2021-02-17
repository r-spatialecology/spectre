#ifndef MINCONF_H
#define MINCONF_H
#include <random>
#include <vector>

class MinConf {
public:
  MinConf(const std::vector<unsigned> &alpha_list, const unsigned gamma_div,
          const std::vector<int> &target,
          const std::vector<int> &partial_solution,
          const std::vector<int> &fixed_species, const unsigned long seed,
          const int na_val = -2147483648);

  int optimize(const long max_steps, bool verbose, bool interruptible);

  std::vector<std::vector<int>> solution;
  std::vector<std::vector<int>> commonness;
  std::vector<int> iteration_count;
  std::vector<unsigned> error_vector;
  const int RET_ABORT = -999;
  const int NA;

protected:
  std::mt19937 rng;

  std::vector<std::vector<int>> target;
  const std::vector<unsigned> alpha_list;
  std::vector<std::vector<int>> fixed_species;
  const unsigned gamma_div;
  const unsigned n_sites;

  std::vector<unsigned> calc_missing_species();
  std::vector<unsigned> present_species_index(unsigned site,
                                              bool omit_fixed_species = true);
  std::vector<unsigned> absent_species_index(unsigned site);
  bool remove_random_species(const unsigned site);
  void add_species_min_conf(unsigned site);
  std::vector<unsigned> calc_min_conflict_species(const unsigned site);
  void gen_init_solution();
  void update_solution_commonness();
  unsigned calc_error();
};

#endif // MINCONF_H
