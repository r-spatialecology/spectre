#ifndef MINCONF_H
#define MINCONF_H
#include "constraint_satisfaction_problem.h"
#include <random>


class MinConf : public Constraint_satisfaction_problem
{
public:
    using Constraint_satisfaction_problem::Constraint_satisfaction_problem;
    int optimize0(const long max_steps_ = 5000, const double max_energy = 0.0, long long seed = 0, const unsigned patience = 2000);
    std::vector<int> iteration_count;
    std::vector<double> energy_vector;
    bool solution_has_best_energy = true;
    double calc_energy_random_solution(const unsigned n = 10);

    long long getSeed() const;
    void setSeed(long long value);

protected:
    std::mt19937 rng;
    long long seed;
    const double epsilon = 0.00001;
    const std::vector<long long> factorials {1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800, 39916800, 479001600};

    std::vector<std::vector<int> > gen_random_solution();
    void set_fixed_species();
    void set_fixed_species(unsigned site);
    void add_species_min_conf(unsigned site,
                              const std::vector<std::vector<int> > &target);
    bool remove_species_max_conf(unsigned site,
                                 const std::vector<std::vector<int> > &target);
    std::vector<unsigned> calc_min_conflict_species(const unsigned site,
                                                    const std::vector<unsigned> free_species,
                                                    const std::vector<std::vector<int> > &target);
    std::vector<unsigned> calc_max_conflict_species(const unsigned site,
                                                    const std::vector<unsigned> pesent_species,
                                                    const std::vector<std::vector<int> > &target);
//    std::vector<unsigned> worst_sites(const std::vector<std::vector<int> > &commonness,
//                                      const std::vector<std::vector<int> > &target);
    int next_site(const std::vector<unsigned> &missing_species);
    bool add_missing_species(std::vector<unsigned> &missing_species);
};

#endif // MINCONF_H
