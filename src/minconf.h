#ifndef MINCONF_H
#define MINCONF_H
#include "constraint_satisfaction_problem.h"
#include <random>


class MinConf : public Constraint_satisfaction_problem
{
public:
    using Constraint_satisfaction_problem::Constraint_satisfaction_problem;
    int optimize1(long max_steps_ = 5000, double max_energy = 0.0, long long seed = 0);
    int optimize2(long max_steps_ = 5000, double max_energy = 0.0, long long seed = 0);
    std::vector<int> iteration_count;
    std::vector<double> energy_vector;

protected:
    std::mt19937 rng;

    void gen_init_solution();
    void set_fixed_species();
    void add_species_min_conf(unsigned site, const std::vector<std::vector<int> > &target);
    bool remove_species_max_conf(unsigned site, const std::vector<std::vector<int> > &target);
    std::vector<unsigned> calc_min_conflict_species(const unsigned site,
                                                    const std::vector<unsigned> free_species,
                                                    const std::vector<std::vector<int> > &target);
    std::vector<unsigned> calc_max_conflict_species(const unsigned site,
                                                    const std::vector<unsigned> pesent_species,
                                                    const std::vector<std::vector<int> > &target);

};

#endif // MINCONF_H
