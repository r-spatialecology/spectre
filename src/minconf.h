#ifndef MINCONF_H
#define MINCONF_H
#include "constraint_satisfaction_problem.h"
#include <random>


class MinConf : public Constraint_satisfaction_problem
{
public:
    using Constraint_satisfaction_problem::Constraint_satisfaction_problem;
    int optimize(long max_steps_ = 5000, double max_energy = 0.0, long long seed = 0);

private:
    std::mt19937 rng;

    std::vector<std::vector<int> > target_cur(int n_common);
    int add_species_min_conf(unsigned site, const std::vector<std::vector<int> > &target);
    std::vector<unsigned> calc_min_conflict_species(const unsigned site,
                                                    const std::vector<unsigned> free_species,
                                                    const std::vector<std::vector<int> > &target);
    std::vector<std::vector<int> > calculate_commonness();
    double calc_energy(const std::vector<std::vector<int> > &commonness,
                       const std::vector<std::vector<int> > &target);
};

#endif // MINCONF_H
