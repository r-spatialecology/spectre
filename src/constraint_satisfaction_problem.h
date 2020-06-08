#ifndef CONSTRAINT_SATISFACTION_PROBLEM_H
#define CONSTRAINT_SATISFACTION_PROBLEM_H
#include <vector>

class Constraint_satisfaction_problem
{
public:
    Constraint_satisfaction_problem(const std::vector<unsigned> &alpha_list,
                                    const unsigned gamma_div,
                                    const std::vector<int> &target,
                                    const std::vector<int> &fixed_species_ = std::vector<int>());

    std::vector<std::vector<int> > solution;

protected:
    long max_steps = 0;
    std::vector<std::vector<int> > fixed_species;
    std::vector<std::vector<int> > fixed_species_idx;
    std::vector<std::vector<int> > target;
    const std::vector<unsigned> alpha_list;
    const unsigned gamma_div;
    const unsigned n_sites;
    double calc_energy(const std::vector<std::vector<int> > &commonness,
                       const std::vector<std::vector<int> > &target);
    std::vector<unsigned> present_species_index(unsigned site);
    std::vector<unsigned> present_species_index(unsigned site, std::vector<std::vector<int> > &mat,
                                                std::vector<std::vector<int> > fixed_species = std::vector<std::vector<int> >());
    std::vector<unsigned> absent_species_index(unsigned site);
    std::vector<std::vector<int> > calculate_commonness();
    void update_solution_commonness_site(const std::vector<std::vector<int> > &solution_matrix,
                                         std::vector<std::vector<int> > &solution_commonness,
                                         const unsigned n_sites,
                                         const unsigned n_species,
                                         const unsigned site);
};

#endif // CONSTRAINT_SATISFACTION_PROBLEM_H
