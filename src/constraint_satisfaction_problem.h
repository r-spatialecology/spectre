#ifndef CONSTRAINT_SATISFACTION_PROBLEM_H
#define CONSTRAINT_SATISFACTION_PROBLEM_H
#include <vector>

class Constraint_satisfaction_problem
{
public:
    Constraint_satisfaction_problem(const std::vector<unsigned> &alpha_list,
                                    const unsigned gamma_div,
                                    const std::vector<int> &target,
                                    const std::vector<int> &fixed_species = std::vector<int>());

    std::vector<std::vector<int> > solution;
    std::vector<unsigned> solved_sites;
    int solved_species = 0;

protected:
    long max_steps = 0;

    std::vector<std::vector<int> > target;
    const std::vector<unsigned> alpha_list;
    const unsigned gamma_div;
    const unsigned n_sites;
};

#endif // CONSTRAINT_SATISFACTION_PROBLEM_H
