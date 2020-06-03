#include "constraint_satisfaction_problem.h"

Constraint_satisfaction_problem::Constraint_satisfaction_problem(const std::vector<unsigned> &alpha_list,
                                                                 const unsigned gamma_div,
                                                                 const std::vector<int> &target_,
                                                                 const std::vector<int> &fixed_species)
    : alpha_list(alpha_list), gamma_div(gamma_div), n_sites(alpha_list.size())
{
    solution.resize(n_sites);
    target.resize(n_sites);

    for (unsigned site = 0; site < n_sites; site++) {
        solution[site].resize(gamma_div);

        // convert target matrix to a more convenient format
        target[site].resize(n_sites);
        for (unsigned other_site = 0; other_site < n_sites; other_site++) {
            if (site == other_site) {
                target[site][other_site] = -1; // i.e. NA
            } else if (target_[other_site * n_sites + site] < 0) { // i.e. NA
                target[site][other_site] = target_[site * n_sites + other_site];
            } else {
                target[site][other_site] = target_[other_site * n_sites + site];
            }
        }
    }
}
