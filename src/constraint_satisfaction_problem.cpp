#include "constraint_satisfaction_problem.h"
#include <iostream>
#include <algorithm>

Constraint_satisfaction_problem::Constraint_satisfaction_problem(const std::vector<unsigned> &alpha_list,
                                                                 const unsigned gamma_div,
                                                                 const std::vector<int> &target_,
                                                                 const std::vector<int> &fixed_species_,
                                                                 const std::vector<int> &partial_solution)
    : alpha_list(alpha_list), gamma_div(gamma_div), n_sites(alpha_list.size())
{
    solution.resize(n_sites);
    target.resize(n_sites);
    //  tabu_list.resize(tabu, std::numeric_limits<unsigned>::max()); // species list should be smaller than that...

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

    if (partial_solution.size() > 1) {
        if(partial_solution.size() != n_sites * gamma_div) {
            std::cerr << "The size of the partial_solution vector does not match n_sites * gamma_div. "
                      << "partial_solution ignored."
                      << std::endl;
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

    if (fixed_species_.size() > 1) {
        if(fixed_species_.size() != n_sites * gamma_div) {
            std::cerr << "The size of the fixed_species vector does not match n_sites * gamma_div. "
                      << "fixed_species ignored."
                      << std::endl;
        } else {
            fixed_species.resize(n_sites);
            for (unsigned site = 0; site < n_sites; site++) {
                fixed_species[site].resize(gamma_div);
                for (unsigned species = 0; species < gamma_div; species++) {
                    if(fixed_species_[site * gamma_div + species]) {
                        fixed_species[site][species] = 1;
                    }
                }
            }
        }
    }
}

std::vector<std::vector<int> > Constraint_satisfaction_problem::calculate_commonness()
{
    return calculate_commonness(solution);
}

std::vector<std::vector<int> > Constraint_satisfaction_problem::calculate_commonness(const std::vector<std::vector<int> > &solution) {
    std::vector<std::vector<int> > result(n_sites, std::vector<int>(n_sites));

#pragma omp parallel for
    for (unsigned site = 0; site < n_sites; site++) {
        update_solution_commonness_site(solution, result, n_sites, gamma_div, site);
    }

    return result;
}

double Constraint_satisfaction_problem::calc_energy(const std::vector<std::vector<int> > &commonness,
                                                    const std::vector<std::vector<int> > &target,
                                                    const std::string severity,
                                                    int omit_site)
{
    long long sum_target = 0;
    long long sum_diff = 0;
    for (unsigned site = 0; site < n_sites; site++) {
        for (unsigned other_site = 0; other_site < n_sites; other_site++) {
            if (target[site][other_site] < 0
                    || static_cast<int>(site) == omit_site
                    || static_cast<int>(other_site) == omit_site) {// i.e. NA
                continue;
            }
            if (severity == "commonness") {
                // be strict: there can't be more species in common than allowed in the target matrix
                // (less is ok)
                const long long diff = commonness[site][other_site] -
                        target[site][other_site];
                if (diff > 0) {
                    return std::numeric_limits<double>::max();
                }
                sum_diff += std::abs(diff);
            }
            sum_diff += std::abs(commonness[site][other_site] -
                                 target[site][other_site]);
            sum_target += target[site][other_site];
        }
    }

    return static_cast<double>(sum_diff) / static_cast<double>(sum_target);
}

std::vector<unsigned> Constraint_satisfaction_problem::present_species_index(unsigned site, bool omit_fixed_species)
{
    std::vector<unsigned> species_idx;
    for (unsigned species = 0; species < gamma_div; species++) {
        if (omit_fixed_species && fixed_species.size()) {
            if (fixed_species[site][species] == solution[site][species]) {
                continue;
            }
        }
        if (solution[site][species]) {
            species_idx.push_back(species);
        }
    }

    return  species_idx;
}

std::vector<unsigned> Constraint_satisfaction_problem::present_species_index(unsigned site,
                                                                             const std::vector<std::vector<int> > partial_solution)
{
    std::vector<unsigned> species_idx;
    for (unsigned species = 0; species < gamma_div; species++) {
        if (partial_solution[site][species]) {
            species_idx.push_back(species);
        }
    }

    return  species_idx;
}

std::vector<unsigned> Constraint_satisfaction_problem::absent_species_index(unsigned site)
{
    std::vector<unsigned> absent_species_idx;
    for (unsigned species = 0; species < gamma_div; species++) {
        if (!solution[site][species]) {
            absent_species_idx.push_back(species);
        }
    }
    return absent_species_idx;
}

void Constraint_satisfaction_problem::Constraint_satisfaction_problem::update_solution_commonness_site(const std::vector<std::vector<int> > &solution_matrix,
                                                                                                       std::vector<std::vector<int> > &solution_commonness,
                                                                                                       const unsigned n_sites,
                                                                                                       const unsigned n_species,
                                                                                                       const unsigned site)
{
    for (unsigned other_site = 0; other_site < n_sites; other_site++) {
        if (site == other_site) {
            solution_commonness[site][other_site] = -1; // i.e. NA
            continue;
        } else {
            for (unsigned species = 0; species < n_species; species++) {
                if (!solution_matrix[site][species]) { // no species at current site
                    continue;
                } else if (!solution_matrix[other_site][species]) { // no species at other site
                    continue;
                } else {
                    solution_commonness[site][other_site]++;
                }
            }
        }
    }
}






