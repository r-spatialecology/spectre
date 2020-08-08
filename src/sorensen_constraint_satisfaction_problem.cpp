#include "sorensen_constraint_satisfaction_problem.h"
#include <iostream>
#include <algorithm>
#include <cmath>

Sorensen_constraint_satisfaction_problem::Sorensen_constraint_satisfaction_problem(const std::vector<unsigned> &alpha_list,
                                                                 const unsigned gamma_div,
                                                                 const std::vector<double> &target_,
                                                                 const std::vector<int> &fixed_species_,
                                                                 const std::vector<int> &partial_solution, const std::string norm)
    : alpha_list(alpha_list), gamma_div(gamma_div), n_sites(alpha_list.size()), norm(norm)
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

// MSP start
std::vector<std::vector<double> > Sorensen_constraint_satisfaction_problem::calculate_sorensen()
{
    return calculate_sorensen(solution);
}

std::vector<std::vector<double> > Sorensen_constraint_satisfaction_problem::calculate_sorensen(const std::vector<std::vector<int> > &solution) {
    
    // create commonness matrix 
    std::vector<std::vector<int> > commonness_result(n_sites, std::vector<int>(n_sites));
    
#pragma omp parallel for
    for (unsigned site = 0; site < n_sites; site++) {
        update_solution_commonness_site(solution, commonness_result, n_sites, gamma_div, site);
    }
    
    
    // create sorensen matrix
    std::vector<std::vector<double> > sorensen_result(n_sites, std::vector<double>(n_sites));
    
#pragma omp parallel for
    for (unsigned site = 0; site < n_sites; site++) {
        update_solution_sorensen_site(solution, commonness_result, sorensen_result, n_sites, gamma_div, site);
    }
    
    return sorensen_result;
}
// MSP end 

// MSP start :  calc energy for sorensen, too 

double Sorensen_constraint_satisfaction_problem::calc_energy_sorensen(const std::vector<std::vector<double> > &sorensen,
                                                             const std::vector<std::vector<double> > &target_sorensen,
                                                             const bool use_custom_norm,
                                                             int omit_site)
{
    double retval = 0;
    if (use_custom_norm && norm != "sum") {
        if (norm == "euclid") {
            retval = calc_energy_euclid_sorensen(sorensen, target_sorensen, omit_site);
        } else if (norm == "max") {
            retval = calc_energy_max_sorensen(sorensen, target_sorensen, omit_site);
        } else if (norm == "p") {
            retval = calc_energy_p_sorensen(sorensen, target_sorensen, omit_site);
        }
    } else {
        retval = calc_energy_sum_sorensen(sorensen, target_sorensen, omit_site);
    }
    return retval;
}

double Sorensen_constraint_satisfaction_problem::calc_energy_sum_sorensen(const std::vector<std::vector<double> > &sorensen,
                                                                 const std::vector<std::vector<double> > &target_sorensen,
                                                                 int omit_site)
{
    double sum_diff = 0.0; // MSP: long double better here??? 
    for (unsigned site = 0; site < n_sites; site++) {
        for (unsigned other_site = 0; other_site < n_sites; other_site++) {
            if (target_sorensen[site][other_site] < 0
                    || static_cast<int>(site) == omit_site
                    || static_cast<int>(other_site) == omit_site) {// i.e. NA
                    continue;
            }
            sum_diff += std::abs(sorensen[site][other_site] -
                target_sorensen[site][other_site]);
        }
    }
    
    return sum_diff;
}

double Sorensen_constraint_satisfaction_problem::calc_energy_euclid_sorensen(const std::vector<std::vector<double> > &sorensen,
                                                                    const std::vector<std::vector<double> > &target_sorensen,
                                                                    int omit_site)
{
    double sum_diff = 0.0;
    for (unsigned site = 0; site < n_sites; site++) {
        for (unsigned other_site = 0; other_site < n_sites; other_site++) {
            if (target_sorensen[site][other_site] < 0
                    || static_cast<int>(site) == omit_site
                    || static_cast<int>(other_site) == omit_site) {// i.e. NA
                    continue;
            }
            const auto diff = sorensen[site][other_site] - target_sorensen[site][other_site];
            sum_diff += diff * diff;
        }
    }
    
    return std::sqrt(sum_diff);
}

double Sorensen_constraint_satisfaction_problem::calc_energy_p_sorensen(const std::vector<std::vector<double> > &sorensen,
                                                               const std::vector<std::vector<double> > &target_sorensen,
                                                               const int omit_site)
{
    double sum_diff = 0;
    for (unsigned site = 0; site < n_sites; site++) {
        for (unsigned other_site = 0; other_site < n_sites; other_site++) {
            if (target_sorensen[site][other_site] < 0
                    || static_cast<int>(site) == omit_site
                    || static_cast<int>(other_site) == omit_site) {// i.e. NA
                    continue;
            }
            const auto diff = std::abs(sorensen[site][other_site] - target_sorensen[site][other_site]);
            sum_diff += std::pow(diff, p);
        }
    }
    
    return std::pow(sum_diff, 1.0 / p);
}

double Sorensen_constraint_satisfaction_problem::calc_energy_max_sorensen(const std::vector<std::vector<double> > &sorensen,
                                                                 const std::vector<std::vector<double> > &target_sorensen,
                                                                 int omit_site)
{
    long max = 0.0;
    for (unsigned site = 0; site < n_sites; site++) {
        for (unsigned other_site = 0; other_site < n_sites; other_site++) {
            if (target_sorensen[site][other_site] < 0
                    || static_cast<int>(site) == omit_site
                    || static_cast<int>(other_site) == omit_site) {// i.e. NA
                    continue;
            }
            const auto diff = std::abs(sorensen[site][other_site] - target_sorensen[site][other_site]);
            if (diff > max) {
                max = diff;
            }
        }
    }
    return max;
}

// MSP end 


std::vector<unsigned> Sorensen_constraint_satisfaction_problem::present_species_index(unsigned site, bool omit_fixed_species)
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

std::vector<unsigned> Sorensen_constraint_satisfaction_problem::present_species_index(unsigned site,
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

std::vector<unsigned> Sorensen_constraint_satisfaction_problem::absent_species_index(unsigned site)
{
    std::vector<unsigned> absent_species_idx;
    for (unsigned species = 0; species < gamma_div; species++) {
        if (!solution[site][species]) {
            absent_species_idx.push_back(species);
        }
    }
    return absent_species_idx;
}


void Sorensen_constraint_satisfaction_problem::update_solution_commonness_site(const std::vector<std::vector<int> > &solution_matrix,
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

// MSP 
void Sorensen_constraint_satisfaction_problem::update_solution_sorensen_site(const std::vector<std::vector<int> > &solution_matrix,
                                                                    std::vector<std::vector<int> > &commonness_result,
                                                                    std::vector<std::vector<double> > &sorensen_result,
                                                                    const unsigned n_sites,
                                                                    const unsigned n_species,
                                                                    const unsigned site)
{
    for (unsigned other_site = 0; other_site < n_sites; other_site++) {
        if (site == other_site) {
            sorensen_result[site][other_site] = -1; // i.e. NA
            continue;
        } else {
            // sorensen richness variables for site a and site b
            int site_a = 0;
            int site_b = 0;
            for (unsigned species = 0; species < n_species; species++) {
                if (solution_matrix[site][species] ) { // present species at current site
                    site_a ++;
                }
                if (solution_matrix[other_site][species]) { // present species at other site
                    site_b ++;
                } 
            }
            // sorensen commonness variable 
            int c_temp = commonness_result[site][other_site]; // tricky! 
            sorensen_result[site][other_site] = 1 - 2.0 * c_temp / (site_a + site_b);
        }
    }
}


// MSP end

