#include "minconf.h"
#include <algorithm>
#include <chrono>
#include <iostream>
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>

MinConf::MinConf(const std::vector<unsigned> &alpha_list,
                 const unsigned gamma_div,
                 const std::vector<int> &target_,
                 const std::vector<int> &fixed_species_,
                 const std::vector<int> &partial_solution)
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

int MinConf::optimize(const long max_steps_, const double max_energy, long long seed, bool verbose, bool interruptible)
{
    Progress p(max_steps_, verbose);

    // Random number generator
    if (seed == 0) {
        seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    }

    this->seed = seed;
    rng = std::mt19937(seed);
    std::uniform_int_distribution<unsigned> site_dist(0, n_sites - 1);
    auto iter = max_steps_;

    if (fixed_species.size()) {
        set_fixed_species();
    }

    // optimize
    auto missing_species = calc_missing_species();

    const auto commonness = calculate_commonness();
    auto current_best_solution = solution;
    double energy = calc_energy(commonness, target);
    double min_energy = std::numeric_limits<double>::max();
    iteration_count.push_back(0);
    energy_vector.push_back(energy);

    while(energy > max_energy) {
        p.increment(); // update progress
        if(interruptible) {
            if (Progress::check_abort() ) { return RET_ABORT; }
        }

        if(iter-- == 0) {
            break;
        }

        const auto site = site_dist(rng); // choose random site

        if (!missing_species[site]) {
            // remove a random species at this site
            std::vector<unsigned> species_idx = present_species_index(site, true);
            if (species_idx.size() == 0) {
                continue; // all species fixed, nothing to remove
            }
            std::shuffle(species_idx.begin(), species_idx.end(), rng);
            const unsigned species = species_idx.back();
            solution[site][species] = 0;
        } else {
            missing_species[site]--;
        }

        // add min conf species
        add_species_min_conf(site, target);
        const auto commonness = calculate_commonness();
        energy = calc_energy(commonness, target);

        if (min_energy > energy) {
            if (correct_alpha_div(missing_species)) {
                current_best_solution = solution;
                min_energy = energy;
            }
        }

        iteration_count.push_back(max_steps_ - iter);
        energy_vector.push_back(energy);
    }

    // add species if there are some still missing after the iterations
    add_missing_species(missing_species);
    energy = calc_energy(calculate_commonness(solution), target);
    if (min_energy > energy) {
        current_best_solution = solution;
        min_energy = energy;
    }

    iteration_count.push_back(iteration_count.back() + 1);
    energy_vector.push_back(energy);
    solution = current_best_solution;

    return iter;
}

double MinConf::calc_energy_random_solution(const unsigned n)
{
    double avg_energy = 0;
    for (unsigned i = 0; i < n; i++) {
        const auto random_solution = gen_random_solution();
        const auto commonness = calculate_commonness(random_solution);
        avg_energy += calc_energy(commonness, target); // MSP
    }

    return avg_energy / n;
}

long long MinConf::getSeed() const
{
    return seed;
}

void MinConf::setSeed(long long value)
{
    seed = value;
    rng.seed(seed);
}

std::vector<std::vector<int> > MinConf::gen_random_solution()
{
    const auto n_sites = alpha_list.size(); // number of cols
    std::vector<std::vector<int> > random_solution(n_sites, std::vector<int>(gamma_div));
    for (unsigned site = 0; site < n_sites; site++) {
        const unsigned alpha = static_cast<unsigned>(alpha_list[site]);
        for (unsigned species = 0; species < alpha; species++) {
            random_solution[site][species] = 1;
        }
        std::shuffle(random_solution[site].begin(), random_solution[site].end(), rng);
    }

    return random_solution;
}

void MinConf::set_fixed_species()
{
    for (unsigned site = 0; site < n_sites; site++) {
        set_fixed_species(site);
    }
}

void MinConf::set_fixed_species(unsigned site)
{
    auto site_present_species = present_species_index(site);
    const auto fixed_species_idx = present_species_index(site, fixed_species);

    // remove fixed species from present_species to avoid removal of that species
    for (unsigned idx = 0; idx < fixed_species_idx.size(); idx++) {
        const unsigned fixed_species = fixed_species_idx[idx];
        std::vector<unsigned>::iterator it = std::find(site_present_species.begin(),
                                                       site_present_species.end(),
                                                       fixed_species);
        if (it != site_present_species.end()) {
            site_present_species.erase(it);
        }
    }

    // Add remaing fixed species
    if (site_present_species.size()) {
        std::shuffle(site_present_species.begin(),
                     site_present_species.end(), rng);
    }

    for (unsigned idx = 0; idx < fixed_species_idx.size(); idx++) {
        const unsigned fixed_species = fixed_species_idx[idx];
        solution[site][fixed_species] = 1;
        if (site_present_species.size()) {
            solution[site][site_present_species.back()] = 0;
            site_present_species.pop_back();
        }
    }
}

std::vector<unsigned> MinConf::calc_missing_species()
{
    auto missing_species = alpha_list;
    for (unsigned site = 0; site < n_sites; site++) {
        const auto present_species = present_species_index(site, false).size();
        missing_species[site] = (missing_species[site] < present_species) ? 0 : missing_species[site] - present_species;
    }

    return missing_species;
}

std::vector<unsigned> MinConf::present_species_index(unsigned site, bool omit_fixed_species)
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

std::vector<unsigned> MinConf::present_species_index(unsigned site, const std::vector<std::vector<int> > partial_solution)
{
    std::vector<unsigned> species_idx;
    for (unsigned species = 0; species < gamma_div; species++) {
        if (partial_solution[site][species]) {
            species_idx.push_back(species);
        }
    }

    return  species_idx;
}

std::vector<unsigned> MinConf::absent_species_index(unsigned site)
{
    std::vector<unsigned> absent_species_idx;
    for (unsigned species = 0; species < gamma_div; species++) {
        if (!solution[site][species]) {
            absent_species_idx.push_back(species);
        }
    }
    return absent_species_idx;
}

void MinConf::add_species_min_conf(unsigned site,
                                   const std::vector<std::vector<int> > &target)
{
    // Get index of all non-present species of this site
    std::vector<unsigned> absent_species_idx = absent_species_index(site);

    // calculate the best-fitting species for this site
    auto min_conflict_species = calc_min_conflict_species(site,
                                                          absent_species_idx,
                                                          target);

    if (min_conflict_species.size() < 1) {
        std::cerr << "no species found to add at add_species_min_conf, site: "
                  << site
                  << std::endl; // something odd happened
    } else {
        std::shuffle(min_conflict_species.begin(), min_conflict_species.end(), rng);
        const auto species = min_conflict_species[0];
        solution[site][species] = 1;
    }
}

std::vector<unsigned> MinConf::calc_min_conflict_species(const unsigned site,
                                                         const std::vector<unsigned> free_species,
                                                         const std::vector<std::vector<int> > &target)
{
    double energy = std::numeric_limits<double>::max(); // makes sure that the first energy_ is smaller
    std::vector<unsigned> min_conflict_species;

    // try for each species at this site where the enery would be minimal

    for (unsigned species_idx = 0; species_idx < free_species.size(); species_idx++) {
        const unsigned species = free_species[species_idx];
        solution[site][species] = 1; // assign species (will be un-done later)
        const auto commoness_new = calculate_commonness();

        double energy_ = calc_energy(commoness_new, target);

        if (energy_ < energy) {
            min_conflict_species.clear(); // found better fitting species delete other
            min_conflict_species.push_back(species);
            energy = energy_;
        } else if (fabs(energy_ - energy) <= epsilon) {
            min_conflict_species.push_back(species); // as good as others, add this species
        }
        solution[site][species] = 0; // undo assign species
    }

    return min_conflict_species;
}

bool MinConf::add_missing_species(std::vector<unsigned> &missing_species)
{
    bool retval = false;
    for (unsigned site = 0; site < n_sites; site++) {
        while (missing_species[site]) {
            add_species_min_conf(site, target);
            missing_species[site]--;
            retval = true;
        }
    }
    return retval;
}

std::vector<std::vector<int> > MinConf::calculate_commonness()
{
    return calculate_commonness(solution);
}

std::vector<std::vector<int> > MinConf::calculate_commonness(const std::vector<std::vector<int> > &solution) {
    std::vector<std::vector<int> > result(n_sites, std::vector<int>(n_sites));

    // #pragma omp parallel for
    for (unsigned site = 0; site < n_sites; site++) {
        update_solution_commonness_site(solution, result, n_sites, gamma_div, site);
    }

    return result;
}

void MinConf::update_solution_commonness_site(const std::vector<std::vector<int> > &solution_matrix, std::vector<std::vector<int> > &solution_commonness, const unsigned n_sites, const unsigned n_species, const unsigned site)
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

double MinConf::calc_energy(const std::vector<std::vector<int> > &commonness,
                            const std::vector<std::vector<int> > &target)
{
    long long sum_diff = 0;
    for (unsigned site = 0; site < n_sites; site++) {
        for (unsigned other_site = 0; other_site < n_sites; other_site++) {
            if (target[site][other_site] < 0) {// i.e. NA
                continue;
            }
            sum_diff += std::abs(commonness[site][other_site] -
                                 target[site][other_site]);
        }
    }

    return static_cast<double>(sum_diff);
}

bool MinConf::correct_alpha_div(const std::vector<unsigned> &missing_species)
{
    for (unsigned site = 0; site < n_sites; site++) {
        if (missing_species[site]) {
            return false;
        }
    }
    return true;
}
