#include "minconf.h"
#include <algorithm>
#include <chrono>
#include <iostream>


int MinConf::optimize1(long max_steps_, double max_energy, long long seed)
{
    // Random number generator
    if (seed == 0) {
        seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    }
    rng = std::mt19937(seed);
    std::uniform_int_distribution<unsigned> site_dist(0, n_sites - 1);

    auto iter = max_steps_;
    gen_init_solution();
    if (fixed_species.size()) {
        set_fixed_species();
    }

    // optimize
    auto energy = calc_energy(calculate_commonness(), target);
    iteration_count.push_back(0);
    energy_vector.push_back(energy);

    while(energy > max_energy) {
        if(iter-- == 0) {
            break;
        }

        const auto site = site_dist(rng);
        // remove a random species at this site
        //        std::vector<unsigned> species_idx = present_species_index(site);
        //        std::shuffle(species_idx.begin(), species_idx.end(), rng);
        //        const unsigned species = species_idx.back();
        //        solution[site][species] = 0;
        // new: remove max_conf species:
        if (remove_species_max_conf(site, target)) {
            // add min conf species
            add_species_min_conf(site, target);
            const auto commonness = calculate_commonness();
            energy = calc_energy(commonness, target);
        }
        iteration_count.push_back(max_steps_ - iter);
        energy_vector.push_back(energy);
    }

    return iter;
}

int MinConf::optimize2(long max_steps_, double max_energy, long long seed)
{
    // Random number generator
    if (seed == 0) {
        seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    }
    rng = std::mt19937(seed);
    std::uniform_int_distribution<unsigned> site_dist(0, n_sites - 1);

    auto iter = max_steps_;
    auto missing_species = alpha_list;

    // random first site
    {
        const unsigned site = 0;
        for (unsigned species = 0; species < alpha_list[site]; species++) {
            solution[site][species] = 1;
        }
        std::shuffle(solution[site].begin(), solution[site].end(), rng);
    }

    if (fixed_species.size()) {
        set_fixed_species();
        for (unsigned site = 0; site < n_sites; site++) {
            const auto species_present = present_species_index(site);
            missing_species[site] = alpha_list[site] - species_present.size();
        }
    }

    // solve all dependencies site-by-site

    for (unsigned site = 0; site < (n_sites - 1); site++) {
        for (unsigned other_site = 0; other_site < n_sites; other_site++) {
            const auto commonness = calculate_commonness();
            int common_species = target[site][other_site] - commonness[site][other_site];
            if (common_species > 0) {
                while (missing_species[other_site] &&
                       common_species--) {
                    add_species_min_conf(other_site, target);
                    missing_species[other_site]--;
                }
            }
        }

        // assign all remaining species to the next site
        const unsigned next_site = site + 1;
        while (missing_species[next_site]) {
            add_species_min_conf(next_site, target);
            missing_species[next_site]--;
        }
    }

    iteration_count.push_back(0);
    energy_vector.push_back(calc_energy(calculate_commonness(), target));

    // optimize
    auto energy = calc_energy(calculate_commonness(), target);
    while(energy > max_energy) {
        if(iter-- == 0) {
            break;
        }

        const auto site = site_dist(rng);
        if (remove_species_max_conf(site, target)) {
            add_species_min_conf(site, target);
            const auto commonness = calculate_commonness();
            energy = calc_energy(commonness, target);
        }
        iteration_count.push_back(max_steps_ - iter);
        energy_vector.push_back(energy);
    }

    return iter;
}

void MinConf::gen_init_solution()
{
    const auto n_sites = alpha_list.size(); // number of cols

    for (unsigned site = 0; site < n_sites; site++) {
        const unsigned alpha = static_cast<unsigned>(alpha_list[site]);
        for (unsigned species = 0; species < alpha; species++) {
            solution[site][species] = 1;
        }
        std::random_shuffle(solution[site].begin(), solution[site].end()); ///BUG: seed is ignored, here
    }
}

void MinConf::set_fixed_species()
{
    for (unsigned site = 0; site < n_sites; site++) {
        if (fixed_species_idx[site].size()) {
            auto site_present_species = present_species_index(site);

            // remove fixed species from present_species to avoid removal of that species
            for (unsigned idx = 0; idx < fixed_species_idx[site].size(); idx++) {
                const unsigned fixed_species = fixed_species_idx[site][idx];
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
            for (unsigned idx = 0; idx < fixed_species_idx[site].size(); idx++) {
                const unsigned fixed_species = fixed_species_idx[site][idx];
                solution[site][fixed_species] = 1;
                if (site_present_species.size()) {
                    solution[site][site_present_species.back()] = 0;
                    site_present_species.pop_back();
                }
            }
        }
    }
}

void MinConf::add_species_min_conf(unsigned site, const std::vector<std::vector<int> > &target)
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

bool MinConf::remove_species_max_conf(unsigned site, const std::vector<std::vector<int> > &target)
{
    // Get index of all present species of this site
    std::vector<unsigned> site_present_species = present_species_index(site);

    if (fixed_species_idx.size()) {
        // remove fixed_species from present_species_index
        for (unsigned idx = 0; idx < fixed_species_idx[site].size(); idx++) {
            const unsigned fixed_species = fixed_species_idx[site][idx];
            std::vector<unsigned>::iterator it = std::find(site_present_species.begin(),
                                                           site_present_species.end(),
                                                           fixed_species);
            if (it != site_present_species.end()) {
                site_present_species.erase(it);
            }
        }
    }

    if (site_present_species.size()) {
        // calculate the best-fitting species for this site
        auto max_conflict_species = calc_max_conflict_species(site,
                                                              site_present_species,
                                                              target);

        if (max_conflict_species.size() < 1) {
            std::cerr << "no species found to remove at remove_species_max_conf, site: "
                      << site
                      << std::endl; // something odd happened
        } else {
            std::shuffle(max_conflict_species.begin(), max_conflict_species.end(), rng);
            const auto species = max_conflict_species[0];
            solution[site][species] = 0;
        }
        return true;
    }
    return false;
}

std::vector<unsigned> MinConf::calc_min_conflict_species(const unsigned site,
                                                         const std::vector<unsigned> free_species,
                                                         const std::vector<std::vector<int> > &target)
{
    const double epsilon = 0.00001;
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

std::vector<unsigned> MinConf::calc_max_conflict_species(const unsigned site,
                                                         const std::vector<unsigned> pesent_species,
                                                         const std::vector<std::vector<int> > &target)
{
    const double epsilon = 0.00001;
    double energy = std::numeric_limits<double>::max(); // makes sure that the first energy_ is smaller
    std::vector<unsigned> max_conflict_species;

    // try for each present species at this site where the energy would be minimal if it was absent

    for (unsigned species_idx = 0; species_idx < pesent_species.size(); species_idx++) {
        const unsigned species = pesent_species[species_idx];
        solution[site][species] = 0; // assign species (will be un-done later)
        const auto commoness_new = calculate_commonness();

        double energy_ = calc_energy(commoness_new, target);

        if (energy_ < energy) {
            max_conflict_species.clear(); // found better fitting species delete other
            max_conflict_species.push_back(species);
            energy = energy_;
        } else if (fabs(energy_ - energy) <= epsilon) {
            max_conflict_species.push_back(species); // as good as others, add this species
        }
        solution[site][species] = 1; // undo assign species
    }

    return max_conflict_species;
}


