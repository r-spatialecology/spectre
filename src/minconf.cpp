#include "minconf.h"
#include <algorithm>
#include "calculate_solution_commonness.h"
#include <chrono>


int MinConf::optimize(long max_steps_, double max_energy, long long seed)
{
    // Random number generator
    RNGScope scope; // needed for debugging in Qt Creator
    if (seed == 0) {
        seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    }
    rng = std::mt19937(seed);
    std::uniform_int_distribution<unsigned> site_dist(0, n_sites - 1);

    auto iter = max_steps_;
    auto missing_species = alpha_list;
    std::vector<unsigned> max_common_species(n_sites);
    std::vector<int> currently_added_species(n_sites);

    // find the highest commonness for each site and overall
    unsigned max_common_species_all = 0;
    for (unsigned site = 0; site < n_sites; site++) {
        max_common_species[site] = *(std::max_element(target[site].begin(), target[site].end()));
        if (max_common_species_all < max_common_species[site]) {
            max_common_species_all = max_common_species[site];
        }
    }

    // assign all species that are necessary to fulfill commonness first
    for (unsigned common_species = 1; common_species <= max_common_species_all; common_species++) {
        // set a new target with max. common_species
        const auto target_tmp = target_cur(common_species);

        // add species
        for (unsigned site = 0; site < n_sites; site++) {
            if (max_common_species[site] >= common_species) {
                currently_added_species[site] = add_species_min_conf(site, target_tmp);
                missing_species[site]--;
            }
        }

        // optimize last added species
        auto energy = calc_energy(calculate_commonness(), target_tmp);
        while(energy > max_energy) {
            if(iter-- == 0) {
                break;
            }

            const auto site = site_dist(rng);
            solution[site][currently_added_species[site]] = 0;
            currently_added_species[site] = add_species_min_conf(site, target_tmp);
            const auto commonness = calculate_commonness();
            energy = calc_energy(commonness, target_tmp);
        }
        if (iter < 1) {
            break;
        }
        solved_species++;
    }

    // get the indices of sites where still species are missing
    std::vector<unsigned> sites_idx;
    for (unsigned site = 0; site < n_sites; site++) {
        if(missing_species[site] > 0) {
            sites_idx.push_back(site);
        }    // unsigned max_alpha = *(std::max_element(alpha_list.begin(), alpha_list.end()));

    }

    // assign all remaining species
    while (sites_idx.size()) {
        std::shuffle(sites_idx.begin(), sites_idx.end(), rng);
        const auto site = sites_idx.back();

        if(missing_species[site]--) {
            currently_added_species[site] = add_species_min_conf(site, target);
        } else {
            sites_idx.pop_back();
        }
    }

    return iter;
}

std::vector<std::vector<int> > MinConf::target_cur(int n_common)
{
    auto result = target;
    for (unsigned site = 0; site < n_sites; site++) { // could be optimized due to symmetry
        for (unsigned other_site = 0; other_site < n_sites; other_site++) {
            if (result[site][other_site] > n_common) {
                result[site][other_site] = n_common;
            }
        }
    }
    return result;
}

int MinConf::add_species_min_conf(unsigned site, const std::vector<std::vector<int> > &target)
{
    // Get index of all non-present species of this site
    std::vector<unsigned> free_species_idx;
    for (unsigned species = 0; species < gamma_div; species++) {
        if (!solution[site][species]) {
            free_species_idx.push_back(species);
        }
    }

    // calculate the best-fitting species for this site
    auto min_conflict_species = calc_min_conflict_species(site,
                                                          free_species_idx,
                                                          target);

    if (min_conflict_species.size() < 1) {
        return -1; // something odd happened
    }

    std::shuffle(min_conflict_species.begin(), min_conflict_species.end(), rng);
    const auto species = min_conflict_species[0];
    solution[site][species] = 1;

    return species;
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

std::vector<std::vector<int> > MinConf::calculate_commonness()
{
    std::vector<std::vector<int> > result(n_sites, std::vector<int>(n_sites));

#pragma omp parallel for
    for (unsigned site = 0; site < n_sites; site++) {
        update_solution_commonness_site(solution, result, n_sites, gamma_div, site);
    }

    return result;
}

double MinConf::calc_energy(const std::vector<std::vector<int> > &commonness,
                            const std::vector<std::vector<int> > &target)
{
    long long sum_target = 0;
    long long sum_diff = 0;
    for (unsigned site = 0; site < n_sites; site++) {
        for (unsigned other_site = 0; other_site < n_sites; other_site++) {
            sum_diff += std::abs(commonness[site][other_site] -
                                 target[site][other_site]);
            sum_target += target[site][other_site];
        }
    }

    return static_cast<double>(sum_diff) / static_cast<double>(sum_target);
}
