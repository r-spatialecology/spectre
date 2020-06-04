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

    // random first site
    {
        const unsigned site = 0;
        for (unsigned species = 0; species < alpha_list[site]; species++) {
            solution[site][species] = 1;
            missing_species[site]--;
        }
        std::shuffle(solution[site].begin(), solution[site].end(), rng);
    }

    // solve all dependencies site-by-site
    for (unsigned site = 0; site < (n_sites - 1); site++) {
        for (unsigned other_site = 0; other_site < n_sites; other_site++) {
            int common_species = target[site][other_site];
            if (common_species > 0) {
                if (missing_species[other_site]) {
                    while (common_species--) {
                        currently_added_species[other_site] = add_species_min_conf(other_site, target);
                        missing_species[other_site]--;
                    }
                }
            }
        }

        // assign all remaining species to the next site
        const unsigned next_site = site + 1;
        while (missing_species[next_site]) {
            currently_added_species[next_site] = add_species_min_conf(site + 1, target);
            missing_species[next_site]--;
        }
    }

    // optimize
    auto energy = calc_energy(calculate_commonness(), target);
    while(energy > max_energy) {
        if(iter-- == 0) {
            break;
        }

        const auto site = site_dist(rng);
        solution[site][currently_added_species[site]] = 0; // would be better to use species with least common species
        currently_added_species[site] = add_species_min_conf(site, target);
        const auto commonness = calculate_commonness();
        energy = calc_energy(commonness, target);
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
