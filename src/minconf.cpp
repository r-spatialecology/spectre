#include "minconf.h"
#include <algorithm>
#include <chrono>
#include <iostream>

int MinConf::optimize0(long max_steps_, double max_energy, long long seed)
{
    // Random number generator
    if (seed == 0) {
        seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    }
    this->seed = seed;
    rng = std::mt19937(seed);
    std::uniform_int_distribution<unsigned> site_dist(0, n_sites - 1);
    auto iter = max_steps_;
    //gen_init_solution();
    if (fixed_species.size()) {
        set_fixed_species();
    }

    // optimize
    int patience_counter = 0;
    auto missing_species = alpha_list;
    for (unsigned site = 0; site < n_sites; site++) {
        missing_species[site] -= present_species_index(site, false).size();
    }

    const auto commonness = calculate_commonness();
    auto current_best_solution = solution;
    auto energy = calc_energy(commonness, target);
    double min_energy = energy;
    unsigned n_resets = 1;
    iteration_count.push_back(0);
    energy_vector.push_back(energy);

    while(energy > max_energy) {
        if(iter-- == 0) {
            break;
        }

        if (patience_counter >= 2000) {
            if (min_energy <= energy_vector.back()) {
                n_resets++;
                if (n_resets > n_sites) {
                    n_resets = n_sites;
                }
            } else {
                n_resets = 1;
            }

            for (unsigned i = 0; i < n_resets; i++) {
                const auto site = site_dist(rng);
                for (unsigned species = 0; species < gamma_div; species++) {
                    solution[site][species] = 0;
                }
                if (fixed_species.size()) {
                    set_fixed_species(site);
                }
                missing_species[site] = alpha_list[site] - present_species_index(site, false).size();
            }
            patience_counter = -1000;
        }

        const auto site = site_dist(rng);

        if (!missing_species[site]) {
            // remove a random species at this site
            std::vector<unsigned> species_idx = present_species_index(site, true);
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
            current_best_solution = solution;
            min_energy = energy;
        }

        if (energy >= energy_vector.back()) {
            patience_counter++;
        } else {
            patience_counter = 0;
        }

        iteration_count.push_back(max_steps_ - iter);
        energy_vector.push_back(energy);
    }

    if(add_missing_species(missing_species)) {
        solution_has_best_energy = false;
    }

    energy = calc_energy(calculate_commonness(solution), target);

    if (calc_energy(calculate_commonness(current_best_solution), target) <
            energy) {
        solution = current_best_solution;
    }

    if (iteration_count.back() < (max_steps_ - iter)) {
        iteration_count.push_back(max_steps_ - iter);
        energy_vector.push_back(energy_vector.back());
    }

    return iter;
}

int MinConf::optimize1(long max_steps_, double max_energy, long long seed)
{
    // Random number generator
    if (seed == 0) {
        seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    }
    this->seed = seed;
    rng = std::mt19937(seed);
    std::uniform_int_distribution<unsigned> site_dist(0, n_sites - 1);

    auto iter = max_steps_;
    solution = gen_random_solution();
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
    this->seed = seed;
    rng = std::mt19937(seed);
    std::uniform_int_distribution<unsigned> site_dist(0, n_sites - 1);

    auto iter = max_steps_;
    auto missing_species = alpha_list;
    for (unsigned site = 0; site < n_sites; site++) {
        missing_species[site] -= present_species_index(site).size();
    }

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
            const bool omit_fixed_species = false;
            const auto species_present = present_species_index(site, omit_fixed_species);
            missing_species[site] = alpha_list[site] - species_present.size();
        }
    } else {
        missing_species[0] = 0; // first site is solved anyway
    }

    iteration_count.push_back(0);
    energy_vector.push_back(calc_energy(calculate_commonness(), target));

    // solve dependencies site-by-site
    for (int site = 0; static_cast<unsigned>(site) < (n_sites - 1); site++) {
        for (int other_site = site + 1; static_cast<unsigned>(other_site) < n_sites; other_site++) {
            if (target[site][other_site] <= 0) { // no common species or NA
                continue;
            }
            auto commonness = calculate_commonness();
            int common_species = target[site][other_site] - commonness[site][other_site];

            while (common_species > missing_species[other_site]) {
                if (common_species > 0) {
                    remove_species_max_conf(other_site, target);
                    commonness = calculate_commonness();
                    common_species = target[site][other_site] - commonness[site][other_site];
                    missing_species[other_site]++;
                } else {
                    remove_species_max_conf(site, target);
                    commonness = calculate_commonness();
                    common_species = target[site][other_site] - commonness[site][other_site];
                    missing_species[site]++;
                }
                iter--;
            }

            while (missing_species[other_site] &&
                   common_species &&
                   iter > 0) {
                add_species_min_conf(other_site, target);
                iteration_count.push_back(max_steps_ - iter);
                energy_vector.push_back(calc_energy(calculate_commonness(), target));
                missing_species[other_site]--;
                common_species--;
                iter--;
            }
        }

        // assign all remaining species to the next site
        while (missing_species[site] &&
               iter > 0) {
            add_species_min_conf(site, target);
            iteration_count.push_back(max_steps_ - iter);
            energy_vector.push_back(calc_energy(calculate_commonness(), target));
            missing_species[site]--;
            iter--;
        }
    }

    // if we could't solve above, at least make alpha diversity fit
    for (unsigned site = 0; site < n_sites; site++) {
        while (missing_species[site] > 0) {
            add_species_min_conf(site, target);
            iteration_count.push_back(max_steps_ - iter);
            energy_vector.push_back(calc_energy(calculate_commonness(), target));
            missing_species[site]--;
        }
    }

    return optimize0(iter, max_energy, seed);;
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
        //        if (tabu_species_list.size()) {
        //            tabu_species_list.erase(tabu_species_list.begin());
        //            tabu_species_list.push_back(species);
        //        }
    }
}

bool MinConf::remove_species_max_conf(unsigned site,
                                      const std::vector<std::vector<int> > &target)
{
    // Get index of all present species of this site
    std::vector<unsigned> site_present_species = present_species_index(site);
    if (site_present_species.size()) {
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
        //        if (tabu_species_list.size()) {
        //            // remove tabu_species from present_species_index
        //            for (unsigned idx = 0; idx < tabu_species_list.size(); idx++) {
        //                const unsigned tabu_species = tabu_species_list[idx];
        //                std::vector<unsigned>::iterator it = std::find(site_present_species.begin(),
        //                                                               site_present_species.end(),
        //                                                               tabu_species);
        //                if (it != site_present_species.end()) {
        //                    site_present_species.erase(it);
        //                }
        //            }
        //        }


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
    double energy = std::numeric_limits<double>::max(); // makes sure that the first energy_ is smaller
    std::vector<unsigned> min_conflict_species;

    // try for each species at this site where the enery would be minimal

    for (unsigned species_idx = 0; species_idx < free_species.size(); species_idx++) {
        const unsigned species = free_species[species_idx];
        solution[site][species] = 1; // assign species (will be un-done later)
        const auto commoness_new = calculate_commonness();

        double energy_ = calc_energy(commoness_new, target, true); // may use a different matrix norm

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

        double energy_ = calc_energy(commoness_new, target, true); // may use a different matrix norm

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


//// worst site
//std::vector<unsigned> MinConf::worst_sites(const std::vector<std::vector<int> > &commonness,
//                                           const std::vector<std::vector<int> > &target) {
//    std::vector<unsigned> worst_site(n_sites);
//    std::vector<std::pair<double, unsigned> > energy_site(n_sites);

//    for (unsigned site = 0; site < n_sites; site++) {
//        //        if (tabu_sites_list.size()) { // omit this site if on the tabu_sites_list
//        //            std::vector<unsigned>::iterator it = std::find(tabu_sites_list.begin(),
//        //                                                           tabu_sites_list.end(),
//        //                                                           site);
//        //            if (it != tabu_sites_list.end()) {
//        //                continue;
//        //            }
//        //        }
//        // calculate energy w/o the current site
//        const double energy = calc_energy(commonness, target, true, site); // may use a different matrix norm
//        energy_site[site] = std::make_pair(energy, site);
//    }
//    std::sort(energy_site.begin(), energy_site.end());
//    for (unsigned site = 0; site < n_sites; site++) {
//        worst_site[site] = energy_site[site].second;
//    }

//    return worst_site;
//}

int MinConf::next_site(const std::vector<unsigned> &missing_species)
{
    long long min_values = std::numeric_limits<long long>::max();
    std::vector<unsigned> mrv_site;

    for (unsigned site = 0; site < n_sites; site++) {
        if (!missing_species[site]) {
            continue;
        }
        const auto absent_species = absent_species_index(site).size();
        unsigned n_minus_r = absent_species - missing_species[site];
        long long values = std::numeric_limits<long long>::max();
        if (n_minus_r == 0) {
            values = 1;
        } else if (n_minus_r == 1) {
            values = absent_species;
        } else if (absent_species < factorials.size()) {
            values = factorials[absent_species] /
                    (factorials[missing_species[site]] * factorials[n_minus_r]);
        }
        if (values < min_values) {
            mrv_site.clear();
            mrv_site.push_back(site);
            min_values = values;
        } else if (values == min_values) {
            mrv_site.push_back(site);
        }
    }

    std::shuffle(mrv_site.begin(), mrv_site.end(), rng);
    if (mrv_site.size()) {
        return mrv_site[0];
    }
    return -1;
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
