#include "sorensen_minconf.h"
#include <algorithm>
#include <cmath>
#include <chrono>
#include <iostream>

// MSP start
int Sorensen_MinConf::optimize0_sorensen(const long max_steps_, const double max_energy, long long seed, const unsigned patience)
{
    // Random number generator
    if (seed == 0) {
        seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    }
    this->seed = seed;
    rng = std::mt19937(seed);
    std::uniform_int_distribution<unsigned> site_dist(0, n_sites - 1);
    auto iter = max_steps_;
    
    if (fixed_species.size()) {
        set_fixed_species_sorensen();
    }
    
    // optimize
    int patience_counter = 0;
    auto missing_species = alpha_list;
    for (unsigned site = 0; site < n_sites; site++) {
        missing_species[site] -= present_species_index(site, false).size();
    }
    
    const auto sorensen = calculate_sorensen();
    auto current_best_solution = solution;
    auto energy = calc_energy_sorensen(sorensen, target);
    double min_energy = energy;
    unsigned n_resets = 1;
    iteration_count.push_back(0);
    energy_vector.push_back(energy);
    
    while(energy > max_energy) {
        if(iter-- == 0) {
            break;
        }
        
        if (patience_counter >= patience) {
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
                    set_fixed_species_sorensen(site);
                }
                missing_species[site] = alpha_list[site] - present_species_index(site, false).size();
            }
            patience_counter = 0;
        }
        
        const auto site = site_dist(rng);
        
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
        add_species_min_conf_sorensen(site, target);
        const auto sorensen = calculate_sorensen();
        energy = calc_energy_sorensen(sorensen, target);
        
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
    
    // add species if there are some still missing
    add_missing_species_sorensen(missing_species);
    energy = calc_energy_sorensen(calculate_sorensen(solution), target);
    iteration_count.push_back(iteration_count.back() + 1);
    energy_vector.push_back(energy);
    
    if (calc_energy_sorensen(calculate_sorensen(current_best_solution), target) <
        energy) {
        solution = current_best_solution;
    }
    
    missing_species = alpha_list;
    for (unsigned site = 0; site < n_sites; site++) {
        missing_species[site] -= present_species_index(site, false).size();
        if (missing_species[site]) {
            solution_has_best_energy = false;
            break;
        }
    }
    
    return iter;
}



// MSP end 


// MSP start

double Sorensen_MinConf::calc_energy_random_solution_sorensen(const unsigned n)
{
    double avg_energy = 0;
    for (unsigned i = 0; i < n; i++) {
        const auto random_solution = gen_random_solution_sorensen();
        const auto sorensen = calculate_sorensen(random_solution);
        avg_energy += calc_energy_sorensen(sorensen, target); // MSP
    }
    
    return avg_energy / n;
}


// MSP end



// MSP start
std::vector<unsigned> Sorensen_MinConf::calc_min_conflict_species_sorensen(const unsigned site,
                                                                           const std::vector<unsigned> free_species,
                                                                           const std::vector<std::vector<double> > &target)
{
    double energy = std::numeric_limits<double>::max(); // makes sure that the first energy_ is smaller
    std::vector<unsigned> min_conflict_species;
    
    // try for each species at this site where the enery would be minimal
    
    for (unsigned species_idx = 0; species_idx < free_species.size(); species_idx++) {
        const unsigned species = free_species[species_idx];
        solution[site][species] = 1; // assign species (will be un-done later)
        
       
        const auto sorensen_new = calculate_sorensen(solution);
        
        double energy_ = calc_energy_sorensen(sorensen_new, target, true); // may use a different matrix norm
        
       
        
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

// MSP end

// MSP

bool Sorensen_MinConf::add_missing_species_sorensen(std::vector<unsigned> &missing_species)
{
    bool retval = false;
    for (unsigned site = 0; site < n_sites; site++) {
        while (missing_species[site]) {
            add_species_min_conf_sorensen(site, target);
            missing_species[site]--;
            retval = true;
        }
    }
    return retval;
}

// MSP end 



// MSP
void Sorensen_MinConf::add_species_min_conf_sorensen(unsigned site,
                                                     const std::vector<std::vector<double> > &target)
{
    // Get index of all non-present species of this site
    std::vector<unsigned> absent_species_idx = absent_species_index(site);
    
    // calculate the best-fitting species for this site
    auto min_conflict_species = calc_min_conflict_species_sorensen(site,
                                                                   absent_species_idx,
                                                                   target);
    
    if (min_conflict_species.size() < 1) {
        std::cerr << "No species found to add, error. Sorensen_MinConf::add_species_min_conf_sorensen() " << std::endl; // something odd happened
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
// MSP end

void Sorensen_MinConf::set_fixed_species_sorensen()
{
    for (unsigned site = 0; site < n_sites; site++) {
        set_fixed_species_sorensen(site);
    }
}

void Sorensen_MinConf::set_fixed_species_sorensen(unsigned site)
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

std::vector<std::vector<int> > Sorensen_MinConf::gen_random_solution_sorensen()
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