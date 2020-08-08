#ifndef SORENSEN_CONSTRAINT_SATISFACTION_PROBLEM_H
#define SORENSEN_CONSTRAINT_SATISFACTION_PROBLEM_H
#include <vector>
#include <string>

class Sorensen_constraint_satisfaction_problem
{
public:
    Sorensen_constraint_satisfaction_problem(const std::vector<unsigned> &alpha_list,
                                    const unsigned gamma_div,
                                    const std::vector<double> &target, // MSP 
                                    const std::vector<int> &fixed_species_ = std::vector<int>(),
                                    const std::vector<int> &partial_solution = std::vector<int>(),
                                    const std::string norm = "sum");
    
    std::vector<std::vector<int> > solution;
    // MSP
    std::vector<std::vector<double> > solution_sorensen;
    // MSP
    unsigned p = 1;
    
protected:
    long max_steps = 0;
    
    std::vector<std::vector<int> > fixed_species;
    std::vector<std::vector<int> > fixed_species_idx;
    // MSP
    std::vector<std::vector<double> > target;
    // MSP 
    const std::vector<unsigned> alpha_list;
    std::vector<unsigned> tabu_list;
    const unsigned gamma_div;
    const unsigned n_sites;
    const std::string norm;
    
    // MSP start
    double calc_energy_sorensen(const std::vector<std::vector<double> > &commonness,
                                const std::vector<std::vector<double> > &target,
                                const bool use_custom_norm = false,
                                int omit_site = -1);
    double calc_energy_sum_sorensen(const std::vector<std::vector<double> > &commonness,
                                    const std::vector<std::vector<double> > &target,
                                    int omit_site = -1);
    double calc_energy_euclid_sorensen(const std::vector<std::vector<double> > &commonness,
                                       const std::vector<std::vector<double> > &target,
                                       int omit_site = -1);
    double calc_energy_p_sorensen(const std::vector<std::vector<double> > &commonness,
                                  const std::vector<std::vector<double> > &target,
                                  const int omit_site = -1);
    double calc_energy_max_sorensen(const std::vector<std::vector<double> > &commonness,
                                    const std::vector<std::vector<double> > &target,
                                    int omit_site = -1);
    
    // MSP end
    std::vector<unsigned> present_species_index(unsigned site, bool omit_fixed_species = true);
    std::vector<unsigned> present_species_index(unsigned site,
                                                const std::vector<std::vector<int> > fixed_species);
    std::vector<unsigned> absent_species_index(unsigned site);
    
    // MSP start
    std::vector<std::vector<double> > calculate_sorensen();
    std::vector<std::vector<double> > calculate_sorensen(const std::vector<std::vector<int> > &solution);
    // MSP end
    
    void update_solution_commonness_site(const std::vector<std::vector<int> > &solution_matrix, // remove??? 
                                         std::vector<std::vector<int> > &solution_commonness,
                                         const unsigned n_sites,
                                         const unsigned n_species,
                                         const unsigned site);
    
    // MSP start
    void update_solution_sorensen_site(const std::vector<std::vector<int> > &solution,
                                       std::vector<std::vector<int> > &commonness_result,
                                       std::vector<std::vector<double> > &sorensen_result,
                                       const unsigned n_sites,
                                       const unsigned n_species,
                                       const unsigned site);
    // MSP end
};

#endif // SORENSEN_CONSTRAINT_SATISFACTION_PROBLEM_H
