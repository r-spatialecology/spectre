#ifndef MINCONF_H
#define MINCONF_H
#include <random>
#include <vector>
#include <string>


class MinConf
{
public:
    MinConf(const std::vector<unsigned> &alpha_list,
            const unsigned gamma_div,
            const std::vector<int> &target_,
            const unsigned long seed,
            const std::vector<int> &fixed_species_ = std::vector<int>(),
            const std::vector<int> &partial_solution = std::vector<int>());

    int optimize(const long max_steps_,
                 bool verbose,
                 bool interruptible);
    unsigned calc_error_random_solution(const unsigned n = 10);

    std::vector<std::vector<int> > solution;
    std::vector<std::vector<int> > commonness;
    std::vector<int> iteration_count;
    std::vector<unsigned> error_vector;
    bool solution_has_best_error = true;
    const int RET_ABORT = -999;
    static const int NA = -2147483648; // == NA value in Rcpp

    static std::vector<std::vector<int> > calculate_commonness(const std::vector<std::vector<int> > &solution,
                                                               const unsigned n_sites);

protected:
    std::mt19937 rng;
    const double epsilon = 0.00001;

    std::vector<std::vector<int> > target;
    const std::vector<unsigned> alpha_list;
    std::vector<std::vector<int> > fixed_species;
    const unsigned gamma_div;
    const unsigned n_sites;

    void gen_init_solution(std::vector<unsigned> missing_species);
    std::vector<std::vector<int> > gen_random_solution();
    void set_fixed_species();
    void set_fixed_species(unsigned site);
    std::vector<unsigned> calc_missing_species();
    std::vector<unsigned> present_species_index(unsigned site, bool omit_fixed_species = true);
    std::vector<unsigned> present_species_index(unsigned site,
                                                const std::vector<std::vector<int> > partial_solution);
    std::vector<unsigned> absent_species_index(unsigned site);
    void add_species_min_conf(unsigned site,
                              const std::vector<std::vector<int> > &target);
    std::vector<unsigned> calc_min_conflict_species(const unsigned site,
                                                    const std::vector<unsigned> free_species,
                                                    const std::vector<std::vector<int> > &target);

    unsigned calc_error(const std::vector<std::vector<int> > &commonness,
                        const std::vector<std::vector<int> > &target);
    void update_solution_commonness();
};

#endif // MINCONF_H
