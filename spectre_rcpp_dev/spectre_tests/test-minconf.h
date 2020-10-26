#ifndef TESTMINCONF_H
#define TESTMINCONF_H
#include "../../src/minconf.h"


class TestMinConf : public MinConf
{
public:
    using MinConf::MinConf;
    void test_add_species_min_conf(unsigned site, const std::vector<std::vector<int> > &target);
    void test_remove_species_max_conf(unsigned site, const std::vector<std::vector<int> > &target);
    std::vector<unsigned> test_calc_min_conflict_species(const unsigned site,
                                                    const std::vector<unsigned> free_species,
                                                    const std::vector<std::vector<int> > &target);
    std::vector<unsigned> test_calc_max_conflict_species(const unsigned site,
                                                    const std::vector<unsigned> pesent_species,
                                                    const std::vector<std::vector<int> > &target);
    std::vector<std::vector<int> > test_calculate_commonness();
    double test_calc_energy(const std::vector<std::vector<int> > &commonness,
                       const std::vector<std::vector<int> > &target);
    void test_gen_init_solution();
};

#endif // TESTMINCONF_H
