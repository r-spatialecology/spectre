#ifndef TESTMINCONF_H
#define TESTMINCONF_H
#include <testthat.h>
#include "minconf.h"


class TestMinConf : public MinConf
{
public:
    using MinConf::MinConf;
    void test_add_species_min_conf(unsigned site, const std::vector<std::vector<int> > &target);
    std::vector<unsigned> test_calc_min_conflict_species(const unsigned site,
                                                    const std::vector<unsigned> free_species,
                                                    const std::vector<std::vector<int> > &target);
    std::vector<std::vector<int> > test_calculate_commonness();
    double test_calc_energy(const std::vector<std::vector<int> > &commonness,
                       const std::vector<std::vector<int> > &target);
    void test_gen_init_solution();
    std::vector<std::vector<int> > getTarget() const;
    std::vector<unsigned> test_calc_missing_species();
};

#endif // TESTMINCONF_H
